#!/usr/bin/env python3

import argparse
from array import array
import base64
from bisect import bisect_left, bisect_right
from dataclasses import dataclass
import gzip
import html
import json
from pathlib import Path
import struct
import subprocess
import sys
from typing import Iterable, Iterator, TextIO


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(SCRIPT_DIR.parent))

from plotly_html import ensure_plotly_asset  # noqa: E402
from window import make_windows  # noqa: E402


DATA_DIR = REPO_ROOT / "data/1000g_phase3_v5b"
DEFAULT_REF_VCF = DATA_DIR / "ref_target/ref.vcf.gz"
DEFAULT_OBSERVED = DATA_DIR / "dfcp_prep/observed_loci.txt"
DEFAULT_MAP = DATA_DIR / "plink.chr20.GRCh37.map"
DEFAULT_SEQ_FILE = DATA_DIR / "dfcp_prep/ref.bin"


@dataclass(frozen=True)
class GeneticMap:
    positions: list[int]
    cm: list[float]
    format_name: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build an interactive report for choosing overlapping DFCP chr20 windows.",
    )
    parser.add_argument("--ref-vcf", type=Path, default=DEFAULT_REF_VCF)
    parser.add_argument("--observed-loci", type=Path, default=DEFAULT_OBSERVED)
    parser.add_argument("--genetic-map", type=Path, default=DEFAULT_MAP)
    parser.add_argument("--seq-file", type=Path, default=DEFAULT_SEQ_FILE)
    parser.add_argument("--physical-only", action="store_true")
    parser.add_argument("--density-bins", type=int, default=240)
    parser.add_argument("--output", type=Path, default=Path("window.html"))
    args = parser.parse_args()
    if args.density_bins < 20:
        parser.error("--density-bins must be at least 20")
    return args


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def read_variant_positions(vcf: Path) -> tuple[str, list[int]]:
    if not vcf.is_file():
        raise FileNotFoundError(f"reference VCF not found: {vcf}")

    command = ["bcftools", "query", "-f", "%CHROM\\t%POS\\n", str(vcf)]
    try:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except FileNotFoundError as error:
        raise RuntimeError("bcftools is required to read --ref-vcf") from error

    assert process.stdout is not None
    chromosome = ""
    positions: list[int] = []
    for row, line in enumerate(process.stdout, start=1):
        fields = line.split()
        if len(fields) != 2:
            process.kill()
            raise ValueError(f"invalid bcftools output on row {row}: {line.rstrip()}")
        if not chromosome:
            chromosome = fields[0]
        elif fields[0] != chromosome:
            process.kill()
            raise ValueError(f"reference VCF contains chromosomes {chromosome} and {fields[0]}")
        position = int(fields[1])
        if positions and position < positions[-1]:
            process.kill()
            raise ValueError(f"reference VCF positions are reordered at row {row}")
        positions.append(position)

    assert process.stderr is not None
    stderr = process.stderr.read()
    return_code = process.wait()
    if return_code:
        raise RuntimeError(f"bcftools query failed ({return_code}): {stderr.strip()}")
    if not positions:
        raise ValueError(f"reference VCF contains no variants: {vcf}")
    return chromosome, positions


def read_obs_ls(path: Path, n_loci: int) -> list[int]:
    if not path.is_file():
        raise FileNotFoundError(f"observed-loci file not found: {path}")
    loci: list[int] = []
    with path.open() as stream:
        for row, line in enumerate(stream, start=1):
            stripped = line.strip()
            if not stripped:
                continue
            try:
                locus = int(stripped)
            except ValueError as error:
                raise ValueError(f"invalid observed locus on row {row}: {stripped}") from error
            if not 0 <= locus < n_loci:
                raise ValueError(f"observed locus {locus} on row {row} is outside [0, {n_loci})")
            if loci and locus <= loci[-1]:
                raise ValueError("observed loci must be unique and strictly increasing")
            loci.append(locus)
    if not loci:
        raise ValueError(f"observed-loci file is empty: {path}")
    return loci


def normalized_chromosome(chromosome: str) -> str:
    return chromosome.removeprefix("chr").upper()


def parse_map_row(fields: list[str], chromosome: str) -> tuple[int, float, str] | None:
    wanted = normalized_chromosome(chromosome)
    if len(fields) >= 4:
        if normalized_chromosome(fields[0]) != wanted:
            return None
        try:
            return int(fields[3]), float(fields[2]), "PLINK (chromosome, marker, cM, bp)"
        except ValueError:
            return None
    if len(fields) >= 3:
        try:
            return int(fields[0]), float(fields[2]), "SHAPEIT/IMPUTE (bp, cM/Mb, cM)"
        except ValueError:
            return None
    if len(fields) == 2:
        try:
            return int(fields[0]), float(fields[1]), "two-column (bp, cM)"
        except ValueError:
            return None
    return None


def read_genetic_map(path: Path, chromosome: str) -> GeneticMap:
    positions: list[int] = []
    cm: list[float] = []
    formats: set[str] = set()
    with open_text(path) as stream:
        for line in stream:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parsed = parse_map_row(stripped.split(), chromosome)
            if parsed is None:
                continue
            position, genetic_position, format_name = parsed
            positions.append(position)
            cm.append(genetic_position)
            formats.add(format_name)

    if len(positions) < 2:
        raise ValueError(f"genetic map has fewer than two rows for chromosome {chromosome}: {path}")
    if len(formats) != 1:
        raise ValueError(f"genetic map mixes formats: {path}")
    for index in range(1, len(positions)):
        if positions[index] <= positions[index - 1]:
            raise ValueError(f"genetic-map bp positions are not strictly increasing at row {index + 1}")
        if cm[index] < cm[index - 1]:
            raise ValueError(f"genetic-map cM positions decrease at row {index + 1}")
    return GeneticMap(positions, cm, formats.pop())


def interpolate_cm(position: int, genetic_map: GeneticMap, hint: int = 0) -> tuple[float, int]:
    map_positions = genetic_map.positions
    map_cm = genetic_map.cm
    last_segment = len(map_positions) - 2
    segment = min(max(hint, 0), last_segment)
    while segment < last_segment and position > map_positions[segment + 1]:
        segment += 1
    while segment > 0 and position < map_positions[segment]:
        segment -= 1
    left_bp = map_positions[segment]
    right_bp = map_positions[segment + 1]
    fraction = (position - left_bp) / (right_bp - left_bp)
    return map_cm[segment] + fraction * (map_cm[segment + 1] - map_cm[segment]), segment


def coordinate_stream(positions: Iterable[int], genetic_map: GeneticMap | None) -> Iterator[float]:
    if genetic_map is None:
        yield from (position / 1_000_000 for position in positions)
        return
    hint = 0
    for position in positions:
        coordinate, hint = interpolate_cm(position, genetic_map, hint)
        yield coordinate


def density_data(
    positions: list[int],
    obs_ls: list[int],
    genetic_map: GeneticMap | None,
    n_bins: int,
) -> dict:
    start = next(coordinate_stream(positions[:1], genetic_map))
    end = next(coordinate_stream(positions[-1:], genetic_map))
    span = end - start
    if span <= 0:
        raise ValueError("variant coordinate span must be positive")
    width = span / n_bins
    reference_counts = [0] * n_bins
    obs_counts = [0] * n_bins
    obs_cursor = 0
    for locus, coordinate in enumerate(coordinate_stream(positions, genetic_map)):
        bin_index = min(int((coordinate - start) / width), n_bins - 1)
        reference_counts[bin_index] += 1
        if obs_cursor < len(obs_ls) and locus == obs_ls[obs_cursor]:
            obs_counts[bin_index] += 1
            obs_cursor += 1
    centers = [start + (index + 0.5) * width for index in range(n_bins)]
    return {
        "x": centers,
        "reference_density": [count / width for count in reference_counts],
        "observed_density": [count / width for count in obs_counts],
        "bin_width": width,
        "axis_unit": "cM" if genetic_map else "Mb",
    }


def packed_base64(typecode: str, values: Iterable[int | float]) -> str:
    packed = array(typecode, values)
    if sys.byteorder != "little":
        packed.byteswap()
    return base64.b64encode(packed.tobytes()).decode()


def read_n_sequences(path: Path) -> int | None:
    if not path.is_file():
        return None
    with path.open("rb") as stream:
        header = stream.read(12)
    if len(header) != 12:
        raise ValueError(f"DFCP sequence header is truncated: {path}")
    magic, n_sequences, _ = struct.unpack("<4sII", header)
    if magic != b"DFCP":
        raise ValueError(f"invalid DFCP sequence magic in {path}")
    return n_sequences


CSS = r"""
:root { color-scheme:light; --ink:#2a3f5f; --muted:#697386; --line:#e5e7eb;
  --blue:#3366cc; --orange:#ef8a35; }
* { box-sizing:border-box; }
body { margin:0; color:var(--ink); background:white;
  font-family:Inter,ui-sans-serif,system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; }
.window-sidebar { position:fixed; inset:0 auto 0 0; z-index:5; width:320px; overflow-y:auto;
  padding:24px 20px; border-right:1px solid var(--line); background:#fafbfc; }
.window-sidebar h1 { margin:0 0 7px; font-size:20px; font-weight:650; }
.source { margin:0; color:var(--muted); font:11px/1.5 ui-monospace,monospace; overflow-wrap:anywhere; }
.sidebar-section { margin-top:22px; }
.sidebar-section h2 { margin:0 0 9px; font-size:14px; font-weight:650; }
.control { margin:0 0 16px; }
.control-head { display:flex; align-items:center; justify-content:space-between; gap:10px;
  margin-bottom:5px; font-size:12px; }
.control-head input { width:92px; padding:4px 6px; border:1px solid #b8c0cc; border-radius:4px;
  color:var(--ink); background:white; font:600 12px ui-monospace,monospace; text-align:right; }
.control input[type=range] { width:100%; accent-color:var(--blue); }
.note { margin:6px 0 0; color:var(--muted); font-size:11px; line-height:1.45; }
.metric-table { width:100%; border-collapse:collapse; font-size:11px; }
.metric-table th,.metric-table td { padding:5px 2px; border-bottom:1px solid var(--line); }
.metric-table th { width:58%; text-align:left; color:var(--muted); font-weight:500; }
.metric-table td { text-align:right; font:600 11px ui-monospace,monospace; }
.report-main { width:calc(100% - 320px); margin-left:320px; padding:14px 28px 60px; }
.intro { margin:0 8px 14px; }
.intro h2 { margin:0 0 5px; font-size:20px; }
.intro p { max-width:940px; margin:0; color:var(--muted); font-size:12px; line-height:1.55; }
.grid { display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); gap:30px 26px; }
.panel { min-width:0; }
.wide { grid-column:1/-1; }
@media (max-width:1050px) {
  .window-sidebar { width:280px; }
  .report-main { width:calc(100% - 280px); margin-left:280px; padding-inline:18px; }
  .grid { grid-template-columns:1fr; }
  .wide { grid-column:auto; }
}
"""


JAVASCRIPT = r"""
(() => {
  const data = JSON.parse(document.getElementById('report-data').textContent);
  const config = {responsive: true, displaylogo: false};
  const colors = {blue: '#3366cc', orange: '#ef8a35', red: '#d1495b', green: '#2a9d8f'};

  function decodeBase64(encoded, Type) {
    if (!encoded) return new Type(0);
    const binary = atob(encoded);
    const bytes = new Uint8Array(binary.length);
    for (let index = 0; index < binary.length; ++index) bytes[index] = binary.charCodeAt(index);
    return new Type(bytes.buffer);
  }

  const positions = decodeBase64(data.positions_b64, Uint32Array);
  const observed = decodeBase64(data.observed_b64, Uint32Array);
  const mapPositions = decodeBase64(data.map_positions_b64, Uint32Array);
  const mapCm = decodeBase64(data.map_cm_b64, Float64Array);
  const defaultStarts = decodeBase64(data.default_starts_b64, Uint32Array);
  const defaultEnds = decodeBase64(data.default_ends_b64, Uint32Array);
  const hasMap = mapPositions.length > 1;

  function lowerBound(values, target) {
    let left = 0;
    let right = values.length;
    while (left < right) {
      const middle = (left + right) >>> 1;
      if (values[middle] < target) left = middle + 1;
      else right = middle;
    }
    return left;
  }

  function geneticPosition(bp) {
    if (!hasMap) return bp / 1e6;
    let segment;
    if (bp <= mapPositions[0]) segment = 0;
    else if (bp >= mapPositions[mapPositions.length - 1]) segment = mapPositions.length - 2;
    else segment = lowerBound(mapPositions, bp) - 1;
    const leftBp = mapPositions[segment];
    const rightBp = mapPositions[segment + 1];
    const fraction = (bp - leftBp) / (rightBp - leftBp);
    return mapCm[segment] + fraction * (mapCm[segment + 1] - mapCm[segment]);
  }

  function coordinateAt(locus) {
    return geneticPosition(positions[locus]);
  }

  function makeWindows(nWindows, overlap) {
    const memberships = positions.length + (nWindows - 1) * overlap;
    const baseSize = Math.floor(memberships / nWindows);
    const remainder = memberships % nWindows;
    if (baseSize <= overlap && nWindows > 1) throw new Error('overlap leaves no new loci per window');
    const windows = [];
    let start = 0;
    for (let windowIndex = 0; windowIndex < nWindows; ++windowIndex) {
      const size = baseSize + (windowIndex < remainder ? 1 : 0);
      const end = start + size;
      windows.push({windowIndex, start, end, size});
      start = end - overlap;
    }
    if (windows[windows.length - 1].end !== positions.length) {
      throw new Error('internal window construction error');
    }
    if (nWindows === data.default_n_windows && overlap === data.default_overlap) {
      for (let index = 0; index < windows.length; ++index) {
        if (windows[index].start !== defaultStarts[index] || windows[index].end !== defaultEnds[index]) {
          throw new Error('visualization and file-writer window calculations disagree');
        }
      }
    }
    return windows;
  }

  function stats(values) {
    if (!values.length) return null;
    let minimum = values[0];
    let maximum = values[0];
    let total = 0;
    for (const value of values) {
      minimum = Math.min(minimum, value);
      maximum = Math.max(maximum, value);
      total += value;
    }
    const mean = total / values.length;
    let variance = 0;
    for (const value of values) variance += (value - mean) ** 2;
    return {minimum, mean, maximum, cv: mean ? Math.sqrt(variance / values.length) / mean : 0};
  }

  function integer(value) {
    return Math.round(value).toLocaleString('en-US');
  }

  function decimal(value, digits = 2) {
    return Number(value).toLocaleString('en-US', {
      minimumFractionDigits: digits,
      maximumFractionDigits: digits,
    });
  }

  function rangeText(summary, digits = 1) {
    if (!summary) return 'n/a';
    return `${decimal(summary.minimum, digits)} / ${decimal(summary.mean, digits)} / ` +
      `${decimal(summary.maximum, digits)}`;
  }

  function integerRange(summary) {
    if (!summary) return 'n/a';
    return `${integer(summary.minimum)} / ${decimal(summary.mean, 1)} / ${integer(summary.maximum)}`;
  }

  function setText(id, value) {
    document.getElementById(id).textContent = value;
  }

  function commonLayout(title, xTitle, yTitle) {
    return {
      title: {text: title, x: 0.03, font: {size: 16}},
      margin: {l: 72, r: 30, t: 54, b: 58},
      height: 390,
      paper_bgcolor: 'white',
      plot_bgcolor: 'white',
      font: {color: '#2a3f5f', size: 11},
      xaxis: {title: {text: xTitle}, gridcolor: '#edf0f3', zeroline: false, automargin: true},
      yaxis: {title: {text: yTitle}, gridcolor: '#edf0f3', zeroline: false, automargin: true},
      legend: {orientation: 'h', y: 1.08, x: 1, xanchor: 'right'},
      hovermode: 'closest',
    };
  }

  function render(nWindows, overlap) {
    const windows = makeWindows(nWindows, overlap);
    const sizes = [];
    const observedCounts = [];
    const physicalSpans = [];
    const coordinateSpans = [];
    const overlapObserved = [];
    const overlapSpans = [];
    const windowNumbers = [];

    for (const window of windows) {
      sizes.push(window.size);
      observedCounts.push(lowerBound(observed, window.end) - lowerBound(observed, window.start));
      physicalSpans.push((positions[window.end - 1] - positions[window.start]) / 1e6);
      coordinateSpans.push(coordinateAt(window.end - 1) - coordinateAt(window.start));
      windowNumbers.push(window.windowIndex + 1);
      if (window.windowIndex) {
        const overlapStart = window.start;
        const overlapEnd = windows[window.windowIndex - 1].end;
        overlapObserved.push(lowerBound(observed, overlapEnd) - lowerBound(observed, overlapStart));
        overlapSpans.push(
          overlap ? coordinateAt(overlapEnd - 1) - coordinateAt(overlapStart) : 0,
        );
      }
    }

    const sizeStats = stats(sizes);
    const observedStats = stats(observedCounts);
    const physicalStats = stats(physicalSpans);
    const coordinateStats = stats(coordinateSpans);
    const overlapObservedStats = stats(overlapObserved);
    const overlapSpanStats = stats(overlapSpans);
    const membershipFactor = sizes.reduce((total, value) => total + value, 0) / positions.length;

    setText('metric-window-loci', integerRange(sizeStats));
    setText('metric-observed', integerRange(observedStats));
    setText('metric-physical-span', rangeText(physicalStats, 2));
    setText('metric-coordinate-span', rangeText(coordinateStats, 2));
    setText('metric-coordinate-cv', `${decimal(100 * coordinateStats.cv, 1)}%`);
    setText('metric-overlap-observed', integerRange(overlapObservedStats));
    setText('metric-overlap-span', rangeText(overlapSpanStats, 3));
    setText('metric-membership', `${decimal(membershipFactor, 3)}× (+${decimal(
      100 * (membershipFactor - 1), 1,
    )}%)`);

    if (data.n_sequences) {
      const maximumLoci = sizeStats.maximum;
      const packedMib = data.n_sequences * Math.ceil(maximumLoci / 64) * 8 / 2 ** 20;
      const assignmentsGib = data.n_sequences * (2 * maximumLoci - 1) * 8 / 2 ** 30;
      setText('metric-packed-memory', `${decimal(packedMib, 1)} MiB`);
      setText('metric-assignment-memory', `${decimal(assignmentsGib, 2)} GiB`);
    }

    const densityLayout = commonLayout(
      `Chr${data.chromosome}: loci density and current windows`,
      hasMap ? `Genetic position on chr${data.chromosome} (cM)` :
        `Physical position on chr${data.chromosome} (Mb)`,
      `Observed loci per ${data.density.axis_unit}`,
    );
    densityLayout.height = 430;
    densityLayout.margin.r = 78;
    densityLayout.yaxis2 = {
      title: {text: `Reference loci per ${data.density.axis_unit}`},
      overlaying: 'y',
      side: 'right',
      showgrid: false,
      automargin: true,
    };
    densityLayout.shapes = windows.map((window, index) => ({
      type: 'rect',
      xref: 'x',
      yref: 'paper',
      x0: coordinateAt(window.start),
      x1: coordinateAt(window.end - 1),
      y0: 0,
      y1: 1,
      fillcolor: index % 2 ? '#3366cc' : '#ef8a35',
      opacity: 0.055,
      line: {width: 0},
      layer: 'below',
    }));
    Plotly.react('density-plot', [
      {
        type: 'scatter',
        mode: 'lines',
        x: data.density.x,
        y: data.density.observed_density,
        name: 'Observed density',
        line: {color: colors.blue, width: 1.5},
        hovertemplate: `%{x:.3f} ${data.density.axis_unit}<br>%{y:.1f} observed / ` +
          `${data.density.axis_unit}<extra></extra>`,
      },
      {
        type: 'scatter',
        mode: 'lines',
        x: data.density.x,
        y: data.density.reference_density,
        yaxis: 'y2',
        name: 'Reference density',
        line: {color: colors.orange, width: 1.5},
        hovertemplate: `%{x:.3f} ${data.density.axis_unit}<br>%{y:.1f} reference / ` +
          `${data.density.axis_unit}<extra></extra>`,
      },
    ], densityLayout, config);

    const observedLayout = commonLayout(
      'Observed loci per window',
      'Observed loci in a window',
      'Windows',
    );
    Plotly.react('observed-histogram', [{
      type: 'histogram',
      x: observedCounts,
      marker: {color: colors.blue},
      nbinsx: Math.max(6, Math.ceil(Math.sqrt(nWindows))),
      hovertemplate: '%{x}<br>%{y} windows<extra></extra>',
    }], observedLayout, config);

    const spanName = hasMap ? 'Genetic span per window' : 'Physical span per window';
    const spanUnit = hasMap ? 'cM' : 'Mb';
    const spanValues = hasMap ? coordinateSpans : physicalSpans;
    const spanLayout = commonLayout(spanName, `Window span (${spanUnit})`, 'Windows');
    Plotly.react('span-histogram', [{
      type: 'histogram',
      x: spanValues,
      marker: {color: colors.orange},
      nbinsx: Math.max(6, Math.ceil(Math.sqrt(nWindows))),
      hovertemplate: `%{x:.4f} ${spanUnit}<br>%{y} windows<extra></extra>`,
    }], spanLayout, config);

    const overlapLayout = commonLayout(
      'Observed loci shared by neighboring windows',
      'Observed loci in the shared region',
      'Neighboring window pairs',
    );
    const overlapTrace = overlapObserved.length ? [{
      type: 'histogram',
      x: overlapObserved,
      marker: {color: colors.green},
      nbinsx: Math.max(6, Math.ceil(Math.sqrt(overlapObserved.length))),
      hovertemplate: '%{x}<br>%{y} window pairs<extra></extra>',
    }] : [];
    Plotly.react('overlap-histogram', overlapTrace, overlapLayout, config);

    const profileLayout = commonLayout(
      'Window-by-window balance',
      'Window number',
      'Observed loci in the window',
    );
    profileLayout.height = 430;
    profileLayout.margin.r = 78;
    profileLayout.yaxis2 = {
      title: {text: hasMap ? 'Genetic span (cM)' : 'Physical span (Mb)'},
      overlaying: 'y',
      side: 'right',
      showgrid: false,
      automargin: true,
    };
    Plotly.react('window-profile', [
      {
        type: 'scatter',
        mode: 'lines+markers',
        x: windowNumbers,
        y: observedCounts,
        name: 'Observed loci',
        line: {color: colors.blue},
        marker: {size: 5},
        hovertemplate: 'window %{x}<br>%{y} observed<extra></extra>',
      },
      {
        type: 'scatter',
        mode: 'lines+markers',
        x: windowNumbers,
        y: spanValues,
        yaxis: 'y2',
        name: `${spanUnit} span`,
        line: {color: colors.orange},
        marker: {size: 5},
        hovertemplate: `window %{x}<br>%{y:.4f} ${spanUnit}<extra></extra>`,
      },
    ], profileLayout, config);
  }

  const nRange = document.getElementById('n-windows-range');
  const nNumber = document.getElementById('n-windows-number');
  const overlapRange = document.getElementById('overlap-range');
  const overlapNumber = document.getElementById('overlap-number');
  let frame = null;

  function scheduleRender() {
    if (frame !== null) cancelAnimationFrame(frame);
    frame = requestAnimationFrame(() => {
      frame = null;
      const nWindows = Math.max(1, Math.min(Math.round(Number(nNumber.value)), Number(nNumber.max)));
      const overlap = Math.max(
        0,
        Math.min(Math.round(Number(overlapNumber.value)), Number(overlapNumber.max)),
      );
      nRange.value = nWindows;
      nNumber.value = nWindows;
      overlapRange.value = overlap;
      overlapNumber.value = overlap;
      render(nWindows, overlap);
    });
  }

  function connect(range, number) {
    range.addEventListener('input', () => {
      number.value = range.value;
      scheduleRender();
    });
    number.addEventListener('input', () => {
      range.value = number.value;
      scheduleRender();
    });
  }

  connect(nRange, nNumber);
  connect(overlapRange, overlapNumber);
  scheduleRender();
})();
"""


def static_metric_row(label: str, value: str) -> str:
    return f"<tr><th>{html.escape(label)}</th><td>{html.escape(value)}</td></tr>"


def write_report(
    args: argparse.Namespace,
    chromosome: str,
    positions: list[int],
    obs_ls: list[int],
    genetic_map: GeneticMap | None,
    n_sequences: int | None,
) -> None:
    args.output.parent.mkdir(parents=True, exist_ok=True)
    density = density_data(positions, obs_ls, genetic_map, args.density_bins)
    physical_span_mb = (positions[-1] - positions[0]) / 1_000_000
    if genetic_map:
        start_cm, _ = interpolate_cm(positions[0], genetic_map)
        end_cm, _ = interpolate_cm(positions[-1], genetic_map, len(genetic_map.positions) - 2)
        genetic_span_cm = end_cm - start_cm
        outside_map = bisect_left(positions, genetic_map.positions[0])
        outside_map += len(positions) - bisect_right(positions, genetic_map.positions[-1])
    else:
        genetic_span_cm = None
        outside_map = 0

    maximum_windows = min(500, len(positions))
    default_n_windows = min(100, maximum_windows)
    maximum_overlap = min(10_000, len(positions) - 1)
    default_overlap = min(500, maximum_overlap)
    default_windows = make_windows(len(positions), default_n_windows, default_overlap)
    report_data = {
        "chromosome": chromosome,
        "positions_b64": packed_base64("I", positions),
        "observed_b64": packed_base64("I", obs_ls),
        "map_positions_b64": packed_base64("I", genetic_map.positions) if genetic_map else "",
        "map_cm_b64": packed_base64("d", genetic_map.cm) if genetic_map else "",
        "density": density,
        "n_sequences": n_sequences,
        "default_starts_b64": packed_base64("I", (window.start for window in default_windows)),
        "default_ends_b64": packed_base64("I", (window.end for window in default_windows)),
        "default_n_windows": default_n_windows,
        "default_overlap": default_overlap,
    }
    data_json = json.dumps(report_data, separators=(",", ":")).replace("</", "<\\/")

    static_rows = [
        static_metric_row("Reference loci", f"{len(positions):,}"),
        static_metric_row("Observed loci", f"{len(obs_ls):,}"),
        static_metric_row("Observed fraction", f"{100 * len(obs_ls) / len(positions):.2f}%"),
        static_metric_row("Covered physical span", f"{physical_span_mb:.2f} Mb"),
    ]
    if genetic_span_cm is not None:
        static_rows.extend([
            static_metric_row("Genetic span", f"{genetic_span_cm:.2f} cM"),
            static_metric_row("Reference loci / cM", f"{len(positions) / genetic_span_cm:,.0f}"),
            static_metric_row("Observed loci / cM", f"{len(obs_ls) / genetic_span_cm:,.1f}"),
            static_metric_row("Map-extrapolated loci", f"{outside_map:,}"),
        ])
    if n_sequences is not None:
        static_rows.append(static_metric_row("Reference sequences", f"{n_sequences:,}"))

    map_description = "Physical coordinates only"
    if genetic_map:
        map_description = f"{args.genetic_map.name}; {genetic_map.format_name}"
    plotly_asset = ensure_plotly_asset(args.output)
    page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>DFCP window selection</title>
<script src="{html.escape(plotly_asset, quote=True)}"></script><style>{CSS}</style></head>
<body><aside class="window-sidebar"><h1>DFCP windows</h1>
<p class="source">{html.escape(str(args.ref_vcf))}</p>
<section class="sidebar-section"><h2>Window controls</h2>
<div class="control"><div class="control-head"><label for="n-windows-range">Windows</label>
<input id="n-windows-number" type="number" min="1" max="{maximum_windows}"
value="{default_n_windows}">
</div><input id="n-windows-range" type="range" min="1" max="{maximum_windows}"
value="{default_n_windows}">
</div><div class="control"><div class="control-head"><label for="overlap-range">Overlap loci</label>
<input id="overlap-number" type="number" min="0" max="{maximum_overlap}" value="{default_overlap}">
</div><input id="overlap-range" type="range" min="0" max="{maximum_overlap}" step="250"
value="{default_overlap}"><p class="note">The number box permits exact values. Neighboring windows share
exactly this many reference loci.</p></div></section>
<section class="sidebar-section"><h2>Current windows</h2><table class="metric-table">
<tr><th>Window loci min / mean / max</th><td id="metric-window-loci"></td></tr>
<tr><th>Observed loci min / mean / max</th><td id="metric-observed"></td></tr>
<tr><th>Physical span (Mb) min / mean / max</th><td id="metric-physical-span"></td></tr>
<tr><th>{'Genetic span (cM)' if genetic_map else 'Physical span (Mb)'} min / mean / max</th>
<td id="metric-coordinate-span"></td></tr>
<tr><th>{'cM' if genetic_map else 'Mb'} span CV</th><td id="metric-coordinate-cv"></td></tr>
<tr><th>Observed loci shared by neighbors</th><td id="metric-overlap-observed"></td></tr>
<tr><th>Shared-region {'cM' if genetic_map else 'Mb'} min / mean / max</th>
<td id="metric-overlap-span"></td></tr>
<tr><th>Total locus work</th><td id="metric-membership"></td></tr>
<tr><th>Packed reference bits / largest window</th><td id="metric-packed-memory">n/a</td></tr>
<tr><th>R/Q assignment pointers / largest window</th><td id="metric-assignment-memory">n/a</td></tr>
</table><p class="note">Ranges are min / mean / max. Physical span is the bp distance between a window's
first and last loci; equal locus counts can span different Mb because variant density varies. Shared loci are
observed target markers present in both neighboring windows. Memory rows estimate only packed sequence
bits and the two assignment-pointer arrays, excluding clusters, messages, parameters, and allocator
overhead.</p></section>
<section class="sidebar-section"><h2>Chromosome</h2><table class="metric-table">
{''.join(static_rows)}</table><p class="note">Covered physical span is the bp distance from the first
reference locus to the last, not the full chromosome assembly length.</p></section>
<section class="sidebar-section"><h2>Genetic map</h2>
<p class="note">{html.escape(map_description)}. Map-extrapolated loci fall just outside the map's first or
last bp marker and use the nearest segment's recombination rate.</p></section></aside>
<main class="report-main"><header class="intro"><h2>Equal-locus overlapping windows</h2>
<p>For exact overlap O, L = sum(window sizes) − (n − 1)O. Integer remainders are distributed so
window sizes differ by at most one locus. Counts include overlap memberships; the chromosome itself
is covered exactly once after duplicate overlaps are removed.</p></header><div class="grid">
<section class="panel wide"><div id="density-plot"></div></section>
<section class="panel"><div id="observed-histogram"></div></section>
<section class="panel"><div id="span-histogram"></div></section>
<section class="panel wide"><div id="overlap-histogram"></div></section>
<section class="panel wide"><div id="window-profile"></div></section></div></main>
<script id="report-data" type="application/json">{data_json}</script><script>{JAVASCRIPT}</script>
</body></html>
"""
    args.output.write_text(page)


def main() -> None:
    args = parse_args()
    chromosome, positions = read_variant_positions(args.ref_vcf)
    obs_ls = read_obs_ls(args.observed_loci, len(positions))
    genetic_map = None
    if not args.physical_only:
        if args.genetic_map.is_file():
            genetic_map = read_genetic_map(args.genetic_map, chromosome)
        else:
            print(f"warning: genetic map not found; using Mb only: {args.genetic_map}", file=sys.stderr)
    n_sequences = read_n_sequences(args.seq_file)
    write_report(args, chromosome, positions, obs_ls, genetic_map, n_sequences)
    map_status = f"{len(genetic_map.positions):,} map rows" if genetic_map else "physical coordinates"
    print(
        f"wrote {args.output} from {len(positions):,} variants, {len(obs_ls):,} observed loci, "
        f"and {map_status}",
    )


if __name__ == "__main__":
    main()
