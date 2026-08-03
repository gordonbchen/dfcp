import os
from pathlib import Path

from plotly.offline.offline import get_plotlyjs


PLOTLY_ASSET = Path("docs/assets/plotly.min.js")


def ensure_plotly_asset(html_path: Path) -> str:
    asset_path = PLOTLY_ASSET
    asset_path.parent.mkdir(parents=True, exist_ok=True)
    plotly_js = get_plotlyjs() + "\n"
    if not asset_path.exists() or asset_path.read_text() != plotly_js:
        asset_path.write_text(plotly_js)
    return Path(os.path.relpath(asset_path.resolve(), html_path.parent.resolve())).as_posix()


def externalize_plotly_html(html_path: Path) -> None:
    asset_url = ensure_plotly_asset(html_path)
    html = html_path.read_text()
    marker = "<script>/**\n* plotly.js v"
    start = html.find(marker)
    if start < 0:
        return
    end = html.find("</script>", start)
    if end < 0:
        raise ValueError(f"unterminated inline Plotly script in {html_path}")
    replacement = f'<script src="{asset_url}"></script>'
    html_path.write_text(html[:start] + replacement + html[end + len("</script>"):])
