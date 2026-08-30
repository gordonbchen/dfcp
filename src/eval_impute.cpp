#include <algorithm>
#include <array>
#include <charconv>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <format>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>
#include "impute_io.hpp"
#include "io.hpp"
#include "seq_array.hpp"


namespace {

namespace fs = std::filesystem;

constexpr std::string_view manifest_header =
    "window\tstart\tend\tplanned_loci\tchrom\tfirst_pos\tlast_pos\tgenerated"
    "\tref_loci\tobserved\tmasked\toverlap_previous";
constexpr std::uint16_t fixed_point_max = std::numeric_limits<std::uint16_t>::max();

struct Window {
    int index;
    int start;
    int end;
    fs::path dir;
};

struct RefCount {
    int ac;
    int an;
};

struct MaskedLocus {
    int column;
    int mac;
    bool alt_is_minor;
};

struct MacStats {
    std::uint64_t n_loci = 0;
    std::uint64_t n_predictions = 0;
    std::uint64_t n_correct = 0;
    long double mean_q = 0.0;
    long double mean_y = 0.0;
    long double m2_q = 0.0;
    long double m2_y = 0.0;
    long double covariance = 0.0;

    void add(std::uint16_t q, int y) {
        ++n_predictions;
        long double dn = n_predictions;
        long double dq = q - mean_q;
        long double dy = y - mean_y;
        mean_q += dq / dn;
        mean_y += dy / dn;
        m2_q += dq * (q - mean_q);
        m2_y += dy * (y - mean_y);
        covariance += dq * (y - mean_y);
        n_correct += (static_cast<int>(q >= (std::uint16_t{1} << 15)) == y);
    }

    double r2() const {
        if (m2_q <= 0.0 || m2_y <= 0.0) { return -1.0; }
        long double value = covariance * covariance / (m2_q * m2_y);
        return static_cast<double>(std::clamp(value, 0.0L, 1.0L));
    }

    double accuracy() const {
        return static_cast<double>(n_correct) / n_predictions;
    }
};

std::vector<std::string_view> split_tabs(const std::string& line) {
    std::vector<std::string_view> fields;
    std::string_view remaining = line;
    while (true) {
        std::size_t tab = remaining.find('\t');
        fields.push_back(remaining.substr(0, tab));
        if (tab == std::string_view::npos) { return fields; }
        remaining.remove_prefix(tab + 1);
    }
}

int parse_int(std::string_view text, std::string_view field) {
    int value;
    auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), value);
    if (error != std::errc{} || end != text.data() + text.size()) {
        throw std::runtime_error(std::format("Invalid {}: {}", field, text));
    }
    return value;
}

std::string shell_quote(const fs::path& path) {
    std::string quoted{"'"};
    for (char c : path.string()) {
        quoted += c == '\'' ? "'\\''" : std::string(1, c);
    }
    return quoted + "'";
}

std::string query_ref_counts(const fs::path& vcf) {
    std::string command = "bcftools query -f '%AC\\t%AN\\n' " + shell_quote(vcf);
    FILE* pipe = popen(command.c_str(), "r");
    if (pipe == nullptr) {
        throw std::runtime_error("Failed to run bcftools query.");
    }

    std::string output;
    std::array<char, 4096> buffer{};
    while (fgets(buffer.data(), buffer.size(), pipe) != nullptr) {
        output += buffer.data();
    }
    bool read_failed = ferror(pipe);
    int status = pclose(pipe);
    if (read_failed || status != 0) {
        throw std::runtime_error(std::format("Failed to query reference allele counts: {}", vcf.string()));
    }
    return output;
}

std::vector<RefCount> read_ref_counts(const fs::path& vcf) {
    std::istringstream input(query_ref_counts(vcf));
    std::vector<RefCount> counts;
    std::string line;
    while (std::getline(input, line)) {
        std::vector<std::string_view> fields = split_tabs(line);
        if (fields.size() != 2 || fields[0].find(',') != std::string_view::npos) {
            throw std::runtime_error(
                std::format("Invalid biallelic AC/AN row in {}: {}", vcf.string(), line)
            );
        }
        int ac = parse_int(fields[0], "AC");
        int an = parse_int(fields[1], "AN");
        if (ac < 0 || an <= 0 || ac > an) {
            throw std::runtime_error(std::format("Invalid AC/AN in {}: {}/{}", vcf.string(), ac, an));
        }
        counts.push_back({ac, an});
    }
    return counts;
}

std::vector<Window> read_windows(const fs::path& root) {
    std::ifstream input(root / "windows.tsv");
    if (!input.is_open()) {
        throw std::runtime_error("Failed to open windows.tsv.");
    }

    std::string line;
    std::getline(input, line);
    if (line != manifest_header) {
        throw std::runtime_error("Unexpected windows.tsv header.");
    }

    std::vector<Window> windows;
    int row = 1;
    while (std::getline(input, line)) {
        ++row;
        std::vector<std::string_view> fields = split_tabs(line);
        if (fields.size() != 12) {
            throw std::runtime_error(std::format("Invalid windows.tsv row {}.", row));
        }
        int index = parse_int(fields[0], "window index");
        if (fields[7] != "0" && fields[7] != "1") {
            throw std::runtime_error(std::format("Invalid generated value on windows.tsv row {}.", row));
        }
        fs::path dir = root / std::format("window_{:04d}", index);
        if (fields[7] == "0" || !fs::exists(dir / "probs.bin")) { continue; }

        int start = parse_int(fields[1], "window start");
        int end = parse_int(fields[2], "window end");
        int planned_loci = parse_int(fields[3], "planned loci");
        int ref_loci = parse_int(fields[8], "reference loci");
        if (start < 0 || end <= start || planned_loci != end - start || ref_loci != planned_loci) {
            throw std::runtime_error(std::format("Inconsistent dimensions on windows.tsv row {}.", row));
        }
        windows.push_back({index, start, end, std::move(dir)});
    }
    if (!input.eof()) {
        throw std::runtime_error("Failed while reading windows.tsv.");
    }
    if (windows.empty()) {
        throw std::runtime_error("No generated windows contain probs.bin.");
    }
    std::ranges::sort(windows, {}, &Window::index);
    for (std::size_t i = 1; i < windows.size(); ++i) {
        if (windows[i-1].index == windows[i].index) {
            throw std::runtime_error("windows.tsv contains a duplicate window index.");
        }
    }
    return windows;
}

std::vector<int> assign_loci(const std::vector<Window>& windows) {
    int max_end = std::ranges::max_element(windows, {}, &Window::end)->end;
    std::vector<int> owner(max_end, -1);
    for (std::size_t w = 0; w < windows.size(); ++w) {
        const Window& window = windows[w];
        for (int l = 0; l < window.end - window.start; ++l) {
            int global_l = window.start + l;
            int depth = std::min(l, window.end - window.start - 1 - l);
            int current = owner[global_l];
            if (current == -1) {
                owner[global_l] = static_cast<int>(w);
                continue;
            }
            const Window& current_window = windows[current];
            int current_l = global_l - current_window.start;
            int current_depth = std::min(current_l, current_window.end - 1 - global_l);
            if (depth > current_depth) {
                owner[global_l] = static_cast<int>(w);
            }
        }
    }
    return owner;
}

std::vector<bool> read_observed_loci(const fs::path& path, int n_loci, int expected) {
    std::ifstream input(path);
    if (!input.is_open()) {
        throw std::runtime_error(std::format("Failed to open observed loci: {}", path.string()));
    }
    std::vector<bool> observed(n_loci, false);
    int previous = -1;
    int count = 0;
    int l;
    while (input >> l) {
        if (l <= previous || l < 0 || l >= n_loci) {
            throw std::runtime_error(std::format("Invalid observed locus in {}: {}", path.string(), l));
        }
        observed[l] = true;
        previous = l;
        ++count;
    }
    if (!input.eof()) {
        throw std::runtime_error(std::format("Invalid observed loci file: {}", path.string()));
    }
    if (count != expected) {
        throw std::runtime_error(std::format("Observed-locus count does not match {}.", path.string()));
    }
    return observed;
}

void evaluate_window(
    const Window& window, int window_slot, const std::vector<int>& owner, std::vector<MacStats>& stats
) {
    ImputeProbReader prob_reader((window.dir / "probs.bin").c_str());
    const ImputeProbHeader& header = prob_reader.header();
    SeqArray truth{read_seq_file((window.dir / "target_masked_true.bin").c_str())};
    if (header.n_sequences != static_cast<std::uint32_t>(truth.N)
        || header.n_loci != static_cast<std::uint32_t>(truth.L)) {
        throw std::runtime_error(
            std::format("Probability and truth dimensions differ in {}.", window.dir.string())
        );
    }

    std::vector<RefCount> ref_counts = read_ref_counts(window.dir / "ref.vcf.gz");
    int n_ref_loci = window.end - window.start;
    if (ref_counts.size() != static_cast<std::size_t>(n_ref_loci)) {
        throw std::runtime_error(std::format("Reference VCF dimensions differ in {}.", window.dir.string()));
    }
    std::vector<bool> observed = read_observed_loci(
        window.dir / "observed_loci.txt", n_ref_loci, n_ref_loci - truth.L
    );

    std::vector<MaskedLocus> retained;
    int masked_l = 0;
    for (int l = 0; l < n_ref_loci; ++l) {
        if (observed[l]) { continue; }
        if (owner[window.start + l] == window_slot) {
            RefCount count = ref_counts[l];
            int ref_count = count.an - count.ac;
            int mac = std::min(count.ac, ref_count);
            retained.push_back({masked_l, mac, count.ac <= ref_count});
            if (mac >= static_cast<int>(stats.size())) { stats.resize(mac + 1); }
            ++stats[mac].n_loci;
        }
        ++masked_l;
    }

    std::vector<std::uint16_t> prob_row(header.n_loci);
    for (int i = 0; i < truth.N; ++i) {
        prob_reader.read_row(prob_row);
        for (const MaskedLocus& locus : retained) {
            std::uint16_t q = prob_row[locus.column];
            int y = truth(i, locus.column);
            if (!locus.alt_is_minor) {
                q = static_cast<std::uint16_t>(fixed_point_max - q);
                y = 1 - y;
            }
            stats[locus.mac].add(q, y);
        }
    }
    std::cerr << "window=" << std::format("{:04d}", window.index)
        << " retained_loci=" << retained.size() << '\n';
}

void write_tsv(const fs::path& path, const std::vector<MacStats>& stats) {
    AtomicBinaryWriter output(path.c_str());
    std::ostream& stream = output.stream();
    stream << "mac\tn_loci\tn_predictions\tr2\taccuracy\n" << std::setprecision(10);
    for (std::size_t mac = 0; mac < stats.size(); ++mac) {
        const MacStats& value = stats[mac];
        if (value.n_predictions == 0) { continue; }
        stream << mac << '\t' << value.n_loci << '\t' << value.n_predictions << '\t'
            << value.r2() << '\t' << value.accuracy() << '\n';
    }
    output.finish();
}

}


int main(int argc, char* argv[]) {
    if (argc != 3) {
        throw std::invalid_argument("Usage: eval_impute WINDOWS_DIR OUTPUT.tsv");
    }

    fs::path root = argv[1];
    std::vector<Window> windows = read_windows(root);
    std::vector<int> owner = assign_loci(windows);
    std::vector<MacStats> stats;
    for (std::size_t w = 0; w < windows.size(); ++w) {
        evaluate_window(windows[w], static_cast<int>(w), owner, stats);
    }
    write_tsv(argv[2], stats);

    std::uint64_t n_loci = 0;
    std::uint64_t n_predictions = 0;
    for (const MacStats& value : stats) {
        n_loci += value.n_loci;
        n_predictions += value.n_predictions;
    }
    std::cerr << "windows=" << windows.size() << " retained_loci=" << n_loci
        << " predictions=" << n_predictions << " output=" << argv[2] << '\n';
    return 0;
}
