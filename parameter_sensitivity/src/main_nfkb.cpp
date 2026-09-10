#include "nfkb_model.hpp"

#include <cstdint>

using nfkb::NFkBModel;
using nfkb::SimulationConfig;

namespace {

bool has_arg(int argc, char** argv, const std::string& key) {
    for (int i = 1; i < argc; ++i) {
        if (argv[i] == key) return true;
    }
    return false;
}

std::string get_arg(int argc, char** argv, const std::string& key,
                    const std::string& default_value) {
    for (int i = 1; i + 1 < argc; ++i) {
        if (argv[i] == key) return argv[i + 1];
    }
    return default_value;
}

int to_int(const std::string& value) { return std::stoi(value); }
double to_double(const std::string& value) { return std::stod(value); }
std::uint64_t to_uint64(const std::string& value) { return std::stoull(value); }

bool to_bool(const std::string& value) {
    return value == "1" || value == "true" || value == "TRUE" ||
           value == "on" || value == "yes";
}

void print_help() {
    std::cout << R"HELP(
NF-kB revision simulator

This executable is the paper's NF-kB bounded-SDE code with explicit controls for
activation-rate, deactivation-rate, noise-intensity, and Waller-Kraft sensitivity.

Usage:
  nfkb_revision_sim [options]

Core options:
  --scenario canonical|noncanonical|combined|none   default: canonical
  --steps N                                         default: 1000
  --dt DT                                           default: 0.01
  --save-every N                                    default: 1
  --mc N                                            default: 1
  --progress-every N                                default: 0
  --seed N                                          default: 12345
  --out-prefix PREFIX                               default: results/nfkb

Sensitivity parameters:
  --alpha-scale X                                   global multiplier, default: 1
  --beta-scale X                                    global multiplier, default: 1
  --sigma S                                         baseline noise, default: 0.01
  --sigma-scale X                                   multiplier, default: 1
  --rb R                                            Waller-Kraft r in [0.5,1], default: 0.75

Model switches:
  --weights fixed|random                            default: fixed
  --fixed-weight W                                  default: 0.5
  --feedback 0|1                                    default: 0
  --clamp-receptors 0|1                             default: 1
  --projected-em 0|1                                default: 1
  --saturating-inputs 0|1                           default: 1
  --isolate-branches 0|1                            default: 0
  --cross-branch-scale X                            default: 0
  --fucox C                                         default: 0

Output controls:
  --write-timecourse 0|1                            default: 1
  --write-replicates 0|1                            default: 0

Files:
  PREFIX_summary.csv       ensemble final/AUC statistics
  PREFIX_diagnostics.csv   run parameters and projection diagnostics
  PREFIX_replicates.csv    optional replicate-level final/AUC values
  PREFIX_mean.csv          optional mean time course
  PREFIX_sd.csv            optional SD time course
  PREFIX_rep0.csv          optional representative trajectory
  PREFIX_network.csv       network topology

Example:
  ./nfkb_revision_sim --scenario combined --weights random --mc 100 \
    --steps 50000 --dt 0.001 --save-every 50000 --fucox 100 \
    --alpha-scale 1.25 --beta-scale 1 --sigma 0.01 --sigma-scale 1 \
    --rb 0.75 --write-timecourse 0 --write-replicates 1 \
    --out-prefix results/example
)HELP";
}

} // namespace

int main(int argc, char** argv) {
    if (has_arg(argc, argv, "--help") || has_arg(argc, argv, "-h")) {
        print_help();
        return 0;
    }

    try {
        SimulationConfig cfg;
        cfg.scenario = get_arg(argc, argv, "--scenario", cfg.scenario);
        cfg.steps = to_int(get_arg(argc, argv, "--steps", std::to_string(cfg.steps)));
        cfg.dt = to_double(get_arg(argc, argv, "--dt", std::to_string(cfg.dt)));
        cfg.save_every = to_int(get_arg(argc, argv, "--save-every", std::to_string(cfg.save_every)));
        cfg.replicates = to_int(get_arg(argc, argv, "--mc", std::to_string(cfg.replicates)));
        cfg.progress_every = to_int(get_arg(argc, argv, "--progress-every", std::to_string(cfg.progress_every)));
        cfg.seed = to_uint64(get_arg(argc, argv, "--seed", std::to_string(cfg.seed)));
        cfg.out_prefix = get_arg(argc, argv, "--out-prefix", cfg.out_prefix);

        cfg.alpha_scale = to_double(get_arg(argc, argv, "--alpha-scale", std::to_string(cfg.alpha_scale)));
        cfg.beta_scale = to_double(get_arg(argc, argv, "--beta-scale", std::to_string(cfg.beta_scale)));
        cfg.sigma = to_double(get_arg(argc, argv, "--sigma", std::to_string(cfg.sigma)));
        cfg.sigma_scale = to_double(get_arg(argc, argv, "--sigma-scale", std::to_string(cfg.sigma_scale)));
        cfg.rb = to_double(get_arg(argc, argv, "--rb", std::to_string(cfg.rb)));

        cfg.fixed_weight = to_double(get_arg(argc, argv, "--fixed-weight", std::to_string(cfg.fixed_weight)));
        cfg.feedback_on = to_bool(get_arg(argc, argv, "--feedback", cfg.feedback_on ? "1" : "0"));
        cfg.clamp_receptors = to_bool(get_arg(argc, argv, "--clamp-receptors", cfg.clamp_receptors ? "1" : "0"));
        cfg.projected_em = to_bool(get_arg(argc, argv, "--projected-em", cfg.projected_em ? "1" : "0"));
        cfg.saturating_inputs = to_bool(get_arg(argc, argv, "--saturating-inputs", cfg.saturating_inputs ? "1" : "0"));
        cfg.isolate_branches = to_bool(get_arg(argc, argv, "--isolate-branches", cfg.isolate_branches ? "1" : "0"));
        cfg.cross_branch_scale = to_double(get_arg(argc, argv, "--cross-branch-scale", std::to_string(cfg.cross_branch_scale)));
        cfg.fucoxanthin = to_double(get_arg(argc, argv, "--fucox", std::to_string(cfg.fucoxanthin)));
        cfg.write_timecourse = to_bool(get_arg(argc, argv, "--write-timecourse", cfg.write_timecourse ? "1" : "0"));
        cfg.write_replicates = to_bool(get_arg(argc, argv, "--write-replicates", cfg.write_replicates ? "1" : "0"));

        if (cfg.steps <= 0) throw std::invalid_argument("--steps must be positive");
        if (cfg.dt <= 0.0) throw std::invalid_argument("--dt must be positive");
        if (cfg.save_every <= 0) throw std::invalid_argument("--save-every must be positive");
        if (cfg.replicates <= 0) throw std::invalid_argument("--mc must be positive");
        if (cfg.progress_every < 0) throw std::invalid_argument("--progress-every cannot be negative");

        const std::string weight_mode = get_arg(
            argc, argv, "--weights", cfg.random_weights ? "random" : "fixed");
        if (weight_mode == "fixed") {
            cfg.random_weights = false;
        } else if (weight_mode == "random") {
            cfg.random_weights = true;
        } else {
            throw std::invalid_argument("--weights must be fixed or random");
        }

        if (!cfg.random_weights && std::abs(cfg.rb - 0.75) > 1e-14) {
            std::cerr
                << "[warning] With identical fixed scalar edge weights, changing r can be "
                << "numerically degenerate. Use --weights random for the r sensitivity scan.\n";
        }

        NFkBModel model(cfg);
        const auto result = model.run();

        if (cfg.write_timecourse) {
            nfkb::write_timecourse_csv(cfg.out_prefix + "_mean.csv", result.node_names,
                                       result.time, result.mean);
            nfkb::write_timecourse_csv(cfg.out_prefix + "_sd.csv", result.node_names,
                                       result.time, result.sd);
            nfkb::write_timecourse_csv(cfg.out_prefix + "_rep0.csv", result.node_names,
                                       result.time, result.rep0);
        }
        nfkb::write_summary_csv(cfg.out_prefix + "_summary.csv", result);
        nfkb::write_diagnostics_csv(cfg.out_prefix + "_diagnostics.csv", cfg, result);
        if (cfg.write_replicates) {
            nfkb::write_replicates_csv(cfg.out_prefix + "_replicates.csv", result);
        }
        nfkb::write_network_csv(cfg.out_prefix + "_network.csv", model.nodes(), model.edges());

        const double violation_rate = result.diagnostics.predictor_values > 0
            ? static_cast<double>(result.diagnostics.predictor_violations) /
              static_cast<double>(result.diagnostics.predictor_values)
            : 0.0;

        std::cout << "Simulation complete.\n"
                  << "Scenario: " << cfg.scenario << "\n"
                  << "Replicates: " << cfg.replicates
                  << ", steps: " << cfg.steps
                  << ", dt: " << cfg.dt
                  << ", final time: " << cfg.steps * cfg.dt << "\n"
                  << "alpha scale: " << cfg.alpha_scale
                  << ", beta scale: " << cfg.beta_scale
                  << ", sigma effective: " << cfg.sigma * cfg.sigma_scale
                  << ", Waller-Kraft r: " << cfg.rb << "\n"
                  << "Predictor violations before projection: "
                  << result.diagnostics.predictor_violations << " / "
                  << result.diagnostics.predictor_values
                  << " (rate=" << violation_rate << ")\n"
                  << "Maximum lower/upper overshoot: "
                  << result.diagnostics.max_lower_overshoot << " / "
                  << result.diagnostics.max_upper_overshoot << "\n"
                  << "Output prefix: " << cfg.out_prefix << "\n";
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n\n";
        print_help();
        return 1;
    }

    return 0;
}
