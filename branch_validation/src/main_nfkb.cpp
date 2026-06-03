#include "nfkb_model.hpp"

#include <sstream>

using nfkb::NFkBModel;
using nfkb::SimulationConfig;

namespace {

bool has_arg(int argc, char** argv, const std::string& key) {
    for (int i = 1; i < argc; ++i) if (argv[i] == key) return true;
    return false;
}

std::string get_arg(int argc, char** argv, const std::string& key, const std::string& def) {
    for (int i = 1; i + 1 < argc; ++i) {
        if (argv[i] == key) return argv[i + 1];
    }
    return def;
}

int to_int(const std::string& s) { return std::stoi(s); }
double to_double(const std::string& s) { return std::stod(s); }
bool to_bool(const std::string& s) {
    return s == "1" || s == "true" || s == "TRUE" || s == "on" || s == "yes";
}

void print_help() {
    std::cout << R"HELP(
NF-kB bounded stochastic extended Boolean simulator

Usage:
  nfkb_sim [options]

Core options:
  --scenario canonical|noncanonical|combined|none   default: canonical
  --steps N                                         integration steps, default: 1000
  --dt DT                                           integration time step, default: 0.01
  --save-every N                                    save one CSV row every N steps, default: 1
  --mc N                                            Monte Carlo replicates, default: 1
  --progress-every N                                print progress every N replicates, default: 0
  --sigma S                                         node noise intensity, default: 0.01
  --seed N                                          default: 12345
  --out-prefix PREFIX                               default: results/nfkb

Model switches:
  --weights fixed|random                            default: fixed
  --fixed-weight W                                  default: 0.5
  --feedback 0|1                                    default: 0
  --clamp-receptors 0|1                             default: 1
  --projected-em 0|1                                default: 1
  --saturating-inputs 0|1                           default: 1
  --isolate-branches 0|1                            default: 0
  --cross-branch-scale X                            default: 0.0
  --fucox C                                         Fucoxanthin concentration, default: 0

Examples:
  ./nfkb_sim --scenario canonical --sigma 0 --steps 500 --dt 0.02 --out-prefix results/canonical_det
  ./nfkb_sim --scenario noncanonical --isolate-branches 1 --sigma 0 --out-prefix results/noncanonical_iso_det
  ./nfkb_sim --scenario combined --sigma 0.01 --mc 200 --weights random --steps 100000 --dt 0.001 --save-every 100 --progress-every 20 --out-prefix results/combined_long_test
  ./nfkb_sim --scenario combined --fucox 50 --sigma 0.01 --mc 1000 --weights random --steps 200000 --dt 0.001 --save-every 100 --progress-every 50 --out-prefix results/fucox50_long

Output files:
  PREFIX_mean.csv          Monte Carlo mean trajectory, saved sparsely if --save-every > 1
  PREFIX_sd.csv            Monte Carlo standard deviation trajectory
  PREFIX_rep0.csv          representative trajectory from replicate 0
  PREFIX_summary.csv       final value and AUC summaries over the full integration horizon
  PREFIX_network.csv       node-edge list actually used by the simulator
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
        cfg.seed = static_cast<unsigned int>(to_int(get_arg(argc, argv, "--seed", std::to_string(cfg.seed))));
        cfg.out_prefix = get_arg(argc, argv, "--out-prefix", cfg.out_prefix);
        cfg.fixed_weight = to_double(get_arg(argc, argv, "--fixed-weight", std::to_string(cfg.fixed_weight)));
        cfg.feedback_on = to_bool(get_arg(argc, argv, "--feedback", cfg.feedback_on ? "1" : "0"));
        cfg.clamp_receptors = to_bool(get_arg(argc, argv, "--clamp-receptors", cfg.clamp_receptors ? "1" : "0"));
        cfg.projected_em = to_bool(get_arg(argc, argv, "--projected-em", cfg.projected_em ? "1" : "0"));
        cfg.saturating_inputs = to_bool(get_arg(argc, argv, "--saturating-inputs", cfg.saturating_inputs ? "1" : "0"));
        cfg.isolate_branches = to_bool(get_arg(argc, argv, "--isolate-branches", cfg.isolate_branches ? "1" : "0"));
        cfg.cross_branch_scale = to_double(get_arg(argc, argv, "--cross-branch-scale", std::to_string(cfg.cross_branch_scale)));
        cfg.fucoxanthin = to_double(get_arg(argc, argv, "--fucox", std::to_string(cfg.fucoxanthin)));
        cfg.sigma = to_double(get_arg(argc, argv, "--sigma", std::to_string(cfg.sigma)));

        if (cfg.save_every <= 0) throw std::invalid_argument("--save-every must be positive");
        if (cfg.progress_every < 0) throw std::invalid_argument("--progress-every cannot be negative");

        const std::string weight_mode = get_arg(argc, argv, "--weights", cfg.random_weights ? "random" : "fixed");
        if (weight_mode == "fixed") cfg.random_weights = false;
        else if (weight_mode == "random") cfg.random_weights = true;
        else throw std::invalid_argument("--weights must be fixed or random");

        NFkBModel model(cfg);
        auto result = model.run();
        nfkb::write_timecourse_csv(cfg.out_prefix + "_mean.csv", result.node_names, result.time, result.mean);
        nfkb::write_timecourse_csv(cfg.out_prefix + "_sd.csv", result.node_names, result.time, result.sd);
        nfkb::write_timecourse_csv(cfg.out_prefix + "_rep0.csv", result.node_names, result.time, result.rep0);
        nfkb::write_summary_csv(cfg.out_prefix + "_summary.csv", result);
        nfkb::write_network_csv(cfg.out_prefix + "_network.csv", model.nodes(), model.edges());

        std::cout << "Simulation complete.\n";
        std::cout << "Scenario: " << cfg.scenario << "\n";
        std::cout << "Replicates: " << cfg.replicates
                  << ", steps: " << cfg.steps
                  << ", dt: " << cfg.dt
                  << ", save_every: " << cfg.save_every << "\n";
        std::cout << "Saved time points: " << result.time.size() << "\n";
        std::cout << "Predictor boundary violations before projection: "
                  << result.predictor_violations << " / " << result.predictor_values << "\n";
        std::cout << "Output prefix: " << cfg.out_prefix << "\n";
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n\n";
        print_help();
        return 1;
    }

    return 0;
}
