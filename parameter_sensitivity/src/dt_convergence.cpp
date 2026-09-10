#include "nfkb_model.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

struct Config {
    std::vector<double> dts{0.002, 0.001, 0.0005};
    double final_time = 20.0;
    int replicates = 100;
    std::uint64_t seed = 20260907;
    std::string scenario = "combined";
    double fucoxanthin = 0.0;
    double sigma = 0.01;
    double alpha_scale = 1.0;
    double beta_scale = 1.0;
    double rb = 0.75;
    bool random_weights = true;
    double fixed_weight = 0.5;
    bool feedback_on = false;
    bool clamp_receptors = true;
    bool saturating_inputs = true;
    std::string out_prefix = "results/dt_convergence/control";
    int progress_every = 10;
};

struct LevelAccumulator {
    double dt = 0.0;
    int ratio_to_finest = 1;
    std::vector<double> final_sum;
    std::vector<double> final_sumsq;
    std::vector<double> auc_sum;
    std::vector<double> auc_sumsq;
    std::vector<double> final_sqdiff_finest;
    std::vector<double> auc_sqdiff_finest;
    nfkb::StepDiagnostics diagnostics;
};

std::string get_arg(int argc, char** argv, const std::string& key,
                    const std::string& default_value) {
    for (int i = 1; i + 1 < argc; ++i) {
        if (argv[i] == key) return argv[i + 1];
    }
    return default_value;
}

bool has_arg(int argc, char** argv, const std::string& key) {
    for (int i = 1; i < argc; ++i) {
        if (argv[i] == key) return true;
    }
    return false;
}

bool to_bool(const std::string& value) {
    return value == "1" || value == "true" || value == "TRUE" ||
           value == "on" || value == "yes";
}

std::vector<double> parse_dts(const std::string& text) {
    std::vector<double> values;
    std::stringstream stream(text);
    std::string token;
    while (std::getline(stream, token, ',')) {
        if (!token.empty()) values.push_back(std::stod(token));
    }
    if (values.size() < 2) throw std::invalid_argument("provide at least two time steps");
    std::sort(values.begin(), values.end(), std::greater<double>());
    values.erase(std::unique(values.begin(), values.end(), [](double a, double b) {
        return std::abs(a - b) <= 1e-14 * std::max({1.0, std::abs(a), std::abs(b)});
    }), values.end());
    return values;
}

double safe_cv(double mean, double sd) {
    return std::abs(mean) > 1e-12
        ? sd / std::abs(mean)
        : std::numeric_limits<double>::quiet_NaN();
}

void ensure_parent(const std::string& prefix) {
    const fs::path parent = fs::path(prefix).parent_path();
    if (!parent.empty()) fs::create_directories(parent);
}

void print_help() {
    std::cout << R"HELP(
Coupled time-step convergence test for the full NF-kB bounded SDE

Fine Brownian increments are summed to construct every coarser increment.  All
levels therefore use the same interaction weights and the same Brownian paths.

Options:
  --dts LIST              comma-separated time steps, default 0.002,0.001,0.0005
  --time X                final time, default 20
  --mc N                  Monte Carlo replicates, default 100
  --seed N                base seed
  --scenario NAME         canonical|noncanonical|combined|none
  --fucox X               Fucoxanthin level, default 0
  --sigma X               baseline noise, default 0.01
  --alpha-scale X         default 1
  --beta-scale X          default 1
  --rb X                  Waller-Kraft r, default 0.75
  --weights random|fixed  default random
  --fixed-weight X        default 0.5
  --feedback 0|1          default 0
  --clamp-receptors 0|1   default 1
  --saturating-inputs 0|1 default 1
  --progress-every N      default 10
  --out-prefix PATH       output prefix

Every listed dt must be an integer multiple of the finest dt, and the final time
must be an integer multiple of every dt.
)HELP";
}

} // namespace

int main(int argc, char** argv) {
    if (has_arg(argc, argv, "--help") || has_arg(argc, argv, "-h")) {
        print_help();
        return 0;
    }

    try {
        Config cfg;
        cfg.dts = parse_dts(get_arg(argc, argv, "--dts", "0.002,0.001,0.0005"));
        cfg.final_time = std::stod(get_arg(argc, argv, "--time", std::to_string(cfg.final_time)));
        cfg.replicates = std::stoi(get_arg(argc, argv, "--mc", std::to_string(cfg.replicates)));
        cfg.seed = std::stoull(get_arg(argc, argv, "--seed", std::to_string(cfg.seed)));
        cfg.scenario = get_arg(argc, argv, "--scenario", cfg.scenario);
        cfg.fucoxanthin = std::stod(get_arg(argc, argv, "--fucox", std::to_string(cfg.fucoxanthin)));
        cfg.sigma = std::stod(get_arg(argc, argv, "--sigma", std::to_string(cfg.sigma)));
        cfg.alpha_scale = std::stod(get_arg(argc, argv, "--alpha-scale", std::to_string(cfg.alpha_scale)));
        cfg.beta_scale = std::stod(get_arg(argc, argv, "--beta-scale", std::to_string(cfg.beta_scale)));
        cfg.rb = std::stod(get_arg(argc, argv, "--rb", std::to_string(cfg.rb)));
        cfg.fixed_weight = std::stod(get_arg(argc, argv, "--fixed-weight", std::to_string(cfg.fixed_weight)));
        cfg.feedback_on = to_bool(get_arg(argc, argv, "--feedback", cfg.feedback_on ? "1" : "0"));
        cfg.clamp_receptors = to_bool(get_arg(argc, argv, "--clamp-receptors", cfg.clamp_receptors ? "1" : "0"));
        cfg.saturating_inputs = to_bool(get_arg(argc, argv, "--saturating-inputs", cfg.saturating_inputs ? "1" : "0"));
        cfg.progress_every = std::stoi(get_arg(argc, argv, "--progress-every", std::to_string(cfg.progress_every)));
        cfg.out_prefix = get_arg(argc, argv, "--out-prefix", cfg.out_prefix);

        const std::string weights = get_arg(argc, argv, "--weights", "random");
        if (weights == "random") cfg.random_weights = true;
        else if (weights == "fixed") cfg.random_weights = false;
        else throw std::invalid_argument("--weights must be random or fixed");

        if (cfg.final_time <= 0.0 || cfg.replicates <= 0 || cfg.sigma < 0.0) {
            throw std::invalid_argument("time and mc must be positive; sigma must be nonnegative");
        }
        if (cfg.progress_every < 0) throw std::invalid_argument("progress-every cannot be negative");

        const double finest_dt = cfg.dts.back();
        std::vector<int> ratios;
        ratios.reserve(cfg.dts.size());
        for (double dt : cfg.dts) {
            if (dt <= 0.0) throw std::invalid_argument("all time steps must be positive");
            const double ratio_real = dt / finest_dt;
            const int ratio = static_cast<int>(std::llround(ratio_real));
            if (ratio < 1 || std::abs(ratio_real - ratio) > 1e-10) {
                throw std::invalid_argument("every dt must be an integer multiple of the finest dt");
            }
            const double steps_real = cfg.final_time / dt;
            if (std::abs(steps_real - std::llround(steps_real)) > 1e-9) {
                throw std::invalid_argument("final time must be an integer multiple of every dt");
            }
            ratios.push_back(ratio);
        }
        const int finest_steps = static_cast<int>(std::llround(cfg.final_time / finest_dt));

        std::vector<std::unique_ptr<nfkb::NFkBModel>> models;
        std::vector<LevelAccumulator> levels(cfg.dts.size());
        for (std::size_t level = 0; level < cfg.dts.size(); ++level) {
            nfkb::SimulationConfig model_cfg;
            model_cfg.dt = cfg.dts[level];
            model_cfg.steps = static_cast<int>(std::llround(cfg.final_time / cfg.dts[level]));
            model_cfg.scenario = cfg.scenario;
            model_cfg.fucoxanthin = cfg.fucoxanthin;
            model_cfg.sigma = cfg.sigma;
            model_cfg.alpha_scale = cfg.alpha_scale;
            model_cfg.beta_scale = cfg.beta_scale;
            model_cfg.rb = cfg.rb;
            model_cfg.random_weights = cfg.random_weights;
            model_cfg.fixed_weight = cfg.fixed_weight;
            model_cfg.feedback_on = cfg.feedback_on;
            model_cfg.clamp_receptors = cfg.clamp_receptors;
            model_cfg.saturating_inputs = cfg.saturating_inputs;
            model_cfg.projected_em = true;
            model_cfg.replicates = 1;
            models.push_back(std::make_unique<nfkb::NFkBModel>(model_cfg));

            levels[level].dt = cfg.dts[level];
            levels[level].ratio_to_finest = ratios[level];
            const std::size_t n = models[level]->node_count();
            levels[level].final_sum.assign(n, 0.0);
            levels[level].final_sumsq.assign(n, 0.0);
            levels[level].auc_sum.assign(n, 0.0);
            levels[level].auc_sumsq.assign(n, 0.0);
            levels[level].final_sqdiff_finest.assign(n, 0.0);
            levels[level].auc_sqdiff_finest.assign(n, 0.0);
        }

        const std::size_t node_count = models.front()->node_count();
        for (const auto& model : models) {
            if (model->node_count() != node_count) {
                throw std::runtime_error("time-step models have different node counts");
            }
        }

        std::mt19937_64 master_rng(cfg.seed);
        for (int rep = 0; rep < cfg.replicates; ++rep) {
            const std::uint64_t stream_seed = master_rng();
            std::mt19937_64 rng(stream_seed);
            const auto weight_realization = models.front()->draw_weight_realization(rng);
            for (auto& model : models) model->set_weight_realization(weight_realization);

            std::normal_distribution<double> normal(0.0, 1.0);
            std::vector<std::vector<double>> state;
            std::vector<std::vector<double>> auc;
            std::vector<std::vector<double>> accumulated_dW(
                models.size(), std::vector<double>(node_count, 0.0));
            std::vector<int> counters(models.size(), 0);
            std::vector<nfkb::StepDiagnostics> rep_diagnostics(models.size());

            state.reserve(models.size());
            auc.reserve(models.size());
            for (const auto& model : models) {
                state.push_back(model->initial_state());
                auc.emplace_back(node_count, 0.0);
            }

            std::vector<double> fine_dW(node_count, 0.0);
            for (int fine_step = 0; fine_step < finest_steps; ++fine_step) {
                for (std::size_t node = 0; node < node_count; ++node) {
                    fine_dW[node] = std::sqrt(finest_dt) * normal(rng);
                }

                for (std::size_t level = 0; level < models.size(); ++level) {
                    for (std::size_t node = 0; node < node_count; ++node) {
                        accumulated_dW[level][node] += fine_dW[node];
                    }
                    ++counters[level];
                    if (counters[level] == levels[level].ratio_to_finest) {
                        for (std::size_t node = 0; node < node_count; ++node) {
                            auc[level][node] += state[level][node] * levels[level].dt;
                        }
                        state[level] = models[level]->step_with_dW(
                            state[level], accumulated_dW[level], rep_diagnostics[level]);
                        std::fill(accumulated_dW[level].begin(),
                                  accumulated_dW[level].end(), 0.0);
                        counters[level] = 0;
                    }
                }
            }

            for (int counter : counters) {
                if (counter != 0) throw std::runtime_error("internal Brownian aggregation error");
            }

            const std::size_t finest_level = models.size() - 1;
            for (std::size_t level = 0; level < models.size(); ++level) {
                levels[level].diagnostics.add(rep_diagnostics[level]);
                for (std::size_t node = 0; node < node_count; ++node) {
                    const double final_value = state[level][node];
                    const double auc_value = auc[level][node];
                    levels[level].final_sum[node] += final_value;
                    levels[level].final_sumsq[node] += final_value * final_value;
                    levels[level].auc_sum[node] += auc_value;
                    levels[level].auc_sumsq[node] += auc_value * auc_value;

                    const double final_diff = final_value - state[finest_level][node];
                    const double auc_diff = auc_value - auc[finest_level][node];
                    levels[level].final_sqdiff_finest[node] += final_diff * final_diff;
                    levels[level].auc_sqdiff_finest[node] += auc_diff * auc_diff;
                }
            }

            if (cfg.progress_every > 0 &&
                ((rep + 1) % cfg.progress_every == 0 || rep + 1 == cfg.replicates)) {
                std::cout << "Progress: " << rep + 1 << " / " << cfg.replicates << "\n";
            }
        }

        ensure_parent(cfg.out_prefix);
        std::ofstream summary(cfg.out_prefix + "_summary.csv");
        std::ofstream diagnostics(cfg.out_prefix + "_diagnostics.csv");
        std::ofstream settings(cfg.out_prefix + "_settings.csv");
        if (!summary || !diagnostics || !settings) {
            throw std::runtime_error("cannot create convergence output files");
        }

        summary << "dt,finest_dt,node,final_mean,final_sd,final_cv,auc_mean,auc_sd,"
                << "final_mean_abs_difference_vs_finest,final_mean_relative_difference_vs_finest,"
                << "final_pathwise_rmse_vs_finest,auc_mean_abs_difference_vs_finest,"
                << "auc_mean_relative_difference_vs_finest,auc_pathwise_rmse_vs_finest,"
                << "final_cv_abs_difference_vs_finest\n";
        diagnostics << "dt,replicates,predictor_violations,predictor_values,violation_rate,"
                    << "max_lower_overshoot,max_upper_overshoot,mean_overshoot_per_value\n";
        settings << "scenario,fucoxanthin,final_time,replicates,seed,sigma,alpha_scale,beta_scale,"
                 << "r,random_weights,fixed_weight,feedback_on,clamp_receptors,saturating_inputs,dts\n";

        summary << std::setprecision(12);
        diagnostics << std::setprecision(12);
        settings << std::setprecision(12);

        const double R = static_cast<double>(cfg.replicates);
        const std::size_t finest_level = levels.size() - 1;
        std::vector<double> finest_final_mean(node_count, 0.0);
        std::vector<double> finest_auc_mean(node_count, 0.0);
        std::vector<double> finest_final_cv(node_count, 0.0);

        for (std::size_t node = 0; node < node_count; ++node) {
            finest_final_mean[node] = levels[finest_level].final_sum[node] / R;
            finest_auc_mean[node] = levels[finest_level].auc_sum[node] / R;
            const double second = levels[finest_level].final_sumsq[node] / R;
            const double sd = std::sqrt(std::max(
                0.0, second - finest_final_mean[node] * finest_final_mean[node]));
            finest_final_cv[node] = safe_cv(finest_final_mean[node], sd);
        }

        const auto& node_names = models.front()->nodes();
        for (std::size_t level = 0; level < levels.size(); ++level) {
            for (std::size_t node = 0; node < node_count; ++node) {
                const double final_mean = levels[level].final_sum[node] / R;
                const double final_second = levels[level].final_sumsq[node] / R;
                const double final_sd = std::sqrt(std::max(
                    0.0, final_second - final_mean * final_mean));
                const double final_cv = safe_cv(final_mean, final_sd);
                const double auc_mean = levels[level].auc_sum[node] / R;
                const double auc_second = levels[level].auc_sumsq[node] / R;
                const double auc_sd = std::sqrt(std::max(
                    0.0, auc_second - auc_mean * auc_mean));

                const double final_abs = std::abs(final_mean - finest_final_mean[node]);
                const double final_rel = final_abs /
                    std::max(std::abs(finest_final_mean[node]), 1e-12);
                const double final_rmse = std::sqrt(
                    levels[level].final_sqdiff_finest[node] / R);
                const double auc_abs = std::abs(auc_mean - finest_auc_mean[node]);
                const double auc_rel = auc_abs /
                    std::max(std::abs(finest_auc_mean[node]), 1e-12);
                const double auc_rmse = std::sqrt(
                    levels[level].auc_sqdiff_finest[node] / R);
                const double cv_abs = std::abs(final_cv - finest_final_cv[node]);

                summary << levels[level].dt << ',' << finest_dt << ','
                        << node_names[node].name << ',' << final_mean << ',' << final_sd << ','
                        << final_cv << ',' << auc_mean << ',' << auc_sd << ',' << final_abs << ','
                        << final_rel << ',' << final_rmse << ',' << auc_abs << ',' << auc_rel << ','
                        << auc_rmse << ',' << cv_abs << '\n';
            }

            const auto& d = levels[level].diagnostics;
            const double violation_rate = d.predictor_values > 0
                ? static_cast<double>(d.predictor_violations) / d.predictor_values
                : 0.0;
            const double mean_overshoot = d.predictor_values > 0
                ? d.total_overshoot / d.predictor_values
                : 0.0;
            diagnostics << levels[level].dt << ',' << cfg.replicates << ','
                        << d.predictor_violations << ',' << d.predictor_values << ','
                        << violation_rate << ',' << d.max_lower_overshoot << ','
                        << d.max_upper_overshoot << ',' << mean_overshoot << '\n';
        }

        std::ostringstream dt_text;
        for (std::size_t i = 0; i < cfg.dts.size(); ++i) {
            if (i > 0) dt_text << ';';
            dt_text << cfg.dts[i];
        }
        settings << cfg.scenario << ',' << cfg.fucoxanthin << ',' << cfg.final_time << ','
                 << cfg.replicates << ',' << cfg.seed << ',' << cfg.sigma << ','
                 << cfg.alpha_scale << ',' << cfg.beta_scale << ',' << cfg.rb << ','
                 << (cfg.random_weights ? 1 : 0) << ',' << cfg.fixed_weight << ','
                 << (cfg.feedback_on ? 1 : 0) << ',' << (cfg.clamp_receptors ? 1 : 0) << ','
                 << (cfg.saturating_inputs ? 1 : 0) << ',' << dt_text.str() << '\n';

        std::cout << "Saved coupled time-step convergence results to "
                  << cfg.out_prefix << "_*.csv\n";
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << "\n\n";
        print_help();
        return 1;
    }

    return 0;
}
