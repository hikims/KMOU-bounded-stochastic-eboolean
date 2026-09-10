#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

struct Config {
    std::vector<double> inputs{0.01, 0.03, 0.10, 0.30, 1.0, 3.0};
    double dt = 0.001;
    double final_time = 20.0;
    int replicates = 5000;
    int save_every = 100;
    double sigma = 0.7;
    double gain = 1.0;
    double d1 = 1.0;
    double d2 = 1.0;
    std::uint64_t seed = 20260907;
    std::string out_dir = "results/low_activation";
};

struct Diagnostics {
    long long bounded_predictor_violations = 0;
    long long bounded_predictor_values = 0;
    double bounded_max_lower_overshoot = 0.0;
    double bounded_max_upper_overshoot = 0.0;
    int linear_exit_trajectories = 0;
};

struct TimeStats {
    std::vector<std::array<double, 2>> sum;
    std::vector<std::array<double, 2>> sumsq;
};

using MomentState = std::array<double, 5>; // m1, m2, Q11, Q12, Q22

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

std::vector<double> parse_list(const std::string& text) {
    std::vector<double> values;
    std::stringstream stream(text);
    std::string token;
    while (std::getline(stream, token, ',')) {
        if (!token.empty()) values.push_back(std::stod(token));
    }
    if (values.empty()) throw std::invalid_argument("the input list is empty");
    return values;
}

double clamp01(double x) {
    return std::min(1.0, std::max(0.0, x));
}

double safe_cv(double mean, double sd) {
    return std::abs(mean) > 1e-12
        ? sd / std::abs(mean)
        : std::numeric_limits<double>::quiet_NaN();
}

MomentState moment_rhs(const MomentState& y, double input, const Config& cfg) {
    const double m1 = y[0];
    const double m2 = y[1];
    const double q11 = y[2];
    const double q12 = y[3];
    const double q22 = y[4];
    const double a1 = input + cfg.d1;
    const double sigma2 = cfg.sigma * cfg.sigma;

    return {
        input - a1 * m1,
        cfg.gain * m1 - cfg.d2 * m2,
        (-2.0 * a1 + sigma2) * q11 + 2.0 * input * m1,
        cfg.gain * q11 - (a1 + cfg.d2) * q12 + input * m2,
        2.0 * cfg.gain * q12 + (-2.0 * cfg.d2 + sigma2) * q22
    };
}

MomentState add_scaled(const MomentState& a, const MomentState& b, double scale) {
    MomentState out{};
    for (std::size_t i = 0; i < out.size(); ++i) out[i] = a[i] + scale * b[i];
    return out;
}

MomentState rk4_step(const MomentState& y, double dt, double input, const Config& cfg) {
    const auto k1 = moment_rhs(y, input, cfg);
    const auto k2 = moment_rhs(add_scaled(y, k1, 0.5 * dt), input, cfg);
    const auto k3 = moment_rhs(add_scaled(y, k2, 0.5 * dt), input, cfg);
    const auto k4 = moment_rhs(add_scaled(y, k3, dt), input, cfg);

    MomentState out{};
    for (std::size_t i = 0; i < out.size(); ++i) {
        out[i] = y[i] + dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
    }
    return out;
}

void print_help() {
    std::cout << R"HELP(
Low-activation validation for the two-node bounded cascade

The program compares, under identical Brownian increments:
  1. the full bounded two-node SDE,
  2. its first-order affine multiplicative-noise linearization, and
  3. the closed first/second moment ODEs of the linearized system.

Options:
  --inputs LIST       comma-separated input levels (default 0.01,0.03,0.1,0.3,1,3)
  --dt X              integration step (default 0.001)
  --time X            final time (default 20)
  --mc N              Monte Carlo replicates (default 5000)
  --save-every N      save every N steps (default 100)
  --sigma X           multiplicative-noise intensity (default 0.7)
  --gain X            feed-forward gain g (default 1)
  --d1 X              basal loss of node 1 (default 1)
  --d2 X              loss of node 2 (default 1)
  --seed N            random seed
  --out-dir PATH      output directory
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
        cfg.inputs = parse_list(get_arg(argc, argv, "--inputs", "0.01,0.03,0.1,0.3,1,3"));
        cfg.dt = std::stod(get_arg(argc, argv, "--dt", std::to_string(cfg.dt)));
        cfg.final_time = std::stod(get_arg(argc, argv, "--time", std::to_string(cfg.final_time)));
        cfg.replicates = std::stoi(get_arg(argc, argv, "--mc", std::to_string(cfg.replicates)));
        cfg.save_every = std::stoi(get_arg(argc, argv, "--save-every", std::to_string(cfg.save_every)));
        cfg.sigma = std::stod(get_arg(argc, argv, "--sigma", std::to_string(cfg.sigma)));
        cfg.gain = std::stod(get_arg(argc, argv, "--gain", std::to_string(cfg.gain)));
        cfg.d1 = std::stod(get_arg(argc, argv, "--d1", std::to_string(cfg.d1)));
        cfg.d2 = std::stod(get_arg(argc, argv, "--d2", std::to_string(cfg.d2)));
        cfg.seed = std::stoull(get_arg(argc, argv, "--seed", std::to_string(cfg.seed)));
        cfg.out_dir = get_arg(argc, argv, "--out-dir", cfg.out_dir);

        if (cfg.dt <= 0.0 || cfg.final_time <= 0.0) {
            throw std::invalid_argument("dt and final time must be positive");
        }
        if (cfg.replicates <= 0 || cfg.save_every <= 0) {
            throw std::invalid_argument("mc and save-every must be positive");
        }
        if (cfg.sigma < 0.0 || cfg.gain <= 0.0 || cfg.d1 <= 0.0 || cfg.d2 <= 0.0) {
            throw std::invalid_argument("sigma must be nonnegative and gain/loss rates positive");
        }
        for (double input : cfg.inputs) {
            if (input <= 0.0) throw std::invalid_argument("all input levels must be positive");
        }

        const int steps = static_cast<int>(std::llround(cfg.final_time / cfg.dt));
        if (std::abs(steps * cfg.dt - cfg.final_time) > 1e-10 * std::max(1.0, cfg.final_time)) {
            throw std::invalid_argument("final time must be an integer multiple of dt");
        }

        std::vector<int> save_steps;
        for (int step = 0; step <= steps; step += cfg.save_every) save_steps.push_back(step);
        if (save_steps.back() != steps) save_steps.push_back(steps);
        const int saved = static_cast<int>(save_steps.size());

        fs::create_directories(cfg.out_dir);
        std::ofstream time_out(fs::path(cfg.out_dir) / "low_activation_timecourse.csv");
        std::ofstream summary_out(fs::path(cfg.out_dir) / "low_activation_summary.csv");
        std::ofstream error_out(fs::path(cfg.out_dir) / "low_activation_errors.csv");
        std::ofstream diag_out(fs::path(cfg.out_dir) / "low_activation_diagnostics.csv");
        if (!time_out || !summary_out || !error_out || !diag_out) {
            throw std::runtime_error("cannot create one or more output files");
        }

        time_out << "input,time,model,node,mean,sd,cv\n";
        summary_out << "input,model,node,final_mean,final_sd,final_cv\n";
        error_out << "input,node,bounded_mean,linear_mc_mean,moment_mean,"
                  << "linear_vs_bounded_abs_error,linear_vs_bounded_rel_error,"
                  << "moment_vs_linear_abs_error,moment_vs_linear_rel_error,"
                  << "bounded_cv,linear_mc_cv,moment_cv\n";
        diag_out << "input,replicates,dt,final_time,sigma,gain,d1,d2,"
                 << "bounded_predictor_violations,bounded_predictor_values,violation_rate,"
                 << "bounded_max_lower_overshoot,bounded_max_upper_overshoot,"
                 << "linear_exit_trajectories,linear_exit_fraction\n";

        time_out << std::setprecision(12);
        summary_out << std::setprecision(12);
        error_out << std::setprecision(12);
        diag_out << std::setprecision(12);

        for (std::size_t input_index = 0; input_index < cfg.inputs.size(); ++input_index) {
            const double input = cfg.inputs[input_index];
            std::cout << "Input " << input << " (" << input_index + 1 << "/"
                      << cfg.inputs.size() << ")\n";

            TimeStats full{{}, {}};
            TimeStats linear{{}, {}};
            full.sum.assign(saved, {0.0, 0.0});
            full.sumsq.assign(saved, {0.0, 0.0});
            linear.sum.assign(saved, {0.0, 0.0});
            linear.sumsq.assign(saved, {0.0, 0.0});
            Diagnostics diagnostics;

            // The same replicate stream seeds are reused at every input level.
            // This common-random-number design reduces noise in input comparisons.
            std::mt19937_64 master_rng(cfg.seed);
            for (int rep = 0; rep < cfg.replicates; ++rep) {
                std::mt19937_64 rng(master_rng());
                std::normal_distribution<double> normal(0.0, 1.0);

                std::array<double, 2> x{0.0, 0.0};
                std::array<double, 2> y{0.0, 0.0};
                bool linear_exited = false;
                int next_save = 0;

                auto save = [&](int step) {
                    while (next_save < saved && save_steps[next_save] == step) {
                        for (int node = 0; node < 2; ++node) {
                            full.sum[next_save][node] += x[node];
                            full.sumsq[next_save][node] += x[node] * x[node];
                            linear.sum[next_save][node] += y[node];
                            linear.sumsq[next_save][node] += y[node] * y[node];
                        }
                        ++next_save;
                    }
                };
                save(0);

                for (int step = 0; step < steps; ++step) {
                    const double dW1 = std::sqrt(cfg.dt) * normal(rng);
                    const double dW2 = std::sqrt(cfg.dt) * normal(rng);

                    const double full_drift1 = input * (1.0 - x[0]) - cfg.d1 * x[0];
                    const double full_drift2 = cfg.gain * x[0] * (1.0 - x[1]) - cfg.d2 * x[1];
                    const double px1 = x[0] + full_drift1 * cfg.dt
                                     + cfg.sigma * x[0] * (1.0 - x[0]) * dW1;
                    const double px2 = x[1] + full_drift2 * cfg.dt
                                     + cfg.sigma * x[1] * (1.0 - x[1]) * dW2;

                    for (double predictor : {px1, px2}) {
                        ++diagnostics.bounded_predictor_values;
                        const double lower = std::max(0.0, -predictor);
                        const double upper = std::max(0.0, predictor - 1.0);
                        if (lower > 0.0 || upper > 0.0) {
                            ++diagnostics.bounded_predictor_violations;
                        }
                        diagnostics.bounded_max_lower_overshoot =
                            std::max(diagnostics.bounded_max_lower_overshoot, lower);
                        diagnostics.bounded_max_upper_overshoot =
                            std::max(diagnostics.bounded_max_upper_overshoot, upper);
                    }
                    x = {clamp01(px1), clamp01(px2)};

                    const double a1 = input + cfg.d1;
                    const double py1 = y[0] + (input - a1 * y[0]) * cfg.dt
                                     + cfg.sigma * y[0] * dW1;
                    const double py2 = y[1] + (cfg.gain * y[0] - cfg.d2 * y[1]) * cfg.dt
                                     + cfg.sigma * y[1] * dW2;
                    y = {py1, py2};
                    if (py1 < 0.0 || py1 > 1.0 || py2 < 0.0 || py2 > 1.0) {
                        linear_exited = true;
                    }

                    save(step + 1);
                }
                if (linear_exited) ++diagnostics.linear_exit_trajectories;
            }

            // Deterministic moment ODE at the same saved times.
            std::vector<MomentState> moments(saved);
            MomentState moment{0.0, 0.0, 0.0, 0.0, 0.0};
            int next_moment_save = 0;
            auto save_moment = [&](int step) {
                while (next_moment_save < saved && save_steps[next_moment_save] == step) {
                    moments[next_moment_save] = moment;
                    ++next_moment_save;
                }
            };
            save_moment(0);
            for (int step = 0; step < steps; ++step) {
                moment = rk4_step(moment, cfg.dt, input, cfg);
                save_moment(step + 1);
            }

            auto mc_stats = [&](const TimeStats& stats, int s, int node) {
                const double mean = stats.sum[s][node] / cfg.replicates;
                const double second = stats.sumsq[s][node] / cfg.replicates;
                const double sd = std::sqrt(std::max(0.0, second - mean * mean));
                return std::array<double, 3>{mean, sd, safe_cv(mean, sd)};
            };

            for (int s = 0; s < saved; ++s) {
                const double time = save_steps[s] * cfg.dt;
                for (int node = 0; node < 2; ++node) {
                    const auto f = mc_stats(full, s, node);
                    const auto l = mc_stats(linear, s, node);
                    const double mm = moments[s][node];
                    const double qq = (node == 0) ? moments[s][2] : moments[s][4];
                    const double msd = std::sqrt(std::max(0.0, qq - mm * mm));
                    const double mcv = safe_cv(mm, msd);
                    const int node_id = node + 1;

                    time_out << input << ',' << time << ",bounded_mc," << node_id << ','
                             << f[0] << ',' << f[1] << ',' << f[2] << '\n';
                    time_out << input << ',' << time << ",linear_mc," << node_id << ','
                             << l[0] << ',' << l[1] << ',' << l[2] << '\n';
                    time_out << input << ',' << time << ",moment_ode," << node_id << ','
                             << mm << ',' << msd << ',' << mcv << '\n';
                }
            }

            const int last = saved - 1;
            for (int node = 0; node < 2; ++node) {
                const auto f = mc_stats(full, last, node);
                const auto l = mc_stats(linear, last, node);
                const double mm = moments[last][node];
                const double qq = (node == 0) ? moments[last][2] : moments[last][4];
                const double msd = std::sqrt(std::max(0.0, qq - mm * mm));
                const double mcv = safe_cv(mm, msd);
                const int node_id = node + 1;

                summary_out << input << ",bounded_mc," << node_id << ','
                            << f[0] << ',' << f[1] << ',' << f[2] << '\n';
                summary_out << input << ",linear_mc," << node_id << ','
                            << l[0] << ',' << l[1] << ',' << l[2] << '\n';
                summary_out << input << ",moment_ode," << node_id << ','
                            << mm << ',' << msd << ',' << mcv << '\n';

                const double lb_abs = std::abs(l[0] - f[0]);
                const double lb_rel = lb_abs / std::max(std::abs(f[0]), 1e-12);
                const double ml_abs = std::abs(mm - l[0]);
                const double ml_rel = ml_abs / std::max(std::abs(l[0]), 1e-12);
                error_out << input << ',' << node_id << ',' << f[0] << ',' << l[0] << ','
                          << mm << ',' << lb_abs << ',' << lb_rel << ',' << ml_abs << ','
                          << ml_rel << ',' << f[2] << ',' << l[2] << ',' << mcv << '\n';
            }

            const double violation_rate = diagnostics.bounded_predictor_values > 0
                ? static_cast<double>(diagnostics.bounded_predictor_violations) /
                  diagnostics.bounded_predictor_values
                : 0.0;
            const double linear_exit_fraction = static_cast<double>(
                diagnostics.linear_exit_trajectories) / cfg.replicates;
            diag_out << input << ',' << cfg.replicates << ',' << cfg.dt << ','
                     << cfg.final_time << ',' << cfg.sigma << ',' << cfg.gain << ','
                     << cfg.d1 << ',' << cfg.d2 << ','
                     << diagnostics.bounded_predictor_violations << ','
                     << diagnostics.bounded_predictor_values << ',' << violation_rate << ','
                     << diagnostics.bounded_max_lower_overshoot << ','
                     << diagnostics.bounded_max_upper_overshoot << ','
                     << diagnostics.linear_exit_trajectories << ','
                     << linear_exit_fraction << '\n';
        }

        std::cout << "Saved low-activation validation data to " << cfg.out_dir << "\n";
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << "\n\n";
        print_help();
        return 1;
    }

    return 0;
}
