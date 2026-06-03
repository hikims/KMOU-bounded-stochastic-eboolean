#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>
#include <filesystem>

namespace fs = std::filesystem;

struct Options {
    std::string out_prefix = "results/theory";
    int threads = std::max(1u, std::thread::hardware_concurrency());
    std::size_t boundary_mc = 5000;
    std::size_t moment_mc = 20000;
    double boundary_dt = 0.001;
    double boundary_T = 0.8;
    int boundary_save_every = 10;
    double moment_dt = 0.002;
    double moment_T = 8.0;
    int moment_save_every = 20;
};

static void usage(const char* prog) {
    std::cerr << "Usage: " << prog << " [options]\n"
              << "  --out-prefix PATH          default results/theory\n"
              << "  --threads N                default hardware_concurrency\n"
              << "  --boundary-mc N            default 5000\n"
              << "  --moment-mc N              default 20000\n"
              << "  --boundary-dt X            default 0.001\n"
              << "  --boundary-T X             default 0.8\n"
              << "  --moment-dt X              default 0.002\n"
              << "  --moment-T X               default 8.0\n";
}

static Options parse_args(int argc, char** argv) {
    Options opt;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto need = [&](const std::string& name) -> std::string {
            if (i + 1 >= argc) throw std::runtime_error("missing value for " + name);
            return argv[++i];
        };
        if (a == "--out-prefix") opt.out_prefix = need(a);
        else if (a == "--threads") opt.threads = std::max(1, std::stoi(need(a)));
        else if (a == "--boundary-mc") opt.boundary_mc = static_cast<std::size_t>(std::stoull(need(a)));
        else if (a == "--moment-mc") opt.moment_mc = static_cast<std::size_t>(std::stoull(need(a)));
        else if (a == "--boundary-dt") opt.boundary_dt = std::stod(need(a));
        else if (a == "--boundary-T") opt.boundary_T = std::stod(need(a));
        else if (a == "--moment-dt") opt.moment_dt = std::stod(need(a));
        else if (a == "--moment-T") opt.moment_T = std::stod(need(a));
        else if (a == "--help" || a == "-h") { usage(argv[0]); std::exit(0); }
        else throw std::runtime_error("unknown option: " + a);
    }
    return opt;
}

static void ensure_parent(const std::string& path) {
    fs::path p(path);
    if (!p.parent_path().empty()) fs::create_directories(p.parent_path());
}

static std::ofstream open_csv(const std::string& path) {
    ensure_parent(path);
    std::ofstream f(path);
    if (!f) throw std::runtime_error("cannot open output file: " + path);
    f << std::setprecision(10);
    return f;
}

static std::string csv_escape(const std::string& s) {
    if (s.find_first_of(",\"") == std::string::npos) return s;
    std::string out = "\"";
    for (char c : s) out += (c == '"' ? "\"\"" : std::string(1, c));
    out += "\"";
    return out;
}

static void parallel_ranges(std::size_t n, int threads, const std::function<void(std::size_t, std::size_t, int)>& fn) {
    if (n == 0) return;
    int th = std::max(1, std::min<int>(threads, static_cast<int>(n)));
    std::vector<std::thread> workers;
    workers.reserve(th);
    std::size_t block = (n + th - 1) / th;
    for (int tid = 0; tid < th; ++tid) {
        std::size_t begin = static_cast<std::size_t>(tid) * block;
        std::size_t end = std::min(n, begin + block);
        if (begin >= end) continue;
        workers.emplace_back([=, &fn]() { fn(begin, end, tid); });
    }
    for (auto& w : workers) w.join();
}

struct BoundaryResult {
    std::string name;
    std::vector<double> time, minv, mean, maxv;
    double vmax = 0.0;
    double exit_fraction = 0.0;
};

struct BoundaryLocal {
    std::vector<double> sum, minv, maxv;
    double vmax = 0.0;
    std::size_t exit_count = 0;
    explicit BoundaryLocal(std::size_t nsave)
        : sum(nsave, 0.0), minv(nsave, std::numeric_limits<double>::infinity()),
          maxv(nsave, -std::numeric_limits<double>::infinity()) {}
};

enum class BoundaryModel { Original, Tanh, RevisedProjected };

static BoundaryResult simulate_boundary_model(BoundaryModel model, const std::string& name,
                                              const Options& opt, unsigned seed_base) {
    const std::size_t mc = opt.boundary_mc;
    const double dt = opt.boundary_dt;
    const int steps = static_cast<int>(std::llround(opt.boundary_T / dt));
    const int save_every = opt.boundary_save_every;
    const int nsave = steps / save_every + 1;
    const double sigma = 0.25;
    const double sqrt_dt = std::sqrt(dt);

    std::vector<BoundaryLocal> locals;
    locals.reserve(opt.threads);
    for (int i = 0; i < opt.threads; ++i) locals.emplace_back(nsave);

    parallel_ranges(mc, opt.threads, [&](std::size_t begin, std::size_t end, int tid) {
        auto& loc = locals[tid];
        for (std::size_t rep = begin; rep < end; ++rep) {
            std::mt19937_64 rng(seed_base + 1000003ULL * (rep + 1));
            std::normal_distribution<double> normal(0.0, 1.0);
            double x = 0.85;
            bool exit_any = false;

            for (int n = 0; n <= steps; ++n) {
                if (n % save_every == 0) {
                    std::size_t s = static_cast<std::size_t>(n / save_every);
                    loc.sum[s] += x;
                    loc.minv[s] = std::min(loc.minv[s], x);
                    loc.maxv[s] = std::max(loc.maxv[s], x);
                }
                if (n == steps) break;

                double dW = sqrt_dt * normal(rng);
                if (model == BoundaryModel::Original) {
                    double drift = 0.65 - 0.15 * x;
                    double diffusion = sigma * x;
                    x += drift * dt + diffusion * dW;
                } else if (model == BoundaryModel::Tanh) {
                    double drift = std::tanh(0.65 - 0.15 * x);
                    double diffusion = sigma * x;
                    x += drift * dt + diffusion * dW;
                } else {
                    double drift = 0.65 * (1.0 - x) - 0.15 * x;
                    double diffusion = sigma * x * (1.0 - x);
                    x += drift * dt + diffusion * dW;
                    x = std::max(0.0, std::min(1.0, x));
                }

                bool outside = (x < 0.0 || x > 1.0);
                exit_any = exit_any || outside;
                double upper = x - 1.0;
                double lower = -x;
                loc.vmax = std::max(loc.vmax, std::max({upper, lower, 0.0}));
            }
            if (exit_any) loc.exit_count++;
        }
    });

    BoundaryResult out;
    out.name = name;
    out.time.resize(nsave);
    out.minv.assign(nsave, std::numeric_limits<double>::infinity());
    out.maxv.assign(nsave, -std::numeric_limits<double>::infinity());
    out.mean.assign(nsave, 0.0);

    std::size_t total_exit = 0;
    for (std::size_t s = 0; s < static_cast<std::size_t>(nsave); ++s) {
        out.time[s] = static_cast<double>(s * save_every) * dt;
        double sum = 0.0;
        for (const auto& loc : locals) {
            sum += loc.sum[s];
            out.minv[s] = std::min(out.minv[s], loc.minv[s]);
            out.maxv[s] = std::max(out.maxv[s], loc.maxv[s]);
        }
        out.mean[s] = sum / static_cast<double>(mc);
    }
    for (const auto& loc : locals) {
        out.vmax = std::max(out.vmax, loc.vmax);
        total_exit += loc.exit_count;
    }
    out.exit_fraction = static_cast<double>(total_exit) / static_cast<double>(mc);
    return out;
}

static void run_boundary(const Options& opt) {
    std::cout << "[boundary] mc=" << opt.boundary_mc << ", threads=" << opt.threads << "\n";
    std::vector<BoundaryResult> res;
    res.push_back(simulate_boundary_model(BoundaryModel::Original, "Original EM", opt, 11));
    res.push_back(simulate_boundary_model(BoundaryModel::Tanh, "Tanh EM", opt, 12));
    res.push_back(simulate_boundary_model(BoundaryModel::RevisedProjected, "Revised projected EM", opt, 13));

    {
        auto f = open_csv(opt.out_prefix + "_boundary_timeseries.csv");
        f << "model,time,min,mean,max\n";
        for (const auto& r : res) {
            for (std::size_t i = 0; i < r.time.size(); ++i) {
                f << csv_escape(r.name) << ',' << r.time[i] << ',' << r.minv[i] << ',' << r.mean[i] << ',' << r.maxv[i] << '\n';
            }
        }
    }
    {
        auto f = open_csv(opt.out_prefix + "_boundary_summary.csv");
        f << "model,Vmax,exit_fraction\n";
        for (const auto& r : res) {
            f << csv_escape(r.name) << ',' << r.vmax << ',' << r.exit_fraction << '\n';
        }
    }
}

struct MomentSeries {
    std::vector<double> time, mc_m1, mc_m2, mc_cv1, mc_cv2;
    std::vector<double> ode_m1, ode_m2, ode_cv1, ode_cv2;
};

struct MomentLocal {
    std::vector<double> sum1, sum2, sq1, sq2;
    explicit MomentLocal(std::size_t nsave)
        : sum1(nsave, 0.0), sum2(nsave, 0.0), sq1(nsave, 0.0), sq2(nsave, 0.0) {}
};

static void simulate_moment_mc(const Options& opt, MomentSeries& out, double a, double g, double u, double sigma) {
    const std::size_t mc = opt.moment_mc;
    const double dt = opt.moment_dt;
    const int steps = static_cast<int>(std::llround(opt.moment_T / dt));
    const int save_every = opt.moment_save_every;
    const int nsave = steps / save_every + 1;
    const double sqrt_dt = std::sqrt(dt);

    std::vector<MomentLocal> locals;
    locals.reserve(opt.threads);
    for (int i = 0; i < opt.threads; ++i) locals.emplace_back(nsave);

    parallel_ranges(mc, opt.threads, [&](std::size_t begin, std::size_t end, int tid) {
        auto& loc = locals[tid];
        for (std::size_t rep = begin; rep < end; ++rep) {
            std::mt19937_64 rng(21 + 1000003ULL * (rep + 1));
            std::normal_distribution<double> normal(0.0, 1.0);
            double y1 = 0.01;
            double y2 = 0.01;
            for (int n = 0; n <= steps; ++n) {
                if (n % save_every == 0) {
                    std::size_t s = static_cast<std::size_t>(n / save_every);
                    loc.sum1[s] += y1;
                    loc.sum2[s] += y2;
                    loc.sq1[s] += y1 * y1;
                    loc.sq2[s] += y2 * y2;
                }
                if (n == steps) break;
                double old_y1 = y1;
                double dW1 = sqrt_dt * normal(rng);
                double dW2 = sqrt_dt * normal(rng);
                y1 += (u - a * y1) * dt + sigma * y1 * dW1;
                y2 += (g * old_y1 - a * y2) * dt + sigma * y2 * dW2;
                y1 = std::max(y1, 1e-12);
                y2 = std::max(y2, 1e-12);
            }
        }
    });

    out.time.resize(nsave);
    out.mc_m1.resize(nsave); out.mc_m2.resize(nsave); out.mc_cv1.resize(nsave); out.mc_cv2.resize(nsave);
    for (std::size_t s = 0; s < static_cast<std::size_t>(nsave); ++s) {
        double sum1 = 0.0, sum2 = 0.0, sq1 = 0.0, sq2 = 0.0;
        for (const auto& loc : locals) {
            sum1 += loc.sum1[s]; sum2 += loc.sum2[s]; sq1 += loc.sq1[s]; sq2 += loc.sq2[s];
        }
        double m1 = sum1 / static_cast<double>(mc);
        double m2 = sum2 / static_cast<double>(mc);
        double v1 = std::max(0.0, sq1 / static_cast<double>(mc) - m1 * m1);
        double v2 = std::max(0.0, sq2 / static_cast<double>(mc) - m2 * m2);
        out.time[s] = static_cast<double>(s * save_every) * dt;
        out.mc_m1[s] = m1;
        out.mc_m2[s] = m2;
        out.mc_cv1[s] = std::sqrt(v1) / m1;
        out.mc_cv2[s] = std::sqrt(v2) / m2;
    }
}

static void integrate_moment_ode(const Options& opt, MomentSeries& out, double a, double g, double u, double sigma) {
    const double dt = opt.moment_dt;
    const int steps = static_cast<int>(std::llround(opt.moment_T / dt));
    const int save_every = opt.moment_save_every;
    const int nsave = steps / save_every + 1;

    out.ode_m1.resize(nsave); out.ode_m2.resize(nsave); out.ode_cv1.resize(nsave); out.ode_cv2.resize(nsave);
    double m1 = 0.01, m2 = 0.01;
    double q11 = m1 * m1, q12 = m1 * m2, q22 = m2 * m2;
    const double s2 = sigma * sigma;

    for (int n = 0; n <= steps; ++n) {
        if (n % save_every == 0) {
            std::size_t s = static_cast<std::size_t>(n / save_every);
            double v1 = std::max(0.0, q11 - m1 * m1);
            double v2 = std::max(0.0, q22 - m2 * m2);
            out.ode_m1[s] = m1;
            out.ode_m2[s] = m2;
            out.ode_cv1[s] = std::sqrt(v1) / m1;
            out.ode_cv2[s] = std::sqrt(v2) / m2;
        }
        if (n == steps) break;
        // Moment ODE for dY1=(u-aY1)dt+sY1dW1, dY2=(gY1-aY2)dt+sY2dW2
        double dm1 = u - a * m1;
        double dm2 = g * m1 - a * m2;
        double dq11 = -2.0 * a * q11 + 2.0 * u * m1 + s2 * q11;
        double dq12 = g * q11 - 2.0 * a * q12 + u * m2;
        double dq22 = 2.0 * g * q12 - 2.0 * a * q22 + s2 * q22;
        m1 += dt * dm1;
        m2 += dt * dm2;
        q11 += dt * dq11;
        q12 += dt * dq12;
        q22 += dt * dq22;
    }
}

static void run_moment(const Options& opt) {
    std::cout << "[moment] mc=" << opt.moment_mc << ", threads=" << opt.threads << "\n";
    const double a = 1.0, g = 1.0, u = 1.0, sigma = 0.7;
    MomentSeries ser;
    simulate_moment_mc(opt, ser, a, g, u, sigma);
    integrate_moment_ode(opt, ser, a, g, u, sigma);

    {
        auto f = open_csv(opt.out_prefix + "_moment_timeseries.csv");
        f << "time,mc_m1,mc_m2,mc_cv1,mc_cv2,ode_m1,ode_m2,ode_cv1,ode_cv2\n";
        for (std::size_t i = 0; i < ser.time.size(); ++i) {
            f << ser.time[i] << ',' << ser.mc_m1[i] << ',' << ser.mc_m2[i] << ',' << ser.mc_cv1[i] << ',' << ser.mc_cv2[i]
              << ',' << ser.ode_m1[i] << ',' << ser.ode_m2[i] << ',' << ser.ode_cv1[i] << ',' << ser.ode_cv2[i] << '\n';
        }
    }
    double sigma2 = sigma * sigma;
    double cv1 = std::sqrt(sigma2 / (2.0 * a - sigma2));
    double cv2 = std::sqrt(sigma2 * (3.0 * a - sigma2) / std::pow(2.0 * a - sigma2, 2));
    double theory_ratio = cv2 / cv1;
    double mc_ratio = ser.mc_cv2.back() / ser.mc_cv1.back();
    double ode_ratio = ser.ode_cv2.back() / ser.ode_cv1.back();
    {
        auto f = open_csv(opt.out_prefix + "_moment_summary.csv");
        f << "a,g,u,sigma,theory_CV1,theory_CV2,theory_CV2_over_CV1,MC_final_CV2_over_CV1,moment_ODE_final_CV2_over_CV1\n";
        f << a << ',' << g << ',' << u << ',' << sigma << ',' << cv1 << ',' << cv2 << ',' << theory_ratio << ',' << mc_ratio << ',' << ode_ratio << '\n';
    }
}

struct HopfTrace { double L; std::vector<double> time, x1; };

static double hopf_amplitude(double L, bool save, HopfTrace* trace) {
    const double eta = 10.0, T = 220.0, dt = 0.02;
    int steps = static_cast<int>(std::llround(T / dt));
    double z1 = 0.05, z2 = 0.02, z3 = -0.01;
    std::vector<double> last;
    if (save && trace) { trace->L = L; trace->time.clear(); trace->x1.clear(); }

    auto deriv = [&](double y1, double y2, double y3) {
        double r2 = y1*y1 + y2*y2 + y3*y3;
        std::array<double,3> d{
            -y1 - L*y3 - eta*r2*y1,
            y1 - y2 - eta*r2*y2,
            y2 - y3 - eta*r2*y3
        };
        return d;
    };

    for (int n = 0; n <= steps; ++n) {
        double t = n * dt;
        if (save && trace && t >= 120.0) {
            trace->time.push_back(t - 120.0);
            trace->x1.push_back(0.5 + z1);
        }
        if (n > steps/2) last.push_back(z1);
        if (n == steps) break;
        auto k1 = deriv(z1,z2,z3);
        auto k2 = deriv(z1+0.5*dt*k1[0], z2+0.5*dt*k1[1], z3+0.5*dt*k1[2]);
        auto k3 = deriv(z1+0.5*dt*k2[0], z2+0.5*dt*k2[1], z3+0.5*dt*k2[2]);
        auto k4 = deriv(z1+dt*k3[0], z2+dt*k3[1], z3+dt*k3[2]);
        z1 += dt/6.0*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
        z2 += dt/6.0*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
        z3 += dt/6.0*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2]);
    }
    auto [mn, mx] = std::minmax_element(last.begin(), last.end());
    return 0.5 * (*mx - *mn);
}

static void run_hopf(const Options& opt) {
    (void)opt;
    std::cout << "[hopf] deterministic amplitude scan\n";
    const double Lc = 8.0;
    std::vector<double> Ls, amps;
    for (int i = 0; i <= 32; ++i) {
        double L = 6.5 + (10.5 - 6.5) * static_cast<double>(i) / 32.0;
        Ls.push_back(L);
        amps.push_back(hopf_amplitude(L, false, nullptr));
    }
    {
        auto f = open_csv(opt.out_prefix + "_hopf_amplitude.csv");
        f << "L,L_over_Lc,amplitude_z1\n";
        for (std::size_t i = 0; i < Ls.size(); ++i) f << Ls[i] << ',' << Ls[i]/Lc << ',' << amps[i] << '\n';
    }
    std::vector<double> trace_Ls{7.0, 8.2, 10.0};
    {
        auto f = open_csv(opt.out_prefix + "_hopf_timeseries.csv");
        f << "L,time_after_transient,X1\n";
        for (double L : trace_Ls) {
            HopfTrace tr;
            hopf_amplitude(L, true, &tr);
            for (std::size_t i = 0; i < tr.time.size(); ++i) {
                f << L << ',' << tr.time[i] << ',' << tr.x1[i] << '\n';
            }
        }
    }
}

int main(int argc, char** argv) {
    try {
        Options opt = parse_args(argc, argv);
        std::cout << "Theory benchmark simulation\n";
        std::cout << "out_prefix=" << opt.out_prefix << ", threads=" << opt.threads << "\n";
        run_boundary(opt);
        run_moment(opt);
        run_hopf(opt);
        std::cout << "Done. CSV outputs written with prefix " << opt.out_prefix << "\n";
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        usage(argv[0]);
        return 1;
    }
    return 0;
}
