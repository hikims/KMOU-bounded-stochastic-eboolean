#include "nfkb_model.hpp"

#include <limits>
#include <sys/stat.h>
#include <sys/types.h>
#ifdef _WIN32
#include <direct.h>
#endif

namespace nfkb {

namespace {

double clamp01(double x) {
    if (x < 0.0) return 0.0;
    if (x > 1.0) return 1.0;
    return x;
}

std::string edge_type_name(EdgeType t) {
    return (t == EdgeType::Activation) ? "activation" : "inhibition";
}

std::string edge_role_name(EdgeRole r) {
    switch (r) {
        case EdgeRole::Core: return "core";
        case EdgeRole::NoncanonicalToCanonical: return "noncanonical_to_canonical";
        case EdgeRole::CanonicalToNoncanonical: return "canonical_to_noncanonical";
        default: return "unknown";
    }
}

bool dir_exists(const std::string& path) {
    if (path.empty()) return true;
#ifdef _WIN32
    struct _stat info;
    if (_stat(path.c_str(), &info) != 0) return false;
    return (info.st_mode & _S_IFDIR) != 0;
#else
    struct stat info;
    if (stat(path.c_str(), &info) != 0) return false;
    return S_ISDIR(info.st_mode);
#endif
}

void make_dir_if_needed(const std::string& path) {
    if (path.empty() || dir_exists(path)) return;
#ifdef _WIN32
    _mkdir(path.c_str());
#else
    mkdir(path.c_str(), 0755);
#endif
}

void ensure_parent_dir(const std::string& filename) {
    const std::size_t pos = filename.find_last_of("/\\");
    if (pos == std::string::npos) return;

    const std::string parent = filename.substr(0, pos);
    if (parent.empty()) return;

    std::size_t start = 0;
#ifdef _WIN32
    if (parent.size() >= 2 && parent[1] == ':') start = 2;
#endif

    for (std::size_t i = start; i <= parent.size(); ++i) {
        const bool at_sep = (i < parent.size() && (parent[i] == '/' || parent[i] == '\\'));
        const bool at_end = (i == parent.size());
        if (!at_sep && !at_end) continue;

        std::string current = parent.substr(0, i);
#ifndef _WIN32
        if (current.empty() && !parent.empty() && parent[0] == '/') current = "/";
#endif
        if (!current.empty() && current.back() != ':') make_dir_if_needed(current);
    }
}

double safe_cv(double mean, double sd) {
    return (std::abs(mean) > 1e-12)
        ? sd / std::abs(mean)
        : std::numeric_limits<double>::quiet_NaN();
}

} // namespace

void StepDiagnostics::observe(double predictor) {
    ++predictor_values;
    const double lower = std::max(0.0, -predictor);
    const double upper = std::max(0.0, predictor - 1.0);
    if (lower > 0.0 || upper > 0.0) ++predictor_violations;
    max_lower_overshoot = std::max(max_lower_overshoot, lower);
    max_upper_overshoot = std::max(max_upper_overshoot, upper);
    total_overshoot += lower + upper;
}

void StepDiagnostics::add(const StepDiagnostics& other) {
    predictor_violations += other.predictor_violations;
    predictor_values += other.predictor_values;
    max_lower_overshoot = std::max(max_lower_overshoot, other.max_lower_overshoot);
    max_upper_overshoot = std::max(max_upper_overshoot, other.max_upper_overshoot);
    total_overshoot += other.total_overshoot;
}

NFkBModel::NFkBModel(SimulationConfig cfg) : cfg_(std::move(cfg)) {
    if (cfg_.alpha_scale <= 0.0) throw std::invalid_argument("alpha_scale must be positive");
    if (cfg_.beta_scale <= 0.0) throw std::invalid_argument("beta_scale must be positive");
    if (cfg_.sigma < 0.0 || cfg_.sigma_scale < 0.0) {
        throw std::invalid_argument("sigma and sigma_scale must be nonnegative");
    }
    if (cfg_.rb < 0.5 || cfg_.rb > 1.0) {
        throw std::invalid_argument("Waller-Kraft r must lie in [0.5, 1]");
    }

    build_nodes();
    build_edges();
    build_compound_targets();
    set_scenario(cfg_.scenario);
    reset_fixed_weights();
}

void NFkBModel::build_nodes() {
    nodes_.resize(NODE_COUNT);
    auto set_node = [&](int id, const std::string& name, bool receptor=false) {
        nodes_[id].name = name;
        nodes_[id].receptor = receptor;
        nodes_[id].alpha = cfg_.alpha_scale;
        nodes_[id].beta = cfg_.beta_scale;
        nodes_[id].sigma = cfg_.sigma * cfg_.sigma_scale;
        nodes_[id].basal_activation = 0.0;
        nodes_[id].basal_deactivation = 0.05;
    };

    set_node(TNFR, "TNFR", true);
    set_node(IL1R, "IL1R", true);
    set_node(CD40, "CD40", true);
    set_node(LTB, "LTB", true);
    set_node(BAFFR, "BAFFR", true);

    set_node(TRAF25, "TRAF2_5");
    set_node(TRAF6, "TRAF6");
    set_node(TRAF23, "TRAF2_3");
    set_node(BTX, "Btx");
    set_node(TAK1_TAB, "TAK1_TAB");
    set_node(IKK_COMPLEX, "NEMO_IKKa_IKKb");
    set_node(IKBA, "IkBa_release_signal");
    set_node(P50_P65, "p50_p65");
    set_node(TNFA, "TNFa");
    set_node(P100, "p100");
    set_node(PHENOTYPE_CAN, "Phenotype_canonical");

    set_node(NIK, "NIK");
    set_node(IKKA, "IKKa_noncanonical");
    set_node(P100_RELB, "p100_RelB");
    set_node(P52, "p52");
    set_node(BAFF, "BAFF");
    set_node(PHENOTYPE_NONCAN, "Phenotype_noncanonical");

    nodes_[PHENOTYPE_CAN].basal_deactivation = 0.02;
    nodes_[PHENOTYPE_NONCAN].basal_deactivation = 0.02;
}

void NFkBModel::build_edges() {
    auto add = [&](int from, int to, EdgeType type, const std::string& label,
                   bool feedback=false, EdgeRole role=EdgeRole::Core) {
        Edge e;
        e.from = from;
        e.to = to;
        e.type = type;
        e.role = role;
        e.weights = {cfg_.fixed_weight};
        e.feedback = feedback;
        e.label = label;
        edges_.push_back(e);
    };

    add(TNFR,  TRAF25,    EdgeType::Activation, "TNFR_to_TRAF2_5");
    add(IL1R,  TRAF6,     EdgeType::Activation, "IL1R_to_TRAF6");
    add(CD40,  TRAF6,     EdgeType::Activation, "CD40_to_TRAF6",
        false, EdgeRole::NoncanonicalToCanonical);

    add(CD40,  TRAF23,    EdgeType::Activation, "CD40_to_TRAF2_3");
    add(LTB,   TRAF23,    EdgeType::Activation, "LTB_to_TRAF2_3");
    add(BAFFR, TRAF23,    EdgeType::Activation, "BAFFR_to_TRAF2_3");

    add(TRAF25, TAK1_TAB,   EdgeType::Activation, "TRAF2_5_to_TAK1_TAB");
    add(TRAF6,  TAK1_TAB,   EdgeType::Activation, "TRAF6_to_TAK1_TAB");
    add(TAK1_TAB, IKK_COMPLEX, EdgeType::Activation, "TAK1_TAB_to_IKK_complex");
    add(TRAF23, BTX,        EdgeType::Activation, "TRAF2_3_to_Btx",
        false, EdgeRole::NoncanonicalToCanonical);
    add(BTX,    IKK_COMPLEX, EdgeType::Activation, "Btx_to_IKK_complex",
        false, EdgeRole::NoncanonicalToCanonical);

    // IKBA is interpreted as an IkBa-degradation/release signal, as in the paper code.
    add(IKK_COMPLEX, IKBA,     EdgeType::Activation, "IKK_complex_to_IkBa_release_signal");
    add(IKBA,        P50_P65,  EdgeType::Activation, "IkBa_release_signal_to_p50_p65");

    add(P50_P65, TNFA,          EdgeType::Activation, "p50_p65_to_TNFa");
    add(P50_P65, P100,          EdgeType::Activation, "p50_p65_to_p100");
    add(P50_P65, PHENOTYPE_CAN, EdgeType::Activation, "p50_p65_to_canonical_phenotype");

    add(TRAF23, NIK,           EdgeType::Activation, "TRAF2_3_to_NIK");
    add(NIK,    IKKA,          EdgeType::Activation, "NIK_to_IKKa");
    add(IKKA,   P100_RELB,     EdgeType::Activation, "IKKa_to_p100_RelB");
    add(P100,   P100_RELB,     EdgeType::Activation, "p100_to_p100_RelB",
        false, EdgeRole::CanonicalToNoncanonical);
    add(P100_RELB, P52,        EdgeType::Activation, "p100_RelB_to_p52");
    add(P52,    PHENOTYPE_NONCAN, EdgeType::Activation, "p52_to_noncanonical_phenotype");
    add(P100_RELB, BAFF,       EdgeType::Activation, "p100_RelB_to_BAFF");

    add(TNFA, TNFR,  EdgeType::Activation, "TNFa_feedback_to_TNFR", true);
    add(BAFF, BAFFR, EdgeType::Activation, "BAFF_feedback_to_BAFFR", true);
}

void NFkBModel::build_compound_targets() {
    compound_targets_.push_back({IKBA, 0.50, 50.0, 1.0,
                                 "Fucoxanthin_inhibits_IkBa_release_signal"});
    compound_targets_.push_back({P50_P65, 0.80, 50.0, 1.0,
                                 "Fucoxanthin_inhibits_p50_p65"});
    compound_targets_.push_back({IKK_COMPLEX, 0.30, 50.0, 1.0,
                                 "Fucoxanthin_weakly_inhibits_IKK_complex"});
}

void NFkBModel::set_scenario(const std::string& scenario) {
    receptor_stimulus_.clear();
    for (int id : {TNFR, IL1R, CD40, LTB, BAFFR}) receptor_stimulus_[id] = 0.0;

    if (scenario == "canonical") {
        receptor_stimulus_[TNFR] = 1.0;
        receptor_stimulus_[IL1R] = 1.0;
    } else if (scenario == "noncanonical") {
        receptor_stimulus_[CD40] = 1.0;
        receptor_stimulus_[LTB] = 1.0;
        receptor_stimulus_[BAFFR] = 1.0;
    } else if (scenario == "combined") {
        receptor_stimulus_[TNFR] = 1.0;
        receptor_stimulus_[IL1R] = 1.0;
        receptor_stimulus_[CD40] = 1.0;
        receptor_stimulus_[LTB] = 1.0;
        receptor_stimulus_[BAFFR] = 1.0;
    } else if (scenario == "none") {
        // all receptor stimuli remain zero
    } else {
        throw std::invalid_argument("Unknown scenario: " + scenario +
                                    " (use canonical, noncanonical, combined, or none)");
    }
}

void NFkBModel::reset_fixed_weights() {
    for (auto& e : edges_) {
        e.weights.assign(e.weights.size(), cfg_.fixed_weight);
    }
}

WeightRealization NFkBModel::draw_weight_realization(std::mt19937_64& rng) const {
    WeightRealization realization;
    realization.reserve(edges_.size());

    std::uniform_real_distribution<double> unif(0.0, 1.0);
    for (const auto& e : edges_) {
        std::vector<double> weights(e.weights.size(), cfg_.fixed_weight);
        if (cfg_.random_weights) {
            for (double& w : weights) w = unif(rng);
        }
        realization.push_back(std::move(weights));
    }
    return realization;
}

void NFkBModel::set_weight_realization(const WeightRealization& weights) {
    if (weights.size() != edges_.size()) {
        throw std::invalid_argument("weight realization has the wrong number of edges");
    }
    for (std::size_t e = 0; e < edges_.size(); ++e) {
        if (weights[e].size() != edges_[e].weights.size()) {
            throw std::invalid_argument("weight realization has an incompatible edge dimension");
        }
        edges_[e].weights = weights[e];
    }
}

double NFkBModel::wk(const std::vector<double>& weights) const {
    if (weights.empty()) return 0.0;
    auto mm = std::minmax_element(weights.begin(), weights.end());
    return cfg_.rb * (*mm.first) + (1.0 - cfg_.rb) * (*mm.second);
}

double NFkBModel::bounded_response(double z) const {
    if (!cfg_.saturating_inputs) return std::max(0.0, z);
    z = std::max(0.0, z);
    return z / (1.0 + z);
}

double NFkBModel::compound_fraction(double C, double ic50, double hill) const {
    if (C <= 0.0) return 0.0;
    if (ic50 <= 0.0) return 1.0;
    const double num = std::pow(C, hill);
    const double den = std::pow(ic50, hill) + num;
    return (den > 0.0) ? num / den : 0.0;
}

double NFkBModel::edge_scale(const Edge& e) const {
    if (!cfg_.isolate_branches) return 1.0;

    if (cfg_.scenario == "canonical" &&
        e.role == EdgeRole::CanonicalToNoncanonical) {
        return cfg_.cross_branch_scale;
    }
    if (cfg_.scenario == "noncanonical" &&
        e.role == EdgeRole::NoncanonicalToCanonical) {
        return cfg_.cross_branch_scale;
    }
    return 1.0;
}

void NFkBModel::compute_inputs(const std::vector<double>& X,
                               std::vector<double>& A,
                               std::vector<double>& D) const {
    const int N = static_cast<int>(nodes_.size());
    A.assign(N, 0.0);
    D.assign(N, 0.0);

    for (int i = 0; i < N; ++i) {
        A[i] += nodes_[i].basal_activation;
        D[i] += nodes_[i].basal_deactivation;
    }

    if (!cfg_.clamp_receptors) {
        for (const auto& kv : receptor_stimulus_) A[kv.first] += kv.second;
    }

    std::vector<std::vector<double>> outgoing_active_weights(N);
    for (const auto& e : edges_) {
        if (e.feedback && !cfg_.feedback_on) continue;

        const double scale = edge_scale(e);
        if (scale <= 0.0) continue;

        const double b = scale * wk(e.weights);
        if (e.type == EdgeType::Activation) {
            A[e.to] += b * X[e.from];
            outgoing_active_weights[e.from].push_back(b);
        } else {
            D[e.to] += b * X[e.from];
        }
    }

    for (int i = 0; i < N; ++i) {
        if (!outgoing_active_weights[i].empty()) {
            D[i] += wk(outgoing_active_weights[i]);
        }
    }

    for (const auto& target : compound_targets_) {
        const double fraction = compound_fraction(cfg_.fucoxanthin, target.ic50, target.hill);
        D[target.to] += target.rho * fraction;
    }

    for (int i = 0; i < N; ++i) {
        A[i] = bounded_response(A[i]);
        D[i] = bounded_response(D[i]);
    }
}

std::vector<double> NFkBModel::initial_state() const {
    std::vector<double> X(nodes_.size(), 0.0);
    if (cfg_.clamp_receptors) {
        for (const auto& kv : receptor_stimulus_) X[kv.first] = kv.second;
    }
    return X;
}

std::vector<double> NFkBModel::step_with_dW(const std::vector<double>& X,
                                            const std::vector<double>& dW,
                                            StepDiagnostics& diagnostics) const {
    if (X.size() != nodes_.size() || dW.size() != nodes_.size()) {
        throw std::invalid_argument("state and Brownian increment vectors must match node count");
    }

    std::vector<double> A, D;
    compute_inputs(X, A, D);
    std::vector<double> Xnew(X.size(), 0.0);

    for (std::size_t i = 0; i < X.size(); ++i) {
        const double xi = X[i];
        const double drift = nodes_[i].alpha * (1.0 - xi) * A[i]
                           - nodes_[i].beta * xi * D[i];
        const double diffusion = nodes_[i].sigma * xi * (1.0 - xi) * dW[i];
        const double predictor = xi + drift * cfg_.dt + diffusion;

        diagnostics.observe(predictor);
        Xnew[i] = cfg_.projected_em ? clamp01(predictor) : predictor;
    }

    if (cfg_.clamp_receptors) {
        for (const auto& kv : receptor_stimulus_) Xnew[kv.first] = kv.second;
    }

    return Xnew;
}

SimulationResult NFkBModel::run() {
    if (cfg_.steps <= 0) throw std::invalid_argument("steps must be positive");
    if (cfg_.replicates <= 0) throw std::invalid_argument("replicates must be positive");
    if (cfg_.dt <= 0.0) throw std::invalid_argument("dt must be positive");
    if (cfg_.save_every <= 0) throw std::invalid_argument("save_every must be positive");

    const int N = static_cast<int>(nodes_.size());

    std::vector<int> save_steps;
    if (cfg_.write_timecourse) {
        save_steps.reserve(static_cast<std::size_t>(cfg_.steps / cfg_.save_every + 2));
        for (int step_idx = 0; step_idx <= cfg_.steps; step_idx += cfg_.save_every) {
            save_steps.push_back(step_idx);
        }
        if (save_steps.empty() || save_steps.back() != cfg_.steps) save_steps.push_back(cfg_.steps);
    } else {
        save_steps = {0, cfg_.steps};
    }

    const int S = static_cast<int>(save_steps.size());

    SimulationResult result;
    result.save_every = cfg_.save_every;
    result.simulated_steps = cfg_.steps;
    result.replicates = cfg_.replicates;
    result.dt = cfg_.dt;
    result.node_names.reserve(N);
    for (const auto& n : nodes_) result.node_names.push_back(n.name);
    result.time.resize(S);
    result.mean.assign(S, std::vector<double>(N, 0.0));
    result.sd.assign(S, std::vector<double>(N, 0.0));
    result.rep0.assign(S, std::vector<double>(N, 0.0));
    result.final_mean.assign(N, 0.0);
    result.final_sd.assign(N, 0.0);
    result.final_cv.assign(N, 0.0);
    result.auc_mean.assign(N, 0.0);
    result.auc_sd.assign(N, 0.0);
    if (cfg_.write_replicates) result.replicate_metrics.reserve(cfg_.replicates);

    for (int sidx = 0; sidx < S; ++sidx) {
        result.time[sidx] = static_cast<double>(save_steps[sidx]) * cfg_.dt;
    }

    std::vector<std::vector<double>> sum(S, std::vector<double>(N, 0.0));
    std::vector<std::vector<double>> sumsq(S, std::vector<double>(N, 0.0));
    std::vector<double> auc_sum(N, 0.0), auc_sumsq(N, 0.0);
    std::vector<double> final_sum(N, 0.0), final_sumsq(N, 0.0);

    std::mt19937_64 master_rng(cfg_.seed);
    for (int rep = 0; rep < cfg_.replicates; ++rep) {
        const std::uint64_t stream_seed = master_rng();
        std::mt19937_64 rng(stream_seed);
        set_weight_realization(draw_weight_realization(rng));

        std::normal_distribution<double> normal(0.0, 1.0);
        std::vector<double> X = initial_state();
        std::vector<double> auc(N, 0.0);
        StepDiagnostics rep_diagnostics;

        int next_save_index = 0;
        auto save_snapshot = [&](int step_idx) {
            while (next_save_index < S && save_steps[next_save_index] == step_idx) {
                for (int i = 0; i < N; ++i) {
                    sum[next_save_index][i] += X[i];
                    sumsq[next_save_index][i] += X[i] * X[i];
                    if (rep == 0) result.rep0[next_save_index][i] = X[i];
                }
                ++next_save_index;
            }
        };

        save_snapshot(0);
        std::vector<double> dW(N, 0.0);
        for (int step_idx = 0; step_idx < cfg_.steps; ++step_idx) {
            for (int i = 0; i < N; ++i) auc[i] += X[i] * cfg_.dt;
            for (int i = 0; i < N; ++i) dW[i] = std::sqrt(cfg_.dt) * normal(rng);
            X = step_with_dW(X, dW, rep_diagnostics);
            save_snapshot(step_idx + 1);
        }

        result.diagnostics.add(rep_diagnostics);
        for (int i = 0; i < N; ++i) {
            auc_sum[i] += auc[i];
            auc_sumsq[i] += auc[i] * auc[i];
            final_sum[i] += X[i];
            final_sumsq[i] += X[i] * X[i];
        }

        if (cfg_.write_replicates) {
            result.replicate_metrics.push_back({rep, stream_seed, X, auc});
        }

        if (cfg_.progress_every > 0 &&
            ((rep + 1) % cfg_.progress_every == 0 || rep + 1 == cfg_.replicates)) {
            std::cout << "Progress: " << (rep + 1) << " / " << cfg_.replicates
                      << " replicates complete" << std::endl;
        }
    }

    const double R = static_cast<double>(cfg_.replicates);
    for (int sidx = 0; sidx < S; ++sidx) {
        for (int i = 0; i < N; ++i) {
            const double mean = sum[sidx][i] / R;
            const double variance = std::max(0.0, sumsq[sidx][i] / R - mean * mean);
            result.mean[sidx][i] = mean;
            result.sd[sidx][i] = std::sqrt(variance);
        }
    }

    for (int i = 0; i < N; ++i) {
        result.final_mean[i] = final_sum[i] / R;
        const double final_var = std::max(
            0.0, final_sumsq[i] / R - result.final_mean[i] * result.final_mean[i]);
        result.final_sd[i] = std::sqrt(final_var);
        result.final_cv[i] = safe_cv(result.final_mean[i], result.final_sd[i]);

        result.auc_mean[i] = auc_sum[i] / R;
        const double auc_var = std::max(
            0.0, auc_sumsq[i] / R - result.auc_mean[i] * result.auc_mean[i]);
        result.auc_sd[i] = std::sqrt(auc_var);
    }

    return result;
}

void write_timecourse_csv(const std::string& filename,
                          const std::vector<std::string>& node_names,
                          const std::vector<double>& time,
                          const std::vector<std::vector<double>>& data) {
    ensure_parent_dir(filename);
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("Cannot open " + filename);

    out << "time";
    for (const auto& name : node_names) out << ',' << name;
    out << '\n';
    out << std::setprecision(10);
    for (std::size_t t = 0; t < time.size(); ++t) {
        out << time[t];
        for (double value : data[t]) out << ',' << value;
        out << '\n';
    }
}

void write_summary_csv(const std::string& filename,
                       const SimulationResult& result) {
    ensure_parent_dir(filename);
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("Cannot open " + filename);

    out << "node,final_mean,final_sd,final_se,final_cv,auc_mean,auc_sd,auc_se,auc_cv\n";
    out << std::setprecision(10);
    const double root_n = std::sqrt(static_cast<double>(result.replicates));
    for (std::size_t i = 0; i < result.node_names.size(); ++i) {
        out << result.node_names[i] << ','
            << result.final_mean[i] << ','
            << result.final_sd[i] << ','
            << result.final_sd[i] / root_n << ','
            << result.final_cv[i] << ','
            << result.auc_mean[i] << ','
            << result.auc_sd[i] << ','
            << result.auc_sd[i] / root_n << ','
            << safe_cv(result.auc_mean[i], result.auc_sd[i]) << '\n';
    }
}

void write_replicates_csv(const std::string& filename,
                          const SimulationResult& result) {
    ensure_parent_dir(filename);
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("Cannot open " + filename);

    out << "replicate,stream_seed,node,final_value,auc\n";
    out << std::setprecision(10);
    for (const auto& rep : result.replicate_metrics) {
        for (std::size_t i = 0; i < result.node_names.size(); ++i) {
            out << rep.replicate << ',' << rep.stream_seed << ','
                << result.node_names[i] << ','
                << rep.final_value[i] << ',' << rep.auc[i] << '\n';
        }
    }
}

void write_diagnostics_csv(const std::string& filename,
                           const SimulationConfig& cfg,
                           const SimulationResult& result) {
    ensure_parent_dir(filename);
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("Cannot open " + filename);

    out << "scenario,fucoxanthin,steps,dt,final_time,replicates,save_every,seed,"
        << "random_weights,fixed_weight,r,alpha_scale,beta_scale,sigma_base,sigma_scale,"
        << "sigma_effective,feedback_on,clamp_receptors,projected_em,saturating_inputs,"
        << "isolate_branches,cross_branch_scale,predictor_violations,predictor_values,"
        << "violation_rate,max_lower_overshoot,max_upper_overshoot,mean_overshoot_per_value\n";

    const double violation_rate = result.diagnostics.predictor_values > 0
        ? static_cast<double>(result.diagnostics.predictor_violations) /
          static_cast<double>(result.diagnostics.predictor_values)
        : 0.0;
    const double mean_overshoot = result.diagnostics.predictor_values > 0
        ? result.diagnostics.total_overshoot /
          static_cast<double>(result.diagnostics.predictor_values)
        : 0.0;

    out << std::setprecision(12)
        << cfg.scenario << ',' << cfg.fucoxanthin << ',' << cfg.steps << ',' << cfg.dt << ','
        << cfg.steps * cfg.dt << ',' << cfg.replicates << ',' << cfg.save_every << ',' << cfg.seed << ','
        << (cfg.random_weights ? 1 : 0) << ',' << cfg.fixed_weight << ',' << cfg.rb << ','
        << cfg.alpha_scale << ',' << cfg.beta_scale << ',' << cfg.sigma << ',' << cfg.sigma_scale << ','
        << cfg.sigma * cfg.sigma_scale << ',' << (cfg.feedback_on ? 1 : 0) << ','
        << (cfg.clamp_receptors ? 1 : 0) << ',' << (cfg.projected_em ? 1 : 0) << ','
        << (cfg.saturating_inputs ? 1 : 0) << ',' << (cfg.isolate_branches ? 1 : 0) << ','
        << cfg.cross_branch_scale << ',' << result.diagnostics.predictor_violations << ','
        << result.diagnostics.predictor_values << ',' << violation_rate << ','
        << result.diagnostics.max_lower_overshoot << ','
        << result.diagnostics.max_upper_overshoot << ',' << mean_overshoot << '\n';
}

void write_network_csv(const std::string& filename,
                       const std::vector<Node>& nodes,
                       const std::vector<Edge>& edges) {
    ensure_parent_dir(filename);
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("Cannot open " + filename);

    out << "from,to,type,role,feedback,label\n";
    for (const auto& edge : edges) {
        out << nodes[edge.from].name << ',' << nodes[edge.to].name << ','
            << edge_type_name(edge.type) << ',' << edge_role_name(edge.role) << ','
            << (edge.feedback ? 1 : 0) << ',' << edge.label << '\n';
    }
}

} // namespace nfkb
