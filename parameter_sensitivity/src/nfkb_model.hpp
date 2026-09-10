#ifndef REVISION_NFKB_MODEL_HPP
#define REVISION_NFKB_MODEL_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nfkb {

// Node indexing is identical to the NF-kB implementation used for the paper.
enum NodeId : int {
    TNFR = 0,
    IL1R,
    CD40,
    LTB,
    BAFFR,
    TRAF25,
    TRAF6,
    TRAF23,
    BTX,
    TAK1_TAB,
    IKK_COMPLEX,
    IKBA,
    P50_P65,
    TNFA,
    P100,
    PHENOTYPE_CAN,
    NIK,
    IKKA,
    P100_RELB,
    P52,
    BAFF,
    PHENOTYPE_NONCAN,
    NODE_COUNT
};

enum class EdgeType {
    Activation,
    Inhibition
};

enum class EdgeRole {
    Core,
    NoncanonicalToCanonical,
    CanonicalToNoncanonical
};

struct Node {
    std::string name;
    double alpha = 1.0;
    double beta = 1.0;
    double sigma = 0.01;
    double basal_activation = 0.0;
    double basal_deactivation = 0.05;
    bool receptor = false;
};

struct Edge {
    int from = -1;
    int to = -1;
    EdgeType type = EdgeType::Activation;
    EdgeRole role = EdgeRole::Core;
    std::vector<double> weights;
    bool feedback = false;
    std::string label;
};

struct CompoundTarget {
    int to = -1;
    double rho = 0.0;
    double ic50 = 50.0;
    double hill = 1.0;
    std::string label;
};

struct SimulationConfig {
    int steps = 1000;
    int replicates = 1;
    int save_every = 1;
    int progress_every = 0;
    double dt = 0.01;

    // Baseline model parameters and global sensitivity multipliers.
    double rb = 0.75;
    double fixed_weight = 0.5;
    double sigma = 0.01;
    double alpha_scale = 1.0;
    double beta_scale = 1.0;
    double sigma_scale = 1.0;

    bool random_weights = false;
    bool feedback_on = false;
    bool clamp_receptors = true;
    bool projected_em = true;
    bool saturating_inputs = true;
    bool isolate_branches = false;
    double cross_branch_scale = 0.0;

    bool write_timecourse = true;
    bool write_replicates = false;

    double fucoxanthin = 0.0;
    std::uint64_t seed = 12345;
    std::string scenario = "canonical";
    std::string out_prefix = "results/nfkb";
};

struct StepDiagnostics {
    long long predictor_violations = 0;
    long long predictor_values = 0;
    double max_lower_overshoot = 0.0;
    double max_upper_overshoot = 0.0;
    double total_overshoot = 0.0;

    void observe(double predictor);
    void add(const StepDiagnostics& other);
};

struct ReplicateMetrics {
    int replicate = 0;
    std::uint64_t stream_seed = 0;
    std::vector<double> final_value;
    std::vector<double> auc;
};

struct SimulationResult {
    std::vector<std::string> node_names;
    std::vector<double> time;
    std::vector<std::vector<double>> mean;
    std::vector<std::vector<double>> sd;
    std::vector<std::vector<double>> rep0;
    int save_every = 1;
    int simulated_steps = 0;
    int replicates = 0;
    double dt = 0.0;
    std::vector<double> final_mean;
    std::vector<double> final_sd;
    std::vector<double> final_cv;
    std::vector<double> auc_mean;
    std::vector<double> auc_sd;
    StepDiagnostics diagnostics;
    std::vector<ReplicateMetrics> replicate_metrics;
};

using WeightRealization = std::vector<std::vector<double>>;

class NFkBModel {
public:
    explicit NFkBModel(SimulationConfig cfg);

    SimulationResult run();

    const SimulationConfig& config() const { return cfg_; }
    const std::vector<Node>& nodes() const { return nodes_; }
    const std::vector<Edge>& edges() const { return edges_; }
    std::size_t node_count() const { return nodes_.size(); }

    // Public low-level methods are used by the coupled time-step convergence test.
    WeightRealization draw_weight_realization(std::mt19937_64& rng) const;
    void set_weight_realization(const WeightRealization& weights);
    std::vector<double> initial_state() const;
    std::vector<double> step_with_dW(const std::vector<double>& X,
                                     const std::vector<double>& dW,
                                     StepDiagnostics& diagnostics) const;

private:
    SimulationConfig cfg_;
    std::vector<Node> nodes_;
    std::vector<Edge> edges_;
    std::vector<CompoundTarget> compound_targets_;
    std::map<int, double> receptor_stimulus_;

    void build_nodes();
    void build_edges();
    void build_compound_targets();
    void set_scenario(const std::string& scenario);

    double wk(const std::vector<double>& weights) const;
    double bounded_response(double z) const;
    double compound_fraction(double C, double ic50, double hill) const;
    double edge_scale(const Edge& e) const;

    void reset_fixed_weights();
    void compute_inputs(const std::vector<double>& X,
                        std::vector<double>& A,
                        std::vector<double>& D) const;
};

void write_timecourse_csv(const std::string& filename,
                          const std::vector<std::string>& node_names,
                          const std::vector<double>& time,
                          const std::vector<std::vector<double>>& data);

void write_summary_csv(const std::string& filename,
                       const SimulationResult& result);

void write_replicates_csv(const std::string& filename,
                          const SimulationResult& result);

void write_diagnostics_csv(const std::string& filename,
                           const SimulationConfig& cfg,
                           const SimulationResult& result);

void write_network_csv(const std::string& filename,
                       const std::vector<Node>& nodes,
                       const std::vector<Edge>& edges);

} // namespace nfkb

#endif // REVISION_NFKB_MODEL_HPP
