#ifndef NFKB_MODEL_HPP
#define NFKB_MODEL_HPP

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
#include <vector>
#include <utility>

namespace nfkb {

// -----------------------------------------------------------------------------
// Node indexing for the NF-kB canonical / noncanonical pathway in pathway2.png.
// The labels are intentionally ASCII so that CSV output is easy to parse.
// -----------------------------------------------------------------------------
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

// Role labels are used only for optional branch-isolation tests.
// They do not change the biological graph unless --isolate-branches 1 is used.
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
    double rho = 0.0;     // maximum additional inhibition strength
    double ic50 = 50.0;   // same units as concentration C
    double hill = 1.0;
    std::string label;
};

struct SimulationConfig {
    int steps = 1000;
    int replicates = 1;
    int save_every = 1;          // Save one output row every save_every integration steps.
    int progress_every = 0;      // Print progress every N replicates. 0 disables progress output.
    double dt = 0.01;
    double rb = 0.75;                 // Waller-Kraft restriction parameter
    double fixed_weight = 0.5;
    double sigma = 0.01;
    bool random_weights = false;
    bool feedback_on = false;
    bool clamp_receptors = true;
    bool projected_em = true;
    bool saturating_inputs = true;

    // If true, cross-branch edges are weakened or removed in single-branch tests.
    // This is useful for topology validation. Set to false for biological crosstalk simulations.
    bool isolate_branches = false;
    double cross_branch_scale = 0.0;

    double fucoxanthin = 0.0;
    unsigned int seed = 12345;
    std::string scenario = "canonical";
    std::string out_prefix = "results/nfkb";
};

struct SimulationResult {
    std::vector<std::string> node_names;
    std::vector<double> time;
    std::vector<std::vector<double>> mean; // [saved_time][node]
    std::vector<std::vector<double>> sd;   // [saved_time][node]
    std::vector<std::vector<double>> rep0; // representative trajectory, first replicate
    int save_every = 1;                    // copied from config for provenance
    int simulated_steps = 0;
    double dt = 0.0;
    std::vector<double> final_mean;
    std::vector<double> final_sd;
    std::vector<double> final_cv;
    std::vector<double> auc_mean;
    std::vector<double> auc_sd;
    long long predictor_violations = 0;
    long long predictor_values = 0;
};

class NFkBModel {
public:
    explicit NFkBModel(SimulationConfig cfg);

    SimulationResult run();
    const std::vector<Node>& nodes() const { return nodes_; }
    const std::vector<Edge>& edges() const { return edges_; }

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

    void sample_weights(std::mt19937_64& rng);
    void reset_fixed_weights();

    void compute_inputs(const std::vector<double>& X,
                        std::vector<double>& A,
                        std::vector<double>& D) const;

    std::vector<double> initial_state() const;
    std::vector<double> step(const std::vector<double>& X,
                             std::mt19937_64& rng,
                             long long& violations,
                             long long& nvalues) const;
};

// -----------------------------------------------------------------------------
// CSV writers
// -----------------------------------------------------------------------------
void write_timecourse_csv(const std::string& filename,
                          const std::vector<std::string>& node_names,
                          const std::vector<double>& time,
                          const std::vector<std::vector<double>>& data);

void write_summary_csv(const std::string& filename,
                       const SimulationResult& result);

void write_network_csv(const std::string& filename,
                       const std::vector<Node>& nodes,
                       const std::vector<Edge>& edges);

} // namespace nfkb

#endif // NFKB_MODEL_HPP
