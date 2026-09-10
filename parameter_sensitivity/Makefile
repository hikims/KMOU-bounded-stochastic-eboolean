CXX ?= c++
PYTHON ?= python3
CPPFLAGS += -Isrc
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -pedantic
LDFLAGS ?=

NFKB_EXE := nfkb_revision_sim
LOW_EXE := low_activation_validation
DT_EXE := dt_convergence

.PHONY: all check regression-check smoke clean \
        sensitivity-smoke sensitivity-smoke-figures sensitivity-pilot sensitivity-production \
        sensitivity-extend sensitivity-figures \
        low-activation-smoke low-activation-smoke-figures low-activation-production \
        dt-smoke dt-smoke-figures dt-pilot dt-production-control dt-production-high

all: $(NFKB_EXE) $(LOW_EXE) $(DT_EXE)

$(NFKB_EXE): src/main_nfkb.cpp src/nfkb_model.cpp src/nfkb_model.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) src/main_nfkb.cpp src/nfkb_model.cpp -o $@ $(LDFLAGS)

$(LOW_EXE): src/low_activation_validation.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) src/low_activation_validation.cpp -o $@ $(LDFLAGS)

$(DT_EXE): src/dt_convergence.cpp src/nfkb_model.cpp src/nfkb_model.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) src/dt_convergence.cpp src/nfkb_model.cpp -o $@ $(LDFLAGS)

check: all
	$(PYTHON) -m py_compile scripts/*.py tests/*.py
	./$(NFKB_EXE) --help >/dev/null
	./$(LOW_EXE) --help >/dev/null
	./$(DT_EXE) --help >/dev/null

regression-check: $(NFKB_EXE)
	$(MAKE) -C ../fucoxanthin_dose_response nfkb_sim
	$(PYTHON) tests/check_original_equivalence.py

# A short end-to-end test. It is intended to validate the workflow, not to
# generate manuscript-quality Monte Carlo estimates.
smoke: sensitivity-smoke-figures low-activation-smoke-figures dt-smoke-figures

sensitivity-smoke: $(NFKB_EXE)
	rm -rf results/parameter_sensitivity/smoke
	$(PYTHON) scripts/run_parameter_sensitivity.py \
		--profile smoke --run-id smoke --exe ./$(NFKB_EXE)

sensitivity-smoke-figures: sensitivity-smoke
	$(PYTHON) scripts/plot_parameter_sensitivity.py \
		--run-dir results/parameter_sensitivity/smoke \
		--out-dir figures/sensitivity_smoke

sensitivity-pilot: $(NFKB_EXE)
	$(PYTHON) scripts/run_parameter_sensitivity.py \
		--profile pilot --exe ./$(NFKB_EXE)

# RUN_ID may be supplied, for example:
#   make sensitivity-production RUN_ID=production_v1
RUN_ID ?= production_$(shell date +%Y%m%d_%H%M%S)
sensitivity-production: $(NFKB_EXE)
	$(PYTHON) scripts/run_parameter_sensitivity.py \
		--profile production --run-id $(RUN_ID) --exe ./$(NFKB_EXE)

# Extend a completed production run after the grid changed (e.g. the widened
# noise scan in config/sensitivity_grid.csv). Chunk outputs that already exist
# in SOURCE_RUN_ID are copied into RUN_ID and skipped, so only the new
# conditions are simulated; all conditions are then merged and bootstrapped.
#   make sensitivity-extend RUN_ID=production_v2 SOURCE_RUN_ID=production_v1
SOURCE_RUN_ID ?= production_v1
sensitivity-extend: $(NFKB_EXE)
	$(PYTHON) scripts/run_parameter_sensitivity.py \
		--profile production --run-id $(RUN_ID) --exe ./$(NFKB_EXE) \
		--reuse-chunks-from $(SOURCE_RUN_ID)

sensitivity-figures:
	$(PYTHON) scripts/plot_parameter_sensitivity.py \
		--run-dir results/parameter_sensitivity/$(RUN_ID) \
		--out-dir figures/main_text

low-activation-smoke: $(LOW_EXE)
	rm -rf results/low_activation_smoke
	./$(LOW_EXE) --inputs 0.01,0.1,1 --dt 0.01 --time 1 \
		--mc 50 --save-every 10 --out-dir results/low_activation_smoke

low-activation-smoke-figures: low-activation-smoke
	$(PYTHON) scripts/plot_low_activation_validation.py \
		--data-dir results/low_activation_smoke \
		--out-dir figures/low_activation_smoke

low-activation-production: $(LOW_EXE)
	./$(LOW_EXE) --inputs 0.01,0.03,0.1,0.3,1,3 \
		--dt 0.001 --time 20 --mc 5000 --save-every 100 \
		--out-dir results/low_activation
	$(PYTHON) scripts/plot_low_activation_validation.py \
		--data-dir results/low_activation \
		--out-dir figures/low_activation

dt-smoke: $(DT_EXE)
	rm -rf results/dt_convergence_smoke figures/dt_convergence_smoke
	./$(DT_EXE) --dts 0.02,0.01,0.005 --time 0.2 --mc 4 \
		--progress-every 0 --fucox 0 \
		--out-prefix results/dt_convergence_smoke/control
	./$(DT_EXE) --dts 0.02,0.01,0.005 --time 0.2 --mc 4 \
		--progress-every 0 --fucox 100 \
		--out-prefix results/dt_convergence_smoke/high_dose

dt-smoke-figures: dt-smoke
	$(PYTHON) scripts/plot_dt_convergence.py \
		--prefix results/dt_convergence_smoke/control \
		--prefix results/dt_convergence_smoke/high_dose \
		--out-dir figures/dt_convergence_smoke

dt-pilot: $(DT_EXE)
	./$(DT_EXE) --dts 0.002,0.001,0.0005 --time 20 --mc 100 \
		--fucox 0 --out-prefix results/dt_convergence_pilot/control
	./$(DT_EXE) --dts 0.002,0.001,0.0005 --time 20 --mc 100 \
		--fucox 100 --out-prefix results/dt_convergence_pilot/high_dose
	$(PYTHON) scripts/plot_dt_convergence.py \
		--prefix results/dt_convergence_pilot/control \
		--prefix results/dt_convergence_pilot/high_dose \
		--out-dir figures/dt_convergence_pilot

# These full-horizon checks can be expensive. Run them individually or on a
# compute server, then call plot_dt_convergence.py with both output prefixes.
dt-production-control: $(DT_EXE)
	./$(DT_EXE) --dts 0.002,0.001,0.0005 --time 500 --mc 50 \
		--fucox 0 --out-prefix results/dt_convergence/control

dt-production-high: $(DT_EXE)
	./$(DT_EXE) --dts 0.002,0.001,0.0005 --time 500 --mc 50 \
		--fucox 100 --out-prefix results/dt_convergence/high_dose

clean:
	rm -f $(NFKB_EXE) $(LOW_EXE) $(DT_EXE) *.o
	rm -rf __pycache__ scripts/__pycache__ tests/__pycache__
