CXX ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -pedantic

SRC = src/main_nfkb.cpp src/nfkb_model.cpp
HDR = src/nfkb_model.hpp
BIN = nfkb_sim

# Default pilot settings
DT ?= 0.01
STEPS ?= 1000
MC ?= 1000
SIGMA ?= 0.01
SAVE ?= 1
PROGRESS ?= 0

# Long-run settings. Override from command line, e.g.
# make run_long_mc LONG_MC=500 LONG_STEPS=200000 LONG_SAVE=100
LONG_DT ?= 0.001
LONG_STEPS ?= 100000
LONG_MC ?= 200
LONG_SIGMA ?= 0.01
LONG_SAVE ?= 100
LONG_PROGRESS ?= 20

.PHONY: all clean test run_canonical run_noncanonical run_combined run_det run_iso_det run_mc run_fucox run_all \
        run_long_test run_long_mc run_long_fucox run_train_quick run_train_full figures archive_pilot dirs

all: $(BIN)

dirs:
	mkdir -p results logs figure

$(BIN): $(SRC) $(HDR)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(BIN)

clean:
	rm -f $(BIN)
	rm -f results/*.csv
	rm -rf figure
	rm -f logs/*.log

test: $(BIN) dirs
	./$(BIN) --scenario canonical --sigma 0 --steps 1000 --dt 0.01 --out-prefix results/nfkb

run_canonical: $(BIN) dirs
	./$(BIN) --scenario canonical --sigma 0 --steps $(STEPS) --dt $(DT) --save-every $(SAVE) --out-prefix results/canonical_det

run_noncanonical: $(BIN) dirs
	./$(BIN) --scenario noncanonical --sigma 0 --steps $(STEPS) --dt $(DT) --save-every $(SAVE) --out-prefix results/noncanonical_det

run_combined: $(BIN) dirs
	./$(BIN) --scenario combined --sigma 0 --steps $(STEPS) --dt $(DT) --save-every $(SAVE) --out-prefix results/combined_det

run_det: run_canonical run_noncanonical run_combined

run_iso_det: $(BIN) dirs
	./$(BIN) --scenario canonical --sigma 0 --steps $(STEPS) --dt $(DT) --save-every $(SAVE) --isolate-branches 1 --cross-branch-scale 0 --out-prefix results/canonical_iso_det
	./$(BIN) --scenario noncanonical --sigma 0 --steps $(STEPS) --dt $(DT) --save-every $(SAVE) --isolate-branches 1 --cross-branch-scale 0 --out-prefix results/noncanonical_iso_det
	./$(BIN) --scenario combined --sigma 0 --steps $(STEPS) --dt $(DT) --save-every $(SAVE) --isolate-branches 1 --cross-branch-scale 0 --out-prefix results/combined_iso_det

run_mc: $(BIN) dirs
	./$(BIN) --scenario combined --sigma $(SIGMA) --mc $(MC) --weights random --steps $(STEPS) --dt $(DT) --save-every $(SAVE) --progress-every $(PROGRESS) --out-prefix results/combined_mc

run_fucox: $(BIN) dirs
	./$(BIN) --scenario combined --sigma $(SIGMA) --mc $(MC) --weights random --fucox 25 --steps $(STEPS) --dt $(DT) --save-every $(SAVE) --progress-every $(PROGRESS) --out-prefix results/fucox25
	./$(BIN) --scenario combined --sigma $(SIGMA) --mc $(MC) --weights random --fucox 50 --steps $(STEPS) --dt $(DT) --save-every $(SAVE) --progress-every $(PROGRESS) --out-prefix results/fucox50
	./$(BIN) --scenario combined --sigma $(SIGMA) --mc $(MC) --weights random --fucox 100 --steps $(STEPS) --dt $(DT) --save-every $(SAVE) --progress-every $(PROGRESS) --out-prefix results/fucox100

run_all: $(BIN) run_det run_iso_det run_mc run_fucox

# Long-run quick test: good for confirming that save-every works before starting heavy jobs.
run_long_test: $(BIN) dirs
	./$(BIN) --scenario combined --sigma $(LONG_SIGMA) --mc 50 --weights random --steps 20000 --dt $(LONG_DT) --save-every $(LONG_SAVE) --progress-every 10 --out-prefix results/combined_long_test | tee logs/combined_long_test.log

# Long-run Monte Carlo. Override LONG_MC/LONG_STEPS if the train ride is long.
run_long_mc: $(BIN) dirs
	./$(BIN) --scenario combined --sigma $(LONG_SIGMA) --mc $(LONG_MC) --weights random --steps $(LONG_STEPS) --dt $(LONG_DT) --save-every $(LONG_SAVE) --progress-every $(LONG_PROGRESS) --out-prefix results/combined_long_mc | tee logs/combined_long_mc.log

# Long-run Fucoxanthin dose response including dose 0 control.
run_long_fucox: $(BIN) dirs
	./$(BIN) --scenario combined --sigma $(LONG_SIGMA) --mc $(LONG_MC) --weights random --fucox 0 --steps $(LONG_STEPS) --dt $(LONG_DT) --save-every $(LONG_SAVE) --progress-every $(LONG_PROGRESS) --out-prefix results/fucox0_long | tee logs/fucox0_long.log
	./$(BIN) --scenario combined --sigma $(LONG_SIGMA) --mc $(LONG_MC) --weights random --fucox 25 --steps $(LONG_STEPS) --dt $(LONG_DT) --save-every $(LONG_SAVE) --progress-every $(LONG_PROGRESS) --out-prefix results/fucox25_long | tee logs/fucox25_long.log
	./$(BIN) --scenario combined --sigma $(LONG_SIGMA) --mc $(LONG_MC) --weights random --fucox 50 --steps $(LONG_STEPS) --dt $(LONG_DT) --save-every $(LONG_SAVE) --progress-every $(LONG_PROGRESS) --out-prefix results/fucox50_long | tee logs/fucox50_long.log
	./$(BIN) --scenario combined --sigma $(LONG_SIGMA) --mc $(LONG_MC) --weights random --fucox 100 --steps $(LONG_STEPS) --dt $(LONG_DT) --save-every $(LONG_SAVE) --progress-every $(LONG_PROGRESS) --out-prefix results/fucox100_long | tee logs/fucox100_long.log

# Recommended train workflow: first run a quick long test, then a moderate long dose-response.
run_train_quick: $(BIN) run_long_test

run_train_full: $(BIN) run_long_mc run_long_fucox figures

figures: dirs
	python scripts/auto_plot_nfkb.py

archive_pilot:
	mkdir -p archive_pilot
	cp -R results figure archive_pilot/ 2>/dev/null || true
