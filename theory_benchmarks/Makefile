CXX ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -pedantic
PYTHON ?= python

BIN = theory_benchmarks
SRC = src/theory_benchmarks.cpp

THEORY_THREADS ?= 24
BOUNDARY_MC ?= 10000
MOMENT_MC ?= 50000

.PHONY: all build run figures clean test

all: figures

build: $(BIN)

$(BIN): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(BIN)

test: $(BIN)
	mkdir -p results logs figure
	./$(BIN) \
		--out-prefix results/theory \
		--threads 8 \
		--boundary-mc 1000 \
		--moment-mc 3000 \
		2>&1 | tee logs/theory_benchmarks_test.log
	python scripts/plot_theory_benchmarks.py

run: $(BIN)
	mkdir -p results logs figure
	./$(BIN) \
		--out-prefix results/theory \
		--threads $(THEORY_THREADS) \
		--boundary-mc $(BOUNDARY_MC) \
		--moment-mc $(MOMENT_MC) \
		2>&1 | tee logs/theory_benchmarks.log

figures: run
	python scripts/plot_theory_benchmarks.py

clean:
	rm -f $(BIN)
	rm -f results/theory_*.csv
	rm -rf figure/theory_benchmarks
	rm -f logs/theory_benchmarks.log
	rm -f logs/theory_benchmarks_test.log