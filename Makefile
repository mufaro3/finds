.PHONY: help run build logs shell clean docs test lint validate benchmark \
        stats repl jupyter paper derivative
help:
	@echo ""
	@echo "FINDS - AVAILABLE COMMANDS"
	@echo ""
	@echo "  [Building, Running, and Documentation]"
	@echo ""
	@echo "  make run         - Build and run main."
	@echo "  make build       - Build the container image."
	@echo "  make docs        - Build documentation as HTML and PDF."
	@echo "  make paper       - Build the report PDF."
	@echo ""
	@echo "  [Benchmarking, Testing and Validation]"
	@echo ""
	@echo "  make benchmark [numba=]  - Run benchmarking module."
	@echo ""
	@echo "    Use Numba (numba): whether or not to use Numba performance"
	@echo "    optimizations (disabled by default)."
	@echo ""
	@echo "      Example:"
	@echo "        make benchmark numba=enable"
	@echo "        make benchmark numba=disable"
	@echo ""
	@echo "  make derivative [theta=] [n=] - Time a singule derivative"
	@echo "                                  calculation"
	@echo "    Barnes-Hut Minimum Ratio (theta) - Where theta = 0, then"
	@echo "    the Barnes-Hut approximation is not used (0 by default)."
	@echo ""
	@echo "    Number of fish (n) - (10 by default)"
	@echo ""
	@echo "      Example:"
	@echo "        make derivative theta=0.75 n=100000"
	@echo "        make derivative theta=0.5  n=1000"
	@echo ""
	@echo "  make validate         - Run validation module."
	@echo "  make test [tb=] [tf=] - Run tests."
	@echo ""
	@echo "  Traceback options (tb):"
	@echo "    auto   - Default pytest behavior."
	@echo "    long   - Full traceback."
	@echo "    short  - Shorter traceback."
	@echo "    line   - Per-line summary."
	@echo "    native - Python native formatting."
	@echo "    no     - No traceback (default)."
	@echo ""
	@echo "    Example:"
	@echo "      make test tb=short"
	@echo "      make test tb=no"
	@echo ""
	@echo "  Test file input (tf): All tests by default, or specify a file"
	@echo "  under the finds/test/ directory."
	@echo ""
	@echo "    Example:"
	@echo "      make test tf=test_calculations.py"
	@echo ""
	@echo "  [Debugging]"
	@echo ""
	@echo "  make logs             - Open Docker logs."
	@echo "  make shell            - Bash shell."
	@echo "  make repl             - Python REPL."
	@echo "  make jupyter          - Launch Jupyter Lab."
	@echo "  make stats            - Display container performance."
	@echo "  make lint             - Run static code checkers."
	@echo "  make clean            - Remove generated files."
	@echo ""

DRUN := docker compose run --rm
MAIN := $(DRUN) main
PY := $(MAIN) python3 -m

run:
	@$(MAIN)

build:
	@docker compose build

logs:
	@docker compose logs -f

shell:
	@$(MAIN) bash

repl:
	@$(DRUN) repl

docs:
	@$(DRUN) docs

paper:
	@mkdir -p paper/output
	@cd paper && latexmk -pdf -output-directory=output paper.tex
	@[ -f paper/output/paper.pdf ] && mv paper/output/paper.pdf paper/paper.pdf

validate:
	@$(DRUN) validation

# enable | disable
numba ?= disable
BENCHMARK_USENUMBA := $(numba)
benchmark:
	@$(DRUN) benchmark python3 scripts/benchmark.py $(BENCHMARK_USENUMBA)

theta ?= 0
n ?= 10
DERIVATIVE_THETA := $(theta)
DERIVATIVE_N := $(n)
derivative:
	@$(DRUN) benchmark python3 scripts/time_single_derivative.py \
		--theta $(DERIVATIVE_THETA) $(DERIVATIVE_N)

jupyter:
	@docker compose up jupyter

DSTATS_FORMAT := "table {{.Container}}\t{{.CPUPerc}}\t{{.MemUsage}}"
stats:
	@docker stats --format $(DSTATS_FORMAT)

# auto | long | short | line | native | no
tb ?= no
PYTEST_TRACEBACK_MODE := $(tb)
COVERAGE_PACKAGE := finds
COVERAGE_DIR := /finds/output/coverage

# example: test_simulation.py
tf ?= .
TESTING_DIR_OR_FILE := $(tf)

test:
	@$(PY) pytest /finds/test/$(TESTING_DIR_OR_FILE) \
		--tb=$(PYTEST_TRACEBACK_MODE) \
		--cov=$(COVERAGE_PACKAGE) \
		--cov-report=html:$(COVERAGE_DIR)

FLAKE8_IGNORED := \
E201,E202,E203,E221,E222,E225,E226,E231,E241,E251,E252,E402,E502,F541,W504,E731

lint:
	@$(PY) flake8 --ignore=$(FLAKE8_IGNORED) finds/ test/ scripts/
	@$(PY) isort finds/ test/ scripts/
	@$(PY) mypy finds/ test/ scripts/ --ignore-missing-imports

clean:
	@docker compose down -v --remove-orphans
	@rm -rf docs/build/*
	@rm -rf output/*
