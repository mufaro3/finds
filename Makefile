.PHONY: help run build logs shell clean \
        docs test lint validate benchmark stats repl
help:
	@echo ""
	@echo "(FINDS) AVAILABLE COMMANDS"
	@echo ""
	@echo "  [Building, Running, and Documentation]"
	@echo ""
	@echo "  make run              - Build and run main."
	@echo "  make build            - Build the container image."
	@echo "  make docs             - Build documentation as HTML and PDF."
	@echo ""
	@echo "  [Testing and Validation]"
	@echo ""
	@echo "  make benchmark        - Run benchmarking module."
	@echo "  make validate         - Run validation module."
	@echo "  make test [tb=MODE]   - Run all tests."
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
	@echo "  [Debugging]"
	@echo ""
	@echo "  make logs             - Open Docker logs."
	@echo "  make shell            - Bash shell."
	@echo "  make repl             - Python REPL."
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

validate:
	@$(DRUN) validation

benchmark:
	@$(DRUN) benchmark

DSTATS_FORMAT := "table {{.Container}}\t{{.CPUPerc}}\t{{.MemUsage}}"
stats:
	@docker stats --format $(DSTATS_FORMAT)

# auto | long | short | line | native | no
tb ?= no
PYTEST_TRACEBACK_MODE := $(tb)
test:
	@$(PY) pytest /finds/test --tb=$(PYTEST_TRACEBACK_MODE)

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
