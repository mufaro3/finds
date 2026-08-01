export DOCKER_UID := $(shell id -u)
export DOCKER_GID := $(shell id -g)

COMPOSE=docker compose
SERVICE=dev
DOCKER-SHELL=$(COMPOSE) run --user $(DOCKER_UID):$(DOCKER_GID) --rm $(SERVICE)

.PHONY: help update-docker build configure compile clean run shell docs test profile clinfo paper

help:
	@echo ""
	@echo "FINDS - AVAILABLE COMMANDS"
	@echo ""
	@echo "  [ Building, Development, and Documentation ]"
	@echo ""
	@echo "  make update-docker - Rerun the Dockerfile."
	@echo "  make build         - Build the FINDS library."
	@echo "  make shell         - Open the Container shell."
	@echo "  make clean         - Removes all output directories."
	@echo "  make docs          - Build documentation as HTML and PDF."
	@echo "  make paper         - Build the report PDF."
	@echo "  make clinfo        - Prints CLInfo for OpenCL."
	@echo ""
	@echo "  THE FOLLOWING PROGRAMS ALL REQUIRE (make build):"
	@echo ""
	@echo "  [ No-Argument Programs ]"
	@echo ""
	@echo "  make profile   - Run Valgrind to profile memory usage."
	@echo "  make test      - Run all tests."
	@echo "  make validate  - Produce the validation plots."
	@echo "  make run       - Run the simulation for a given config."
	@echo ""
	@echo "  make analyzize file=[FILE] - Analyze a produced dataset file."
	@echo ""
	@echo "  [ Compiled Programs ]"
	@echo ""
	@echo "  ./build/benchmark            - Benchmarks FINDS computation routines."
	@echo "  ./build/benchmark_derivative - Benchmarks a single derivative."
	@echo "  ./build/simulate             - Equivalent to 'make run'."
	@echo "  ./build/validation           - Equivalent to 'make validate'."

# Docker

update-docker:
	$(COMPOSE) build \
		--build-arg USER_UID=$(DOCKER_UID) \
		--build-arg USER_GID=$(DOCKER_GID)

configure:
	@$(DOCKER-SHELL) cmake -B build

compile:
	@$(DOCKER-SHELL) cmake --build build

build: update-docker configure compile

clean:
	@$(COMPOSE) down
	@$rm -rf build output paper/output

shell:
	@$(DOCKER-SHELL)

jupyter:
	@$(COMPOSE) up jupyter

# Project

validate:
	@$(DOCKER-SHELL) xvfb-run ./build/validation

run:
	@$(DOCKER-SHELL) ./build/simulate

file ?= NONE_PROVIDED
analyze:
	@$(DOCKER-SHELL) python3 scripts/process_data.py $(file)

test:
	@$(DOCKER-SHELL) ctest --test-dir build --output-on-failure -V

docs:
	@$(DOCKER-SHELL) cmake --build build --target docs_pdf
	@$(DOCKER-SHELL) mv build/latex/refman.pdf build/finds.pdf

# paper
paper:
	@$(DOCKER-SHELL) mkdir -p paper/output
	@$(DOCKER-SHELL) sh -c \
	"cd paper && latexmk -f -pdf -output-directory=output paper.tex"
	@$(DOCKER-SHELL) [ -f paper/output/paper.pdf ] && \
	mv paper/output/paper.pdf paper/paper.pdf

# Profiling

VALGRIND := @$(DOCKER-SHELL) valgrind \
		--leak-check=full \
        --show-leak-kinds=all \
		--track-origins=yes

exe ?= ./build/simulate
args ?=
profile:
	@$(VALGRIND) $(exe) $(args)

# OpenCL

clinfo:
	@$(DOCKER-SHELL) clinfo
