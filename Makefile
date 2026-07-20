COMPOSE=docker compose
SERVICE=dev

.PHONY: update-docker build configure compile clean run shell docs test profile clinfo

# Docker

update-docker:
	$(COMPOSE) build

build: update-docker
	$(COMPOSE) run --rm $(SERVICE) cmake -B build
	$(COMPOSE) run --rm $(SERVICE) cmake --build build

configure:
	$(COMPOSE) run --rm $(SERVICE) cmake -B build

compile:
	$(COMPOSE) run --rm $(SERVICE) cmake --build build

clean:
	$(COMPOSE) run --rm $(SERVICE) rm -rf build

shell:
	$(COMPOSE) run --rm $(SERVICE)

# Project

run:
	$(COMPOSE) run --rm $(SERVICE) ./build/simulate

file ?= NONE_PROVIDED
analyze:
	$(COMPOSE) run --rm $(SERVICE) python3 scripts/process_data.py $(file)

test:
	$(COMPOSE) run --rm $(SERVICE) \
	ctest --test-dir build --output-on-failure

docs:
	$(COMPOSE) run --rm $(SERVICE) cmake --build build --target docs

# Profiling

profile:
	$(COMPOSE) run --rm $(SERVICE) \
	valgrind \
	--leak-check=full \
	--track-origins=yes \
	./build/simulation

# OpenCL

clinfo:
	$(COMPOSE) run --rm $(SERVICE) clinfo
