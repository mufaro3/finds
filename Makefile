.PHONY: help run dev build down logs shell clean rebuild \
        docs viewdocs test lint validate repl all

FLAKE8_IGNORED := E201,E202,E203,E221,E222,E225,E226,E231,E241,E251,E252,E502,F541,W504,E731
# auto | long | short | line | native | no
tb ?= no
PYTEST_TRACEBACK_MODE := $(tb)

help:
	@echo ""
	@echo "Available commands:"
	@echo "  make run              - Start main"
	@echo "  make dev              - Build + start main"
	@echo "  make build            - Build images"
	@echo "  make down             - Stop containers"
	@echo "  make logs             - Follow logs"
	@echo "  make shell            - Bash shell"
	@echo "  make repl             - Python REPL"
	@echo "  make test [tb=MODE]   - Run pytest"
	@echo ""
	@echo "Traceback options (tb):"
	@echo "  auto   - default pytest behavior"
	@echo "  long   - full traceback"
	@echo "  short  - shorter traceback"
	@echo "  line   - per-line summary"
	@echo "  native - Python native formatting"
	@echo "  no     - no traceback"
	@echo ""
	@echo "Example:"
	@echo "  make test tb=short"
	@echo "  make test tb=no"
	@echo ""
	@echo "  make docs             - Build docs"
	@echo "  make validate         - Run validation"
	@echo "  make lint             - Run linting"
	@echo "  make rebuild          - Rebuild without cache"
	@echo "  make clean            - Remove generated files"
	@echo ""

run:
	@docker compose up main

dev:
	@docker compose run --rm main

build:
	@docker compose build

down:
	@docker compose down

logs:
	@docker compose logs -f

shell:
	@docker compose run --rm main bash

repl:
	@docker compose run --rm repl

docs:
	@docker compose run --rm docs

validate:
	@docker compose run --rm validation

benchmark:
	@docker compose run --rm benchmark

test:
	@docker compose run --rm main \
	python3 -m pytest /finds/test --tb=$(PYTEST_TRACEBACK_MODE)

lint:
	@echo "Running flake8..."
	@docker compose run --rm main \
		python3 -m flake8 --ignore=$(FLAKE8_IGNORED) finds/ test/

	@echo "Checking import order (isort)..."
	@docker compose run --rm main isort finds/

rebuild:
	@docker compose down
	@docker compose build --no-cache

clean:
	@docker compose down -v --remove-orphans
	@rm -rf docs/build/*
	@rm -rf output/*

all: lint build docs test
