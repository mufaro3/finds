.PHONY: help run dev build down logs shell clean rebuild \
        docs viewdocs test lint validate repl all

ENV_FILE := .env
FLAKE8_IGNORED := E201,E202,E203,E221,E222,E225,E226,E231,E241,E251,E252,\
	          E502,F541,W504,E731
SERVICE ?= main

help:
	@echo ""
	@echo "Available commands:"
	@echo "  make run        - Start main"
	@echo "  make dev        - Build + start main"
	@echo "  make build      - Build images"
	@echo "  make down       - Stop containers"
	@echo "  make logs       - Follow logs"
	@echo "  make shell      - Bash shell"
	@echo "  make repl       - Python REPL"
	@echo "  make test       - Run pytest"
	@echo "  make docs       - Build docs"
	@echo "  make validate   - Run validation"
	@echo "  make lint       - Run linting"
	@echo "  make rebuild    - Rebuild without cache"
	@echo "  make clean      - Remove generated files"
	@echo ""

run:
	@docker compose up main

dev:
	@docker compose up --build main

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

test:
	@docker compose run --rm test

validate:
	@docker compose run --rm validation

docs:
	@docker compose run --rm docs

viewdocs: docs
	@zathura docs/build/finds.pdf &

lint:
	@echo "Running flake8..."
	@docker compose run --rm main \
		flake8 --ignore=$(FLAKE8_IGNORED) finds/ test/

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
