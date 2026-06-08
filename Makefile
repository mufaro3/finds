.PHONY: help run build dev down logs shell clean rebuild docs test all lint

help:
	@echo ""
	@echo "Available commands:"
	@echo "  make run      - Start containers (no build)"
	@echo "  make dev      - Build + start containers"
	@echo "  make build    - Build images only"
	@echo "  make down     - Stop containers"
	@echo "  make logs     - Follow container logs"
	@echo "  make shell    - Open shell in main container"
	@echo "  make rebuild  - Full rebuild (no cache)"
	@echo "  make clean    - Remove containers + volumes"
	@echo "  make docs     - Build the documentation"
	@echo "  make test     - Rebuild and run pytest test suite"
	@echo "  make all      - Rebuild, run all tests, and build documentation"
	@echo ""

run:
	@docker compose up
	@docker compose run --rm main chown -R $$(id -u):$$(id -g) docs/build

dev:
	@docker compose up --build
	@docker compose run --rm main chown -R $$(id -u):$$(id -g) docs/build

build:
	@docker compose build

down:
	@docker compose down

logs:
	@docker compose logs -f

shell:
	@docker compose run --rm main bash

rebuild:
	@docker compose down
	@docker compose build --no-cache

clean:
	@docker compose down -v --remove-orphans
	@rm -rf docs/build/latex

docs:
	@docker compose run --rm main make -C docs latexpdf
	@docker compose run --rm main chown -R $$(id -u):$$(id -g) docs/build

viewdocs: docs
	@zathura docs/build/finds.pdf &

test:
	@docker compose run --rm main pytest --tb=long

lint:
	@echo "Running flake8..."
	@docker compose run --rm main flake8 --ignore=E201,E202,E203,E221,E225,E226,E231,E241,E252,E502,F541,W504,E731 finds/

	@echo "Checking import order (isort)..."
	@docker compose run --rm main isort finds/

all: lint build docs test
