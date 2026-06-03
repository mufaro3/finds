.PHONY: help run build dev down logs shell clean rebuild

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
	@echo ""

run:
	@docker compose up

dev:
	@docker compose up --build

build:
	@docker compose build

down:
	@docker compose down

logs:
	@docker compose logs -f

shell:
	@docker compose exec main bash

rebuild:
	@docker compose down
	@docker compose build --no-cache

clean:
	@docker compose down -v --remove-orphans
