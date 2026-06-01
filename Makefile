build:
	docker build -t fish .
run: build
	docker run fish
