run: build
	./donut
build:
	g++ donut.cpp -o donut
clean:
	rm -f donut donut2