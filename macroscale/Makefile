# Builds the triangle counter. Needs a C++ compiler with OpenMP (g++ on Linux; on macOS use e.g. `make CXX=g++-14` from Homebrew).
CXX ?= g++
triangles: triangles.cpp
	$(CXX) -O3 -fopenmp -o $@ $<

clean:
	rm -f triangles
