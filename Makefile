.PHONY: all, run_make, test, clean

test: all
	./bin/epa_test

all: build/CMakeCache.txt run_make

# Run cmake if not yet done or if CMakeLists.txt has changed.
build/CMakeCache.txt: CMakeLists.txt
	@echo "Running cmake"
	@mkdir -p build
	@cd build && cmake ..

run_make: build/CMakeCache.txt
	@echo "Running make"
	@make -C build

clean:
	@echo "Cleaning"
	@rm -rf build
	@rm -rf bin
