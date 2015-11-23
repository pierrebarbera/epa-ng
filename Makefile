
all: build/CMakeCache.txt run_make
.PHONY: all

# Run cmake if not yet done or if CMakeLists.txt has changed.
build/CMakeCache.txt: CMakeLists.txt
	@echo "Running cmake"
	@mkdir -p build
	@cd build && cmake ..

run_make: build/CMakeCache.txt
	@echo "Running make"
	@make -C build
.PHONY: run_make

new_files:
	@touch src/CMakeLists.txt
	@touch test/src/CMakeLists.txt
	@make -C build
.PHONY: new_files

test: new_files
	@./test/bin/epa_test
.PHONY: test

clean:
	@echo "Cleaning"
	@rm -rf build
	@rm -rf bin
	@rm -rf test/bin
.PHONY: clean
