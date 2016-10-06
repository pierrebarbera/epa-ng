PLLMOD=libs/pll-modules
PLL=${PLLMOD}/libs/libpll

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

update:
	@touch src/CMakeLists.txt
	@touch test/src/CMakeLists.txt
	@make -C build
.PHONY: update

test: update
	@./test/bin/epa_test
.PHONY: test

pll:
	@echo "Don't forget to set the fastblo flag!"
	mkdir -p bin
	cd ${PLLMOD} && ./install-with-libpll.sh .. 
.PHONY: pll
		
clean:
	@echo "Cleaning"
	@rm -rf build
	@rm -rf bin
	@rm -rf test/bin
.PHONY: clean
