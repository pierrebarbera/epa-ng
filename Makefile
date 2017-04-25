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
	$(MAKE) -C build 
.PHONY: run_make

update:
	@touch src/CMakeLists.txt
	@touch test/src/CMakeLists.txt
	$(MAKE) -C build
.PHONY: update

unittest: update
	@./test/bin/epa_test
.PHONY: test

pll:
	mkdir -p bin
	cd ${PLLMOD} && ./install-with-libpll.sh .. 
.PHONY: pll

pll_clean:
	cd ${PLL} && make clean
	cd ${PLLMOD} && make clean
.PHONY: pll_clean
		
clean:
	@echo "Cleaning"
	@rm -rf build
	@rm -rf bin
	@rm -rf test/bin
.PHONY: clean

#======================================
#===		Test commands follow				===
#======================================
EPABIN=./bin/epa
TEST=test/data/lucas
TREE=$(TEST)/20k.newick
REF=$(TEST)/1k_reference.fasta
QRY=$(TEST)/1k_query.fasta
OUTDIR=/tmp/epa

BINARY_WRITE= -t $(TREE) -s $(REF) -B -w $(OUTDIR) $(F)
BINARY_READ=-b $(OUTDIR)/epa_binary_file -q $(QRY) -w $(OUTDIR) -g 0.99 --filter-min-lwr 0.0 $(F)
NORM_TEST=-t $(TREE) -s $(REF) -q $(QRY) -w $(OUTDIR) -g 0.99 --chunk-size=100 $(F)

test: update
	mkdir -p $(OUTDIR)
	-rm -f $(OUTDIR)/*
	$(EPABIN) $(NORM_TEST)
.PHONY: test

mpi_test: update
	mkdir -p $(OUTDIR)
	-rm -f $(OUTDIR)/*
	mpirun -n 5 $(EPABIN) $(NORM_TEST)
.PHONY: mpi_test