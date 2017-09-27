PLLMOD=libs/pll-modules
PLL=${PLLMOD}/libs/libpll
GENESIS=libs/genesis

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

genesis:
	mkdir -p bin
	cd ${GENESIS} && $(MAKE)
.PHONY: genesis

genesis_update:
	mkdir -p bin
	cd ${GENESIS} && $(MAKE) update
.PHONY: genesis

genesis_clean:
	cd ${GENESIS} && $(MAKE) clean
.PHONY: genesis

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
QRY=$(TEST)/1k_query_100.fasta.bin
BINFILE=$(TEST)/epa_binary_file
OUTDIR=/tmp/epa

BINARY_WRITE= -t $(TREE) -s $(REF) -B -w $(OUTDIR) --verbose $(F)
BINARY_READ=-b $(BINFILE) -q $(QRY) -w $(OUTDIR) -g 0.99 --verbose $(F)
NORM_TEST=-t $(TREE) -s $(REF) -q $(QRY) -w $(OUTDIR) -g 0.99 --chunk-size=10 --verbose $(F)

test: #update
	mkdir -p $(OUTDIR)
	rm -f $(OUTDIR)/*
	$(EPABIN) $(NORM_TEST) #--threads 4
.PHONY: test

bintest: update
	mkdir -p $(OUTDIR)
	rm -f $(OUTDIR)/*
	$(EPABIN) $(BINARY_WRITE)
	$(EPABIN) $(BINARY_READ)
.PHONY: bintest

mpi_test: update
	mkdir -p $(OUTDIR)
	rm -f $(OUTDIR)/*
	mpirun -n 4 $(EPABIN) $(BINARY_READ)
.PHONY: mpi_test
