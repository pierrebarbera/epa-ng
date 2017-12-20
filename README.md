# EPA-ng

1. **[Introduction](#introduction)**
2. **[Build Instructions](#build-instructions)**
3. **[Installation](#installation)**
4. **[Usage](#usage)**

## IMPORTANT
The by far easiest way to clone this repository is to use the following command
```
git clone --recursive https://github.com/Pbdas/epa.git
```

If that is not an option (perhaps you downloaded the zip file), you can fix the submodules by running

```
git submodule update --init --recursive
```

### SUPPORT

The most reliable way to get in touch with us is to head over to the [RAxML Google Group](https://groups.google.com/forum/#!forum/raxml). You can also search the history of the list for your particular question.

Alternatively I've created a [gitter chat room](https://gitter.im/epa-ng/Lobby) where I can usually be found during office hours.

### DISCLAIMER

This tool is still in an active, *beta* phase of development. Suggestions, bug reports and constructive comments are more than encuraged! Please do so in the [issues section](https://github.com/Pbdas/epa/issues).

## Introduction

`EPA-ng` is a complete rewrite of the [Evolutionary Placement Algorithm (EPA)](https://academic.oup.com/sysbio/article/60/3/291/1667010/Performance-Accuracy-and-Web-Server-for), previously implemented in [RAxML](https://github.com/stamatak/standard-RAxML). It uses [libpll](https://github.com/xflouris/libpll) to perform maximum likelihood-based phylogenetic placement of genetic sequences on a user-supplied reference tree and alignment.

### What can EPA-ng do?

- do phylogenetic placement using the **GTR+GAMMA ML model only** (for now)
- take as input **separated reference and query alignment files**, in the **fasta** format (for now)
- handle **DNA** data (protein data coming VERY soon)
- distribute the work to the **cluster**, with the choice of two different modes:
  - normal: assumes the input tree and alignment's memory footprint is small enough to fit into the memory of each compute node
  - pipeline: for large input trees and alignments: distributes the memory footprint across the compute nodes, at the expense of parallel efficiency (about 20% lower)
- **prepare inputs** for the cluster:
  - precompute reference tree and alignment and save it to a *binary file* that allows random access by the different compute nodes
  - convert query fasta file into a random access, binary encoded file called a `bfast`-file
  - these two precomputations are simple, and reccomended in general! Use them!
- output the placement results in the [jplace format](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0031009) ready for downstream analysis by frameworks such as [genesis](https://github.com/lczech/genesis)

## Build Instructions

After cloning with the above command, ensure the following packages are installed or otherwise available (relevant modules loaded on your cluster):

```
sudo apt-get install autotools-dev libtool flex bison cmake
```

Once these dependencies are available, you need to ensure that your compiler is recent enough, as EPA-ng is built using C++14 features. The minimum required versions are as follows:

| Compiler | Min. Version |
| - | - |
| gcc | 4.9.2  |
| clang  | 3.8  |
| icc | 16 |

Any one of these compilers will be sufficient. gcc is the most wide spread, and current versions of Ubuntu have gcc versions exceeding the minimum.

Now it's time to build the program. First, we build `libpll`:

```
cd epa
make pll
```

This will configure, build, and install the various `libpll` libraries into a subfolder of the `epa-ng` directory.
**Make sure this step has not produced any errors!** Otherwise, we cannot procede.

Next, we will build `EPA-ng` itself:

```
make
```

Thats it! The executable should now be located in the `epa-ng/bin/` folder.

### Apple

Supported in theory, but currently work in progress. In principle same procedure as under Linux.

### Windows

No support yet, will gaugue interest first.

## Installation

**NOTE:** this is only necessary if you want to be able to just type `epa-ng -s ...` from anywhere without having to specify the full path to the executable every time (which is a totally valid way of doing it also).

There are multiple ways of ensuring your system knows where to look for the `EPA-ng` executable.
Here I will assume that the user has very restricted rights, as is usually the case when using a cluster/supercomputer.

### Linux

Often, the system will have the following location added to the search path by default:

```
/home/<username>/bin
```

We can exploit this by navigating there (possibly creating the folder first), and creating a link to where the executable is located:
(here we assume epa was cloned into the base home directory)

```shell
mkdir -p ~/bin
cd ~/bin
ln -s ../epa-ng/bin/epa-ng
```

## Usage

`EPA-ng` is used from the command line, as the main use-case is processing large amounts of data using a supercomputing cluster.

Here is a list of the most basic arguments you will use:

| Flag | Long Flag | Meaning |
| - | - | - |
| -s | --ref-msa | reference MSA (fasta)  |
| -t | --tree | reference Tree (newick)  |
| -q | --query | query sequences (fasta or [bfast](#converting-the-query-file)) |
| -w | --outdir | output directory (default: current directory) |
| -T | --threads | number of threads to use |

For a full overview of command line options either run `EPA-ng` with no input, or with the flag `-h` (or `--help`).

### Basic

On a single computer, an example execution might look like this:

```
epa-ng -s reference_msa.fasta -t reference_tree.newick -q query.fasta -w /some/output/dir/
```

Note that this will use as many threads as specified by the environment variable `OMP_NUM_THREADS`. 
Usually this defaults to the number of cores.
Note however, that no speedup is to be expected from hyperthreads, meaning the number of threads should be set to the number of physical cores.

**NOTE** also, that currently, if you supply a query file in the fasta format, `EPA-ng` will automatically convert the file to the bfast format and write this file to the output directory.

### Advanced

Overview of advanced features:

| Flag | Long Flag | Meaning |
| - | - | - |
| -O | --opt-ref-tree | [optimize the reference tree](#reference-tree-optimization)  |
| -g | --dyn-heur | use dynamic [heuristic](#configuring-the-heuristic-preplacement) (default)  |
| -G | --fix-heur | use fixed [heuristic](#configuring-the-heuristic-preplacement) |
|  | --no-heur | use no [heuristic](#configuring-the-heuristic-preplacement) |
| -B | --dump-binary | [produce binary tree file](#precomputing-the-binary-reference-tree) |
| -c | --bfast | [convert query fasta to binary format](#converting-the-query-file) |
|  | --pipeline | [use pipeline mode](#pipeline-parallelism) |

[The description of basic cluster usage starts here](#running-on-the-cluster)


#### Reference Tree Optimization

When supplying the `-O` (or `--opt-ref-tree`) flag to `EPA-ng`, an initial round of model parameter and tree branch length optimization is performed.

#### Configuring the Heuristic Preplacement

By default, `EPA-ng` performs placement of a sequence in two stages: first selecting promising branches quickly (preplacement), then evaluating the selected branches in greater detail.

`EPA-ng` currently offers two ways of selecting these candidates, based on their Likelihood Weight Ratios (LWR), or ranking of best to worst.

The default is the *accumulated threshold* method, in which branches are added to the set of candidates until the sum of their LWR exceed a user specified threshold.
The flag controlling this mode is `-g` (or `--dyn-heur`), with a default setting of `0.99`, corresponding to a covered likelihood weight of 99%.

There is also a second mode, which functions identically to the candidate selection mode in the original implementation of the `EPA` in `RAxML`.
Here again the branches are sorted by the LWR of the placement of a sequence.
Then, the top x% of the total number of branches are selected into the set of candidates.
Like in `RAxML`, this behavior is controlled via the `-G` (or `--fix-heur`) flag.

Lastly, to disable the heuristic completely, you can simply supply the `--no-heur` flag.
Be warned however: doing so will be significantly more computationally demanding.
Our advice is to use the heuristic, as it sacrifices only insignificant amounts of accuracy for greatly improved speed.


### Cluster usage

Before using the cluster version of `EPA-ng`, the input files must be preprocessed.
In general, these preprocessed files also enable more streamlined re-runs of experiments and are reccomended for use, always.

#### Precomputing the binary reference tree

Internally, when calling `EPA-ng` with `-s reference_msa.fasta -t reference_tree.newick`, both files are parsed and an internal tree structure is created.
Additionally, if the `-O` (or `--opt-ref-tree`) flag was supplied, a likelihood model parameter, and branch length, optimization is performed.
Then, every possible partial likelihood vector (CLV) on the tree is computed and stored, such that all CLV needed to perform phylogenetic placement exist.

To avoid having to repeat these steps between multiple runs using the same reference data, this precomputed tree can be emitted using the `-B` (or `--dump-binary`) flag.
When this flag is supplied, the program will produce a file called `epa_binary_file` in the output directory, and then terminate.

In subequent runs, instead of `-s ... -t ...` you may then simply supply the binary tree file: `-b epa_binary_file`

#### Converting the query file

You may also explicitly convert the input query fasta file to our internal fasta format.
This format is binary encoded (reducing the size by half) and randomly accessible.
Again, using this format is highly reccomended, and required to use the MPI version.

To convert the fasta file, simply run the program with the query file specified thusly:

```
epa-ng -c query.fasta -w /some/output/dir
```

This will produce a file called `query.fasta.bin` in the specified output directory.

#### Running on the cluster

To use distributed parallelism in `EPA-ng`, first we must re-compile the program with MPI enabled.
This requires a version of MPI to be loaded/installed on your system.
The only additional requirement `EPA-ng` has, is that the compiler that is loaded in conjunction with MPI satisfies the [minimum version requirements](#build-instructions).
Often this can be assured by the order in which the relevant modules are loaded on the cluster: first MPI, then the compiler.
However we reccomend you contact your support team should this cause issues for you.

The actual compilation is very straight-forward:
```
make clean && make EPA_HYBRID=1
```

This will attempt to compile the program with both MPI and OpenMP, as the most efficient way to run the program is to map one MPI rank per node (good alternative: one rank per socket!), each rank starting as many threads as there are *physical* cores.

In your job submission script, you can then call the program in a highly similar way to before:

```
mpirun epa-ng -b epa_binary_file -q query.fasta.bin -w ./some/output/dir
```

#### Pipeline parallelism

There is one additional mode for the cluster, which is reccomended for cases when the memory footprint (equivalent to the [`epa_binary_file`](#precomputing-the-binary-reference-tree) size) exceeds the limitations of a single compute node.
You can enable it using the `--pipeline` flag.
Contrary to the *normal* distributed parallel mode, the pipeline mode attempts to distribute both the tree and the query sequences across many nodes.
Doing so enables `EPA-ng` to place on the very largest of trees and alignment sizes, at the expense of parallel efficiency.
Due to the much higher complexity of the code involved this feature can be considered in an **alpha** state.

Development on this feature may or may not continue, as we have some other ideas of handling large trees in the works. 
