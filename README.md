# EPA-ng

## IMPORTANT
The by far easiest way to clone this repository is to use the following command

`git clone --recursive https://github.com/Pbdas/epa.git`

(unless you are brave enough to deal with `git submodules` directly)

## Introduction

EPA-ng is a complete rewrite of the [Evolutionary Placement Algorithm (EPA)](https://academic.oup.com/sysbio/article/60/3/291/1667010/Performance-Accuracy-and-Web-Server-for), previously implemented in [RAxML](https://github.com/stamatak/standard-RAxML). It uses [libpll](https://github.com/xflouris/libpll) to perform maximum likelihood-based phylogenetic placement of genetic sequences on a user-supplied reference tree and alignment.

## Build Instructions

### Linux

After cloning with the above command, ensure the following packages are installed or otherwise available (relevant modules loaded on your cluster):

`sudo apt-get install autotools-dev libtool flex bison cmake`

Once these dependencies are available, you need to ensure that your compiler is recent enough, as EPA-ng is built using C++14 features. The minimum required versions are as follows:

| Compiler | Min. Version |
| -------- | ------------ |
| gcc | 4.9.2  |
| clang  | 3.8  |
| icc | 16 |

Any one of these compilers will be sufficient. gcc is the most wide spread, and current versions of Ubuntu have gcc versions exceeding the minimum.

Now it's time to build the program. First, we build `libpll`:

```
cd epa
make pll
```

This will configure, build, and install the various `libpll` libraries into a subfolder of the `epa` directory.
**Make sure this step has not produced any errors!** Otherwise, we cannot procede.

Next, we will build `EPA-ng` itself:

```
make
```

Thats it! The executable should now be located in the `epa/bin/` folder.

## Usage

### Basic: Single Core


