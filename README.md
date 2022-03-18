# Lattice simulation software for the paper "Simulating rare kaon decays using domain wall lattice QCD with physical light quark masses"
[![License: GPL v2](https://img.shields.io/badge/License-GPL_v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

This repository contains the software source code used to produce the data available [here](https://github.com/aportelli/dp008-rarek-data) for the publication [arXiv:2202.08795](https://arxiv.org/abs/2202.08795). This software is based on [Grid](https://github.com/paboyle/Grid) and [Hadrons](https://github.com/aportelli/Hadrons), which are both free software under GPLv2. Both Grid & Hadrons need to be installed as dependencies, and production binaries were using Grid up to commit [a65a49](https://github.com/paboyle/Grid/tree/a65a497baed751eccf6a0e428b30c98e77570416) and Hadrons up to [76b37d](https://github.com/aportelli/Hadrons/tree/76b37db43095854205648462bb3f2d814c69a904).

The code can be compiled using the sequence of commands below
``` bash
./bootstrap.sh
mkdir build
cd build
../configure --with-grid=<Grid prefix> --with-hadrons=<Hadrons prefix> CXX=<compiler used for Grid>
make
```
