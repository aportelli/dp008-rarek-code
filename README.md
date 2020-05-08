# Rare Kaon production code

Install [Grid](https://github.com/paboyle/Grid) 
and [Hadrons](https://github.com/aportelli/Hadrons) on their `develop` branch.
Then:

``` bash
./bootstrap.sh
mkdir build
cd build
../configure --with-hadrons=<Hadrons install prefix> CXX=<compiler used to compile Grid>
make
```