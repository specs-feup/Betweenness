# Betweenness centrality Code Repository
This repository contains C++ implementation of Betweenness algorithm according to *U. Brandes: On Variants of Shortest-Path Betweenness Centrality and their Generic Computation (2008)*.

## Compile instructions
```
module load intel/2017.00 CMake/3.5.2-intel-2017.00
cd Code
mkdir build && cd build
cmake ..
make
```
## Use instructions

1. Unzip experiment_graphs.zip
```
unzip Graphs/experiment_graphs.zip -d Graphs/
```


2. Run betweenness centrality algorithm
```
cd Code/build/
mpirun -n 24 ./betweenness -f ../../Graphs/graph-cze-brno.csv -v 2 -t 1
```


## Experiments

1. Sequential version on one node

```
qsub -q qexp -l select=1 -I
module load intel/2017.00
cd Code/build/
./betweenness -f ../../Graphs/graph-cze-brno.csv
```

File | Vertices | Edges | Time[s]
------------- |-------------|-------------|-------------
graph-cze-brno.csv    | 16901  | 36349  | 112.36
graph-cze-ostrava.csv | 14140  | 31494  | 78.13
graph-cze-praha.csv   | 48254  | 104033 | 1160.23
graph-prt-lisbon.csv  | 15381  | 30732  | 91.91
graph-prt-port.csv    | 7316   | 14793  | 17.70


2. Our parallel version on 10 nodes = 240 processes and 1 thread per process
```
qsub -q qexp -l select=10 -I
module load intel/2017.00
cd Code/build/
mpirun -n 240 ./betweenness -f ../../Graphs/BetweennessInputs/graph-cze-brno_preprocessed.csv -v 2 -t 1
```

File | Vertices | Edges | Time[s]
------------- |-------------|-------------|-------------
graph-cze.csv         | 969389 | 2287918 | 6434.95
graph-cze-brno.csv    | 16901  | 36349   | 1.36
graph-cze-ostrava.csv | 14140  | 31494   | 1.04
graph-cze-praha.csv   | 48254  | 104033  | 11.95
graph-prt.csv         | 942609 | 2265297 | 6430.89
graph-prt-lisbon.csv  | 15381  | 30732   | 1.20
graph-prt-port.csv    | 7316   | 14793   | 0.39

# LICENSE
This implementation of the Betweenness algorithm is licensed under the **GNU Lesser General Public License (LGPL)**. Full text of the LGPL can be found at [https://www.gnu.org/licenses/lgpl-3.0-standalone.html](https://www.gnu.org/licenses/lgpl-3.0-standalone.html) or [here](../LICENSE.LGPL.md).

This repository contains information from the Open Street Map, which is made available here under the **Open Database License (ODbL)**. Full text of the ODbL is available at [https://opendatacommons.org/licenses/odbl/1.0/](https://opendatacommons.org/licenses/odbl/1.0/) or [here](../LICENSE.ODBL.md).