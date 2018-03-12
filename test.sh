sh compile.sh

#load necessary modules
module load intel/2017.00

#unzip OpenStreetMaps
unzip Graphs/experiment_graphs.zip -d Graphs/

#run BC algorithm on generated files
mpirun -n 24 ./Code/build/betweenness -f Graphs/graph-prt-port.csv -v 2 -t 1

#clean
sh clean.sh