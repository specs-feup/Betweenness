module load intel/2017.00
module load CMake/3.9.0

cd Code
rm -r build
mkdir build
cd build
cmake ..
make

cd ../../
unzip Graphs/experiment_graphs.zip -d Graphs/