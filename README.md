# projektSRiR

```bash
source /opt/nfs/config/source_mpich430.sh
source /opt/nfs/config/source_cuda121.sh
source /opt/nfs/config/source_gaspi.sh
mkdir build
cd build
cmake ..
make
cd ..
gaspi_run -m nodes -n 4 $GPI_BIN_DIR/gaspi_wrapper.sh $(pwd)/build/bellman_ford_gaspi ./input/input_10.txt 2 1
