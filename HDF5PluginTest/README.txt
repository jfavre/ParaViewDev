# test HDF5 reader plugin for ParaView > v5.6

cd build
cmake ..
make && paraview
Load new plugin
try opening ExternalData/Testing/Data/Test1.h5 in your ParaView compile directory
