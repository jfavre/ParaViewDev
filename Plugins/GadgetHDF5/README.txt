This ParaView plugin reads the HDF5 snapshots. A directory of snapshot files
is opened, and al HDF5 files are read. If a parallel ParaView server is used,
files will be spread across servers. If less files than servers are available,
some ParaView servers will be idle.

Contact: Jean M. Favre, Swiss National Supercomputing Center

Code is working with ParaView v5.4 and must be compiled (along with the main
ParaView source code) by setting the following cmake variable, e.g.:

PARAVIEW_EXTERNAL_PLUGIN_DIRS:STRING=/home/jfavre/Projects/Gadget

Last updated: Fri Aug 18 14:16:04 CEST 2017
