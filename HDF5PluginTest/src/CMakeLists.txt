paraview_add_plugin(pvHDF5Reader
  VERSION "1.0"
  REQUIRED_ON_SERVER
  MODULES HDF5
  MODULE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/Reader/vtk.module")
