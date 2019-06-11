
paraview_add_plugin(pvNek5000Reader
  REQUIRED_ON_SERVER
  REQUIRED_ON_CLIENT
  VERSION "1.0"
  MODULES Nek5000Reader
  MODULE_FILES      "${CMAKE_CURRENT_SOURCE_DIR}/Reader/vtk.module")

#if (BUILD_TESTING)
#  add_subdirectory(Testing)
#endif ()

