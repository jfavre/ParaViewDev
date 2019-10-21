set(classes
  vtkNek5000Reader)

set(private_headers
  vtkNek5000Reader.h)

vtk_module_add_module(Nek5000Reader
  CLASSES ${classes}
  PRIVATE_HEADERS ${private_headers})

paraview_add_server_manager_xmls(
  XMLS  Nek5000Reader.xml)

ADD_EXECUTABLE(TestReaderSimple TestReaderSimple.cxx)
TARGET_LINK_LIBRARIES(TestReaderSimple Nek5000Reader)

