set(classes
  vtkHDF5Reader)

set(private_headers
  vtkHDF5Reader.h)
  
vtk_module_add_module(HDF5
  CLASSES ${classes}
  PRIVATE_HEADERS ${private_headers})

#vtk_module_include(HDF5
#  PRIVATE
#  /local/apps/ParaView/dev/include/paraview-5.6/vtkhdf5
#  /local/apps/ParaView/dev/include/paraview-5.6/vtkhdf5/src)
  
paraview_add_server_manager_xmls(
	XMLS  HDF5.xml)
