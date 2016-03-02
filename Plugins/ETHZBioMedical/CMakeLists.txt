
vtk_module_load(vtkhdf5)

INCLUDE_DIRECTORIES(
  ${vtkhdf5_INCLUDE_DIRS}
  ${VTK_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/aimpack64
  )

ADD_PARAVIEW_PLUGIN(ETHZBioMedical "1.0"
  SERVER_MANAGER_SOURCES
    vtkETHZBioMedicalReader.cxx
    vtkETHZBioMedicalHDF5Reader.cxx
    vtkETHZBioMedicalMaskHDF5Reader.cxx
    vtkETHZBioMedicalFAIMReader.cxx
    vtkAIMReader.cxx

  SERVER_SOURCES
    BaseETHZBioMedicalHDF5Reader.cxx
    aimpack64/aimpack.c

  SERVER_MANAGER_XML
    BioMedical_Server.xml

  REQUIRED_ON_SERVER
)

target_link_libraries(ETHZBioMedical LINK_PRIVATE
  ${vtkhdf5_LIBRARIES}
  vtksys)

