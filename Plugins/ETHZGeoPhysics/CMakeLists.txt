cmake_minimum_required(VERSION 2.8)

SET(PLUGIN_NAME ETHZGeoPhysics)
SET(PLUGIN_VERSION "1.0")

FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})

#IF (QT_USE_FILE)
#  INCLUDE(${QT_USE_FILE})
#ENDIF (QT_USE_FILE)

INCLUDE_DIRECTORIES(
  ${VTK_INCLUDE_DIRS}
  ${vtkhdf5_INCLUDE_DIRS}
  )

ADD_PARAVIEW_PLUGIN(${PLUGIN_NAME} ${PLUGIN_VERSION}
    SERVER_MANAGER_XML CSCS_ALPSJackson_Server.xml
    SOURCES
                    BaseGeoPhysicsHDF5SphericalReader.cxx
    SERVER_MANAGER_SOURCES
                    vtkETHZGeoPhysicsHDF5Reader.cxx
                    vtkETHZGeoPhysicsHDF5SphericalReader.cxx
#                    vtkETHZGeoPhysicsHDF5Writer.cxx
#                    vtkETHZGeoPhysicsDerivedQuantities.cxx
#                    vtkETHZGeoPhysicsLongitudinalAverage.cxx
#		    vtkETHZGeoPhysicsVelPhiSGrid.cxx
                    vtkETHZGeoPhysicsSphericalLongitudinalAverage.cxx
                    vtkMakeGeoEarth.cxx
#                    vtkSG2RG.cxx
                    vtkProjMap.cxx
#    GUI_RESOURCE_FILES ETHZGeoPhysics-GUI.xml
    REQUIRED_ON_SERVER
)

vtk_module_load(vtklibproj4)

target_link_libraries(${PLUGIN_NAME} LINK_PRIVATE
  ${vtkhdf5_LIBRARIES} vtkGeovisCore
  vtksys)

ADD_EXECUTABLE(TestSphericalReader TestSphericalReader.cxx)
TARGET_LINK_LIBRARIES(TestSphericalReader vtksys vtkCommonCore vtkRenderingCore ${PLUGIN_NAME})

ADD_EXECUTABLE(TestLeaks TestLeaks.cxx)
TARGET_LINK_LIBRARIES(TestLeaks vtkCommonCore vtkCommonDataModel vtkCommonExecutionModel)

#SUBDIRS(Testing)

#ADD_EXECUTABLE(sp2xyz sp2xyz.cxx)
#TARGET_LINK_LIBRARIES(sp2xyz vtkIO vtkCommon vtkGraphics )

#ADD_EXECUTABLE(SpaceTimePlot SpaceTimePlot.cxx)  
#TARGET_LINK_LIBRARIES(SpaceTimePlot vtkIOXML)
