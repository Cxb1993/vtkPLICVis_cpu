CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

IF(COMMAND cmake_policy)
  CMAKE_POLICY(SET CMP0003 NEW)
  CMAKE_POLICY(SET CMP0023 NEW)
ENDIF(COMMAND cmake_policy)

IF (NOT ParaView_BINARY_DIR)
  FIND_PACKAGE(ParaView REQUIRED)
  INCLUDE(${PARAVIEW_USE_FILE})
ENDIF()

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -std=c++11")

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

ADD_LIBRARY(PLICVis_impl
  SHARED
  plicvis_impl.cxx)

ADD_PARAVIEW_PLUGIN(PLICVis "1.0"
  SERVER_MANAGER_XML PLICVis.xml
  SERVER_MANAGER_SOURCES vtkPLICVis.cxx
  GUI_RESOURCES PLICVis.qrc)

TARGET_LINK_LIBRARIES(PLICVis PRIVATE PLICVis_impl)
