cmake_minimum_required( VERSION 2.6 )

set( CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE ${CMAKE_SOURCE_DIR}/../../libs )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE ${CMAKE_SOURCE_DIR}/../../libs )

set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${CMAKE_SOURCE_DIR}/../../libs )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${CMAKE_SOURCE_DIR}/../../libs )

project( mesh_signed_distance )

find_package( CGAL REQUIRED )
include( ${CGAL_USE_FILE} )

include_directories( . ${CMAKE_SOURCE_DIR}/../../include )

add_library( mesh_signed_distance SHARED mesh_signed_distance.cpp mesh_signed_distance_cgal.h )
