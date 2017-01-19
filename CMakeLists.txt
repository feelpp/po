cmake_minimum_required(VERSION 2.8)
###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
#       Date: 2013-10-14
#
#  Copyright (C) 2013
#
# Distributed under the GPL(GNU Public License):
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
#
project(po)

if ( ${CMAKE_MAJOR_VERSION} EQUAL 3 AND ${CMAKE_MINOR_VERSION} GREATER 3 )
  message(STATUS "[feelpp] use OLD policy CMP0064" )
  cmake_policy(SET CMP0064 OLD)
endif()

# with this it is possible to compile applications/models from an installed
# version of feel++. This is useful also for staged compilation for docker
if ( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} )
  find_package(Feel++ PATHS $ENV{FEELPP_DIR}/share/feel/cmake/modules /usr/share/feel/cmake/modules /usr/local/share/feel/cmake/modules )
  if(NOT FEELPP_FOUND)
      message(FATAL_ERROR "Feel++ was not found on your system. Make sure to install it and specify the FEELPP_DIR to reference the installation directory.")
  endif()
endif()

find_library( PO_H3DREADER_LIB NAMES h3dreader PATHS ${CMAKE_CURRENT_SOURCE_DIR}/lib/h3d/bin/linux64 NO_DEFAULT_PATH )
message(STATUS "[feelpp/po] PO_H3DREADER_LIBRARY: ${PO_H3DREADER_LIB}")
find_library( PO_H3DWRITER_LIB NAMES h3dwriter PATHS ${CMAKE_CURRENT_SOURCE_DIR}/lib/h3d/bin/linux64 NO_DEFAULT_PATH )
message(STATUS "[feelpp/po] PO_H3DWRITER_LIBRARY: ${PO_H3DWRITER_LIB}")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/h3d/include)

#set(FEELPP_LIBRARIES "${PO_H3DWRITER_LIB} ${PO_H3DREADER_LIB} ${FEELPP_LIBRARIES}")

add_subdirectory( lib/h3d/examples )
add_subdirectory( offline )
add_subdirectory( online )
add_subdirectory( test )
add_subdirectory( data_conversion )
