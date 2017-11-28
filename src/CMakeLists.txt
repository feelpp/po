
# with this it is possible to compile applications/models from an installed
# version of feel++. This is useful also for staged compilation for docker
if ( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} )
  find_package(Feel++ PATHS $ENV{FEELPP_DIR}/share/feelpp/feel/cmake/modules /usr/share/feelpp/feel/cmake/modules /usr/local/share/feelpp/feel/cmake/modules )
  if(NOT FEELPP_FOUND)
      message(FATAL_ERROR "Feel++ was not found on your system. Make sure to install it and specify the FEELPP_DIR to reference the installation directory.")
  endif()
endif()


# add your stuff here
feelpp_add_application( converter SRCS converter_database_acusim_new.cpp )
feelpp_add_application( converter_old SRCS converter_database_acusim_old.cpp )
feelpp_add_application( pod SRCS pod.cpp )
feelpp_add_application( dmd SRCS dmd.cpp )
feelpp_add_application( svd_dmd SRCS svd_dmd.cpp )
feelpp_add_application( Qcriterion SRCS Qcriterion_update.cpp )
feelpp_add_application( spdmd SRCS spdmd.cpp )
feelpp_add_application( spdmd_rec SRCS spdmd_rec.cpp )

