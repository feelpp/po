
if ( ACUSIM_FOUND )
  feelpp_add_application(pod SRCS main.cpp LINK_LIBRARIES ${PO_H3DREADER_LIB} ${PO_H3DWRITER_LIB} )
endif()
