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
feelpp_add_application( po_test SRCS test.cpp CFG po_test.cfg )
feelpp_add_application( po_test_curl SRCS test_curl.cpp CFG po_test_curl.cfg )
feelpp_add_application( po_test_grad SRCS test_grad.cpp CFG po_test_grad.cfg )
feelpp_add_application( po_test_cyl SRCS test_cyl.cpp DEFS PO_CYL CFG po_test_cyl.cfg GEO ../mesh/cylCooCyl.geo  )
feelpp_add_application( po_test_car SRCS test_cyl.cpp CFG po_test_car.cfg GEO ../mesh/cylCooCar.geo  )
feelpp_add_application( po_test_sph_sph SRCS test_sph.cpp DEFS PO_SPH CFG po_test_sph_sph.cfg GEO ../mesh/sphCooSph.geo  )
feelpp_add_application( po_test_sph_car SRCS test_sph.cpp CFG po_test_sph_car.cfg GEO ../mesh/sphCooCar.geo  )
feelpp_add_application( po_test_spheromak SRCS spheromak.cpp CFG po_spheromak.cfg GEO ../mesh/shell.geo )
feelpp_add_application( po_test_eigencube SRCS eigencube.cpp CFG po_eigencube.cfg GEO ../mesh/cubepi.geo )
feelpp_add_application( po_eigenmixed_curl SRCS eigenmixed.cpp CFG po_eigenmixed.cfg DEFS PO_UHCURL )
feelpp_add_application( po_eigenmixed_div SRCS eigenmixed.cpp CFG po_eigenmixed.cfg )
feelpp_add_application( po_stokesmixed SRCS stokesmixed.cpp CFG po_stokesmixed.cfg )
