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
feelpp_add_application( po_testcurl SRCS testcurl.cpp CFG po_testcurl.cfg )
feelpp_add_application( po_testgrad SRCS testgrad.cpp CFG po_testgrad.cfg )
feelpp_add_application( po_testD1_cyl SRCS testD1.cpp DEFS PO_CYL CFG po_testD1_cyl.cfg GEO cylCooCyl.geo  )
feelpp_add_application( po_testD1_car SRCS testD1.cpp CFG po_testD1_car.cfg GEO cylCooCar.geo  )
