## sample code for stokes_res_ewald_3f()
## Copyright (C) 2006 Kengo Ichiki <kichiki@users.sourceforge.net>
## $Id: test-stokes.m,v 1.1 2006/10/04 20:39:16 ichiki Exp $
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

# lattice vector
l = [10.0, 10.0, 10.0]

# position of particles (np = 8)
pos = [0.0, 0.0, 0.0,\
	0.0, 0.0, 5.0,\
	0.0, 5.0, 0.0,\
	0.0, 5.0, 5.0,\
	5.0, 0.0, 0.0,\
	5.0, 0.0, 5.0,\
	5.0, 5.0, 0.0,\
	5.0, 5.0, 5.0]

# velocity of particles
u = [1.0, 1.0, 1.0,\
      1.0, 1.0, 1.0,\
      1.0, 1.0, 1.0,\
      1.0, 1.0, 1.0,\
      1.0, 1.0, 1.0,\
      1.0, 1.0, 1.0,\
      1.0, 1.0, 1.0,\
      1.0, 1.0, 1.0]

# solve resistance problem and obtain force
f = stokes_res_ewald_3f(pos, u, l, 60.25)
