# test code for libstokes
# Copyright (C) 2006-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
# $Id: test-stokes.rb,v 1.5 2007/03/08 00:50:58 kichiki Exp $
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

require 'stokes'

sys = Stokes::stokes_init()

np = 8
nm = 8
Stokes::stokes_set_np(sys, np, nm)
# you must call stokes_set_np() because
# this also allocate the memory for pos.

sys.periodic = 1 # periodic boundary condition
lx = 10.0
ly = 10.0
lz = 10.0
Stokes::stokes_set_l(sys, lx, ly, lz)
# you must call stokes_set_l() becuase
# this also initialize parameters other than lx,ly,lz.

tratio = 60.25
xi = Stokes::xi_by_tratio(sys, tratio)
cutlim = 1.0e-12
Stokes::stokes_set_xi(sys, xi, cutlim)

print "xi = ", xi, "\n"

sys.lubcut = 2.0000000001
Stokes::stokes_set_iter(sys, "gmres", 2000, 20, 1.0e-6,
                        1, Stokes::get_stdout())

pos = Stokes::Darray.new(np*3)
u   = Stokes::Darray.new(np*3)
f   = Stokes::Darray.new(np*3)

# set pos in SC lattice configuration
pos[ 0] = 0.0 # x component
pos[ 1] = 0.0 # y component
pos[ 2] = 0.0 # z component

pos[ 3] = 5.0
pos[ 4] = 0.0
pos[ 5] = 0.0

pos[ 6] = 0.0
pos[ 7] = 5.0
pos[ 8] = 0.0

pos[ 9] = 0.0
pos[10] = 0.0
pos[11] = 5.0

pos[12] = 5.0
pos[13] = 5.0
pos[14] = 0.0

pos[15] = 0.0
pos[16] = 5.0
pos[17] = 5.0

pos[18] = 5.0
pos[19] = 0.0
pos[20] = 5.0

pos[21] = 5.0
pos[22] = 5.0
pos[23] = 5.0

for i in 0..(np*3-1)
    u[i] = 1.0
end

print "pos:\n"
for i in 0..(np-1)
    print i, ' ', pos[i*3], ' ', pos[i*3+1], ' ', pos[i*3+2], "\n"
end

print "u:\n"
for i in 0..(np-1)
    print i, ' ', u[i*3], ' ', u[i*3+1], ' ', u[i*3+2], "\n"
end

#sys.pos = pos
Stokes::stokes_set_pos(sys, pos)
Stokes::solve_res_3f(sys, u, f)

nc_f = Stokes::stokes_nc_mob_f_init("test-stokes.res-3f.nc", np)
# f0, x, u are active
Stokes::stokes_nc_set_f0(nc_f, f)
Stokes::stokes_nc_set_time(nc_f, 0, 0.0)
Stokes::stokes_nc_set_x(nc_f, 0, pos)
Stokes::stokes_nc_set_u(nc_f, 0, u)

Stokes::stokes_nc_free(nc_f)


print "f:\n"
for i in 0..(np-1)
    print i, ' ', f[i*3], ' ', f[i*3+1], ' ', f[i*3+2], "\n"
end
