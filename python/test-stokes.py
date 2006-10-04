# test code for libstokes
# Copyright (C) 2006 Kengo Ichiki <kichiki@users.sourceforge.net>
# $Id: test-stokes.py,v 1.1 2006/10/03 21:44:01 ichiki Exp $
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
import stokes

sys = stokes.stokes_init()

np = 8
nm = 8
stokes.stokes_set_np(sys, np, nm)

lx = 10.0
ly = 10.0
lz = 10.0
stokes.stokes_set_ll(sys, lx, ly, lz)

tratio = 60.25
zeta = stokes.zeta_by_tratio(sys, tratio)
cutlim = 1.0e-12
stokes.stokes_set_zeta(sys, zeta, cutlim)

print 'zeta =', zeta

sys.lubcut = 2.0000000001
sys.it = stokes.iter_init ("gmres", 2000, 20, 1.0e-6, 1)

pos = stokes.darray(np*3)
u   = stokes.darray(np*3)
f   = stokes.darray(np*3)

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

for i in range(np*3):
    u[i] = 1.0

print 'pos:'
for i in range(np):
    print i,pos[i*3],pos[i*3+1],pos[i*3+2]

print 'u:'
for i in range(np):
    print i, u[i*3], u[i*3+1], u[i*3+2]

sys.pos = pos

stokes.calc_res_ewald_3f(sys, u, f)

print 'f:'
for i in range(np):
    print i, f[i*3], f[i*3+1], f[i*3+2]