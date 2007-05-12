# stokes-netcdf to pov converter
# Copyright (C) 2006-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
# $Id: stnc2pov.py,v 1.4 2007/05/12 04:37:47 kichiki Exp $
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
import sys
#sys.path.append('/somewhere/ryuon/stokes/python')
import stokes


def bounding_box (np, x):
    (cx,cy,cz) = (0.0, 0.0, 0.0)
    for i in range(np):
        xx = x[i*3]
        yy = x[i*3+1]
        zz = x[i*3+2]

        cx = cx + xx
        cy = cy + yy
        cz = cz + zz

        if i == 0:
            lx0 = lx1 = xx
            ly0 = ly1 = yy
            lz0 = lz1 = zz
        else:
            if lx0 > xx:
                lx0 = xx
            if lx1 < xx:
                lx1 = xx
            if ly0 > yy:
                ly0 = yy
            if ly1 < yy:
                ly1 = yy
            if lz0 > zz:
                lz0 = zz
            if lz1 < zz:
                lz1 = zz
        
    cx = cx / float(np)
    cy = cy / float(np)
    cz = cz / float(np)

    lx = lx1 - lx0
    ly = ly1 - ly0
    lz = lz1 - lz0

    return (cx,cy,cz, lx,ly,lz)


# INPUT
#  f : file
#  lattice = (lx,ly,lz) in simulation coordinates
#  camera = (cx,cy,cz)  
#  lookat = (lax,lay,laz) 
def write_pov_header (f, lattice, camera, lookat):
    # note that in POVRAY,
    # y is the vertical direction
    # z is the depth direction
    # scale factor = 1/100 (0.01 radius = 1 in POV)
    lx = lattice[0]/100.0
    lz = lattice[1]/100.0
    ly = lattice[2]/100.0

    cx = camera[0]/100.0
    cz = camera[1]/100.0
    cy = camera[2]/100.0

    lax = lookat[0]/100.0
    laz = lookat[1]/100.0
    lay = lookat[2]/100.0

    f.write('#include \"colors.inc\"\n')
    f.write('#include "woods.inc"\n\n')
    # place the ground
    f.write('// floor\nplane {\n  y, 0\n  texture { T_Wood6 }\n}\n')
    # place the walls
    f.write('// back wall\nplane {\n  z, 1\n  pigment { color rgb <1,1,0.8> }\n}\n')
    f.write('// ceiling\nplane {\n  y, 5\n  pigment { color White }\n}\n')
    f.write('// right wall\nplane {\n  x, 5\n  pigment { color White }\n}\n')
    f.write('// left wall\nplane {\n  x, -5\n  pigment { color White }\n}\n')
    f.write('// behind wall\nplane {\n  z, -5\n  pigment { color White }\n}\n\n')
    # place the box
    f.write('box {\n  <0, 0, 0>,  // Near lower left corner\n  <%f, %f, %f>   // Far upper right corner\n  pigment { color rgbf <0.9, 0.99, 1, 1> }\n}\n\n'%(lx, ly, lz))
    #
    f.write('camera {\n  location <%f, %f, %f>\n'%(cx, cy, cz))
    f.write('  look_at  <%f, %f,  %f>\n}\n\n'%(lax, lay, laz))
    f.write('light_source { <2, 4.9, -3> color White}\n\n')
    f.write('#declare T_Particle = texture {\n  pigment { color White }\n  //finish { ambient 0.2 diffuse 0 reflection 0.6 }\n  finish {\n    ambient .1\n    diffuse .1\n    specular 1\n    roughness .001\n    metallic\n    reflection {\n      .75\n      metallic\n    }\n  }\n}')

def write_pov_particle (f, x, y, z, a):
    # note that in POVRAY,
    # y is the vertical direction
    # z is the depth direction
    # scale factor = 1/100 (0.01 radius = 1 in POV)
    f.write('sphere {\n')
    f.write('  <%f, %f, %f>, %f\n'%(x/100.0, z/100.0, y/100.0, a/100.0))
    f.write('  texture { T_Particle }\n}\n');


# INPUT
#  f : file
#  lattice = (lx,ly,lz) in simulation coordinates
#  camera = (cx,cy,cz)  
#  lookat = (lax,lay,laz) 
def write_pov_header_open (f, lattice, camera, lookat):
    # note that in POVRAY,
    # y is the vertical direction
    # z is the depth direction
    # scale factor = 1/100 (0.01 radius = 1 in POV)
    lx = lattice[0]/100.0
    lz = lattice[1]/100.0
    ly = lattice[2]/100.0

    cx = camera[0]/100.0
    cz = camera[1]/100.0
    cy = camera[2]/100.0

    lax = lookat[0]/100.0
    laz = lookat[1]/100.0
    lay = lookat[2]/100.0

    f.write('#include \"colors.inc\"\n')
    f.write('#include "woods.inc"\n\n')
    # place the walls
    f.write('// back wall\nplane {\n  z, 2\n  pigment { checker color White, color <0.6, 0.8, 1> }\n  scale 0.1\n}\n')
    f.write('// behind wall\nplane {\n  z, -5\n  pigment { color White }\n}\n\n')
    #
    f.write('camera {\n  location <%f, %f, %f>\n'%(cx, cy, cz))
    f.write('  look_at  <%f, %f,  %f>\n}\n\n'%(lax, lay, laz))
    f.write('light_source { <2, 4.9, -3> color White}\n\n')
    f.write('#declare T_Particle = texture {\n  pigment { color White }\n  //finish { ambient 0.2 diffuse 0 reflection 0.6 }\n  finish {\n    ambient .1\n    diffuse .1\n    specular 1\n    roughness .001\n    metallic\n    reflection {\n      .75\n      metallic\n    }\n  }\n}')
    # M_RYUON
    f.write('#declare M_RYUON = material {\n  texture {\n    pigment {\n      color <0.4, 0.5, 1.0>\n      filter 1\n    }\n    finish {\n      ambient 0\n      diffuse 0\n      reflection .25\n      specular 1\n      roughness .001\n    }\n  } // end of texture\n  interior { ior 1.33 }\n}\n')

def write_pov_particle (f, x, y, z, a):
    # note that in POVRAY,
    # y is the vertical direction
    # z is the depth direction
    # scale factor = 1/100 (0.01 radius = 1 in POV)
    f.write('sphere {\n')
    f.write('  <%f, %f, %f>, %f\n'%(x/100.0, z/100.0, y/100.0, a/100.0))
    #f.write('  texture { T_Particle }\n}\n')
    f.write('  material { M_RYUON }\n}\n')

def write_pov_particle_fixed (f, x, y, z, a):
    # note that in POVRAY,
    # y is the vertical direction
    # z is the depth direction
    # scale factor = 1/100 (0.01 radius = 1 in POV)
    f.write('sphere {\n')
    f.write('  <%f, %f, %f>, %f\n'%(x/100.0, z/100.0, y/100.0, a/100.0))
    f.write('  texture { T_Particle }\n}\n')
    #f.write('  material { M_RYUON }\n}\n')


# now camera angle
# init:   <0.17,  0.50, -1.10>  <0.17,  0.50, 0.0>
# target: <0.17,  0.22, -0.28>  <0.17,  0.15, 0.0>
# those are in POVRAY coordinates
# in simulation coordinates,
# init:   <17, -110,  50>  <17, 0,  50>
# target: <17,  -28,  22>  <17, 0,  15>
# diff    < 0,   82, -28>  < 0, 0, -35>
# let's say targe is reached in the first 200 steps
# d/step= < 0,  .41,-.14>  < 0, 0,-.175>
def move_camera (camera, lookat):
    if (camera[2] <= 22.0): return
    camera[1] += .41
    camera[2] -= .14

    lookat[2] -= .175


def usage():
    print '$Id: stnc2pov.py,v 1.4 2007/05/12 04:37:47 kichiki Exp $'
    print 'USAGE:'
    print '\t-f or --file : stokes-nc-file'
    print '\t-nm          : number of mobile particles (default: all particles)'
    sys.exit ()


def main():
    filename = ''
    nm = 0
    i = 1
    while i < len(sys.argv):
        if sys.argv[i] == '-f' or sys.argv[i] == '--file':
            filename = sys.argv[i+1]
            i += 2
        elif sys.argv[i] == '-nm':
            nm = int(sys.argv[i+1])
            i += 2
        else:
            usage()
    if filename == '': usage()

    nc = stokes.stokes_nc_open (filename)
    #stokes.stokes_nc_print_actives(nc, stokes.get_stdout())
    if nm == 0:
        nm = nc.np
    lattice = stokes.darray(3)
    stokes.stokes_nc_get_l (nc, lattice)

    pos  = stokes.darray(nc.np  * nc.nvec)
    pos0 = stokes.darray(nc.npf * nc.nvec)
    if nc.flag_a != 0:
        a = stokes.darray(nc.np)
        stokes.stokes_nc_get_data0 (nc, "a", a)
    else:
        a = []
    if nc.flag_af != 0:
        af = stokes.darray(nc.npf)
        stokes.stokes_nc_get_data0 (nc, "af", af)
    else:
        af = []

    if lattice[0] != 0.0 or lattice[1] != 0.0 or lattice[2] != 0.0:
        # periodic boundary
        camera = [0.5 * lattice[0], -110, 0.5 * lattice[2]]
        lookat = [0.5 * lattice[0],  0.0, 0.5 * lattice[2]]

    i = 0
    while 1:
        file = 'test%04d.pov'%i
        try:
            f = open(file, 'w')
        except IOError:
            print 'cannot open', file
            sys.exit()

        stokes.stokes_nc_get_data (nc, "x", i, pos)
        stokes.stokes_nc_get_data0 (nc, "xf0", pos0)

        if lattice[0] == 0.0 and lattice[1] == 0.0 and lattice[2] == 0.0:
            # non-periodic boundary
            (cx,cy,cz, lx,ly,lz) = bounding_box (nm, pos)
            if lx > lz:
                l = lx
            else:                  
                l = lz
            camera = [cx, cy-2*l, cz]
            lookat = [cx, cy,   cz]
            write_pov_header_open (f, lattice, camera, lookat)
        else:
            # periodic boundary
            move_camera (camera, lookat)
            write_pov_header (f, lattice, camera, lookat)

        for j in range(nm):
            x = pos[j*3]
            y = pos[j*3+1]
            z = pos[j*3+2]
            write_pov_particle (f, x, y, z, 1.0)
        # in case nc.np contains both mobile and fixed
        #print '# np = %d, nm = %d'%(nm, nc.np)
        for j in range(nm, nc.np):
            x = pos[j*3]
            y = pos[j*3+1]
            z = pos[j*3+2]
            if a != []:
                write_pov_particle_fixed (f, x, y, z, a[j])
            else:
                write_pov_particle_fixed (f, x, y, z, 1.0)
        for j in range(nc.npf):
            x = pos0[j*3]
            y = pos0[j*3+1]
            z = pos0[j*3+2]
            if af != []:
                write_pov_particle_fixed (f, x, y, z, af[j])
            else:
                write_pov_particle_fixed (f, x, y, z, 1.0)
        f.close()

        i += 1


if __name__ == "__main__":
    main()
