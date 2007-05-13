# stokes-netcdf to pov converter
# Copyright (C) 2006-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
# $Id: stnc2pov.py,v 1.5 2007/05/13 23:14:10 kichiki Exp $
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


def write_T_Particle(f):
    f.write('#declare T_Particle = texture {\n'\
            '  pigment { color White }\n'\
            '  //finish { ambient 0.2 diffuse 0 reflection 0.6 }\n'\
            '  finish {\n'\
            '    ambient .1\n'\
            '    diffuse .1\n'\
            '    specular 1\n'\
            '    roughness .001\n'\
            '    metallic\n'\
            '    reflection {\n'\
            '      .75\n'\
            '      metallic\n'\
            '    }\n'\
            '  }\n'\
            '}\n')

def write_M_RYUON(f):
    # M_RYUON
    f.write('#declare M_RYUON = material {\n'\
            '  texture {\n'\
            '    pigment {\n'\
            '      color <0.4, 0.5, 1.0>\n'\
            '      filter 1\n'\
            '    }\n'\
            '    finish {\n'\
            '      ambient 0\n'\
            '      diffuse 0\n'\
            '      reflection .25\n'\
            '      specular 1\n'\
            '      roughness .001\n'\
            '    }\n'\
            '  } // end of texture\n'\
            '  interior { ior 1.33 }\n'
            '}\n')

def write_T_CHECKER(f):
    f.write('#declare T_CHECKER = texture {\n'\
            '  pigment{\n'\
            '    checker\n'\
            '    color <0.4, 0.5, 1.0>\n'\
            '    color White\n'\
            '  }\n'\
            '  scale 0.01\n'\
            '  finish{\n'\
            '    phong 0.9\n'\
            '    metallic\n'\
            '  }\n'\
            '}\n')


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

    f.write('#include "colors.inc"\n')
    #f.write('#include "woods.inc"\n\n')
    # place the ground
    f.write('// floor\nplane {\n'\
            '  y, -0.1\n'\
            '  texture {\n'\
            #'    T_Wood6\n'\
            #'    finish{ ambient 1 }\n'\
            '    pigment { checker color White, color <.7, .7, .7> }\n'\
            '    scale .3\n'\
            '    finish{ ambient .4 }\n'\
            '  }\n'\
            '}\n')
    # place the walls
    f.write('// back wall\n'\
            'plane {\n'\
            '  z, 1\n'\
            '  pigment { color rgb <1,1,0.8> }\n'\
            '  finish{ ambient 0.4 }\n'\
            '}\n')
    f.write('// ceiling\n'\
            'plane {\n'\
            '  y, 5\n'\
            '  pigment { color White }\n'\
            '}\n')
    f.write('// right wall\n'\
            'plane {\n'\
            '  x, 5\n'\
            '  pigment { color White }\n'\
            '}\n')
    f.write('// left wall\n'\
            'plane {\n'\
            '  x, -5\n'\
            '  pigment { color White }\n'\
            '}\n')
    f.write('// behind wall\n'\
            'plane {\n  z, -5\n'\
            '  pigment { color White }\n'\
            '}\n\n')
    # place the box
    f.write('box {\n'\
            '  <0, 0, 0>,  // Near lower left corner\n'\
            '  <%f, %f, %f>   // Far upper right corner\n'\
            '  pigment { color rgbf <0.9, 0.99, 1, 1> }\n'\
            '}\n\n'%(lx, ly, lz))

    f.write('camera {\n  location <%f, %f, %f>\n'%(cx, cy, cz))
    f.write('  look_at  <%f, %f,  %f>\n}\n\n'%(lax, lay, laz))
    f.write('light_source { <2, 4.9, -3> color White}\n\n')

    write_T_Particle(f)
    write_M_RYUON (f)
    write_T_CHECKER(f)


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
    f.write('// back wall\n'\
            'plane {\n'\
            '  z, 2\n'\
            '  pigment { checker color White, color <0.6, 0.8, 1> }\n'\
            '  scale 0.1\n}\n')
    f.write('// behind wall\n'\
            'plane {\n'\
            '  z, -5\n'\
            '  pigment { color White }\n'\
            '}\n\n')
    #
    f.write('camera {\n  location <%f, %f, %f>\n'%(cx, cy, cz))
    f.write('  look_at  <%f, %f,  %f>\n}\n\n'%(lax, lay, laz))
    f.write('light_source { <2, 4.9, -3> color White}\n\n')

    write_T_Particle(f)
    write_M_RYUON (f)
    write_T_CHECKER(f)


def write_pov_particle (f, x, y, z, a):
    # note that in POVRAY,
    # y is the vertical direction
    # z is the depth direction
    # scale factor = 1/100 (0.01 radius = 1 in POV)
    f.write('sphere {\n')
    f.write('  <%f, %f, %f>, %f\n'%(x/100.0, z/100.0, y/100.0, a/100.0))
    f.write('  material { M_RYUON }\n}\n')

def write_pov_particle_fixed (f, x, y, z, a):
    # note that in POVRAY,
    # y is the vertical direction
    # z is the depth direction
    # scale factor = 1/100 (0.01 radius = 1 in POV)
    f.write('sphere {\n')
    f.write('  <%f, %f, %f>, %f\n'%(x/100.0, z/100.0, y/100.0, a/100.0))
    f.write('  texture { T_Particle }\n}\n')


# make transform matrix (3x3) by quaternion
def Q2M (q1,q2,q3,q4):
    m = []

    # parity change
    q4 *= -1.0

    m.append(2.0*(q1*q1 + q4*q4 - .5))
    m.append(2.0*(q1*q2 + q3*q4))
    m.append(2.0*(q1*q3 - q2*q4))

    m.append(2.0*(q1*q2 - q3*q4))
    m.append(2.0*(q2*q2 + q4*q4 - .5))
    m.append(2.0*(q2*q3 + q1*q4))

    m.append(2.0*(q1*q3 + q2*q4))
    m.append(2.0*(q2*q3 - q1*q4))
    m.append(2.0*(q3*q3 + q4*q4 - .5))

    # note that in POVRAY,
    # y is the vertical direction
    # z is the depth direction
    t = [1.0, 0.0, 0.0,\
         0.0, 0.0, 1.0,\
         0.0, 1.0, 0.0]
    x = [0.0, 0.0, 0.0,  0.0, 0.0, 0.0,  0.0, 0.0, 0.0]
    y = [0.0, 0.0, 0.0,  0.0, 0.0, 0.0,  0.0, 0.0, 0.0]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                x[i*3+j] += t[i*3+k] * m[k*3+j]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                y[i*3+j] += x[i*3+k] * t[k*3+j]
    # therefore, Y = T . M . T
    return y

def write_pov_particle_Q (f, x, y, z, a, q):
    # note that in POVRAY,
    # y is the vertical direction
    # z is the depth direction
    # scale factor = 1/100 (0.01 radius = 1 in POV)
    m = Q2M (q[0], q[1], q[2], q[3])
    f.write('sphere {\n')
    f.write('  <0, 0, 0>, %f\n'%(a/100.0))
    f.write('  texture { T_CHECKER }\n')
    f.write('  transform {\n')
    f.write('    matrix <%f, %f, %f,\n'%(m[0], m[3], m[6]))
    f.write('            %f, %f, %f,\n'%(m[1], m[4], m[7]))
    f.write('            %f, %f, %f,\n'%(m[2], m[5], m[8]))
    f.write('            %f, %f, %f> }\n}\n'%(x/100.0, z/100.0, y/100.0))
    

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
    print '$Id: stnc2pov.py,v 1.5 2007/05/13 23:14:10 kichiki Exp $'
    print 'USAGE:'
    print '\t-f or --file : stokes-nc-file'
    sys.exit ()


def main():
    filename = ''
    nm = 0
    i = 1
    while i < len(sys.argv):
        if sys.argv[i] == '-f' or sys.argv[i] == '--file':
            filename = sys.argv[i+1]
            i += 2
        else:
            usage()
    if filename == '': usage()

    nc = stokes.stokes_nc_open (filename)
    #stokes.stokes_nc_print_actives(nc, stokes.get_stdout())

    lattice = stokes.darray(3)
    stokes.stokes_nc_get_l (nc, lattice)

    # x[] : center of particles
    pos  = stokes.darray(nc.np  * nc.nvec)
    # q[] : quaternion
    if nc.flag_q != 0:
        q  = stokes.darray(nc.np  * nc.nquat)
    else:
        q = []
    # a[] : radius of mobile particles
    if nc.flag_a != 0:
        a = stokes.darray(nc.np)
        stokes.stokes_nc_get_data0 (nc, "a", a)
    else:
        a = []
    # af[] : radius of fixed particles
    if nc.flag_af != 0:
        af = stokes.darray(nc.npf)
        stokes.stokes_nc_get_data0 (nc, "af", af)
    else:
        af = []

    # xf0[]
    if nc.npf > 0:
        xf0 = stokes.darray(nc.npf * nc.nvec)
        stokes.stokes_nc_get_data0 (nc, "xf0", xf0)

    if lattice[0] != 0.0 or lattice[1] != 0.0 or lattice[2] != 0.0:
        # periodic boundary
        if lattice[0] > lattice[2]:
            l = lattice[0]
        else:                  
            l = lattice[2]
        #camera = [0.5 * lattice[0], -110, 0.5 * lattice[2]]
        camera = [0.5 * lattice[0], -1.7*l, 0.5 * lattice[2]]
        lookat = [0.5 * lattice[0],  0.0, 0.5 * lattice[2]]

    for i in range(nc.ntime):
        file = 'test%04d.pov'%i
        try:
            f = open(file, 'w')
        except IOError:
            print 'cannot open', file
            sys.exit()

        stokes.stokes_nc_get_data (nc, "x", i, pos)

        if lattice[0] == 0.0 and lattice[1] == 0.0 and lattice[2] == 0.0:
            # non-periodic boundary
            (cx,cy,cz, lx,ly,lz) = bounding_box (nc.np, pos)
            if lx > lz:
                l = lx
            else:                  
                l = lz
            # prevent to go far away
            if l > 50: l = 50

            camera = [cx, cy-2*l, cz]
            lookat = [cx, cy,   cz]
            write_pov_header_open (f, lattice, camera, lookat)

        else:
            # periodic boundary
            #move_camera (camera, lookat)
            write_pov_header (f, lattice, camera, lookat)

        if nc.flag_q != 0:
            # with quaternion
            stokes.stokes_nc_get_data (nc, "q", i, q)
            for j in range(nc.np):
                x = pos[j*3]
                y = pos[j*3+1]
                z = pos[j*3+2]
                if a != []:
                    write_pov_particle_Q (f, x, y, z, a[j],\
                                          [q[j*4+0],q[j*4+1],\
                                           q[j*4+2],q[j*4+3]])
                else:
                    write_pov_particle_Q (f, x, y, z, 1.0,\
                                          [q[j*4+0],q[j*4+1],\
                                           q[j*4+2],q[j*4+3]])
        else:
            # no quaternion
            for j in range(nc.np):
                x = pos[j*3]
                y = pos[j*3+1]
                z = pos[j*3+2]
                if a != []:
                    write_pov_particle (f, x, y, z, a[j])
                else:
                    write_pov_particle (f, x, y, z, 1.0)

        # fixed particles
        for j in range(nc.npf):
            x = xf0[j*3]
            y = xf0[j*3+1]
            z = xf0[j*3+2]
            if af != []:
                write_pov_particle_fixed (f, x, y, z, af[j])
            else:
                write_pov_particle_fixed (f, x, y, z, 1.0)

        # done
        f.close()


if __name__ == "__main__":
    main()
