# configuration generating script
# Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
# $Id: gen-config.py,v 1.2 2007/11/18 05:20:20 kichiki Exp $
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
import math
import random


# OUTPUT
#  0 : no-overlapping
#  1 : overlapping
def check_overlap (pos, x, y, z, a2):
    for i in range(len(pos)/3):
        ix = i*3
        iy = ix+1
        iz = ix+2
        dx = x - pos[ix]
        dy = y - pos[iy]
        dz = z - pos[iz]
        r2 = dx*dx + dy*dy + dz*dz
        if r2 < a2: return 1
    return 0

# INPUT
#  n   : number of particles
#  phi : volume fraction
#  a   : radius of particles
# OUTPUT
#  x[n*3] : configuration
def blob (n, phi, a, seed=0):
    pos = []
    # estimate the radius from phi
    # phi = (4/3) pi a^3 N / (4/3) pi R^3  = N (a/R)^3
    # => R/a = (phi / N)^(-1/3)  = (N / phi)^(1/3)
    R = a * math.pow(float(n) / phi, 1.0/3.0)
    R2 = R*R
    d2 = 4.0*a*a

    random.seed(seed)
    i = 0
    while 1:
        if i >= n: break

        itry = 0
        while 1:
            itry += 1
            if itry > 1000:
                #pos = pos.pop()
                #pos = pos.pop()
                #pos = pos.pop()
                # the above is not working...

                #new = []
                #for j in range(len(pos)-3):
                #    new.append(pos[j])
                #pos = new
                #i -= 1

                # withdraw everything and start again
                pos = []
                i = 0

                break
                
            x = 2.0 * R * (random.random() - 0.5)
            y = 2.0 * R * (random.random() - 0.5)
            z = 2.0 * R * (random.random() - 0.5)
            r2 = x*x + y*y + z*z
            if r2 > R2: continue

            # check overlap with other particles in x[]
            if check_overlap (pos, x, y, z, d2) != 0: continue # overlap

            # add (x,y,z) into pos[]
            pos.append(x)
            pos.append(y)
            pos.append(z)
            i += 1
            break
                
    return pos


# INPUT
#  plxyz : 'x', 'y', 'z' for the plane
#  n   : number of particles
#  phi : area fraction
#  a   : radius of particles
# OUTPUT
#  x[n*3] : configuration
def blob2 (plxyz, n, phi, a, seed=0):
    pos = []
    # estimate the radius from phi
    # phi = 4 pi a^2 N / 4 pi R^2  = N (a/R)^2
    # => R/a = (phi / N)^(-1/2)  = sqrt(N / phi)
    R = a * math.sqrt(float(n) / phi)
    R2 = R*R
    d2 = 4.0*a*a

    random.seed(seed)
    i = 0
    while 1:
        if i >= n: break

        itry = 0
        while 1:
            itry += 1
            if itry > 1000:
                #pos = pos.pop()
                #pos = pos.pop()
                #pos = pos.pop()
                # the above is not working...

                #new = []
                #for j in range(len(pos)-3):
                #    new.append(pos[j])
                #pos = new
                #i -= 1

                # withdraw everything and start again
                pos = []
                i = 0

                break
                
            x = 2.0 * R * (random.random() - 0.5)
            y = 2.0 * R * (random.random() - 0.5)
            z = 2.0 * R * (random.random() - 0.5)
            if plxyz == 'x':   x = 0.0
            elif plxyz == 'y': y = 0.0
            elif plxyz == 'z': z = 0.0
            else:
                print 'invalid plane index'
                sys.exit(1)

            r2 = x*x + y*y + z*z
            if r2 > R2: continue

            # check overlap with other particles in x[]
            if check_overlap (pos, x, y, z, d2) != 0: continue # overlap

            # add (x,y,z) into pos[]
            pos.append(x)
            pos.append(y)
            pos.append(z)
            i += 1
            break

    return pos


# INPUT
#  n   : number of particles
#  Rr  := R / r, where
#        R is the radius of the torous center
#        r is the radius of the torous circle
#  phi : area fraction
#  a   : radius of particles
# OUTPUT
#  x[n*3] : configuration
def torous (n, Rr, phi, a, seed=0):
    L2 = 1.0 / Rr / Rr; # = (r/R)^2
    R = math.pow(2.0 * float(n) * math.pow(a, 3.0)
                 / (3.0 * math.pi * phi * L2),
                 1.0 / 3.0)
    r = R / Rr;
    print '; R = %f, r = %f, lambda = 1/%f'%(R,r, R/r)
    pos = []

    d2 = 4.0*a*a

    random.seed(seed)
    i = 0
    while 1:
        if i >= n: break

        itry = 0
        while 1:
            itry += 1
            if itry > 1000:
                # withdraw everything and start again
                pos = []
                i = 0
                break
                
            x = 2.0 * (R+r) * (random.random() - 0.5)
            y = 2.0 * (R+r) * (random.random() - 0.5)
            z = 2.0 * r     * (random.random() - 0.5)

            rr = math.sqrt(x*x + y*y)
            if rr > (R+r) or rr < (R-r): continue
            rr2 = (rr-R)*(rr-R) + z*z
            if rr2 > r*r: continue

            # check overlap with other particles in x[]
            if check_overlap (pos, x, y, z, d2) != 0: continue # overlap

            # add (x,y,z) into pos[]
            pos.append(x)
            pos.append(y)
            pos.append(z)
            i += 1
            break

    return pos


def rotate(direction, degree, x):
    theta = degree*math.pi/180.0 # radian
    s = math.sin(theta)
    c = math.cos(theta)

    # make rotation matrix
    if direction == 'x':
        # x direction
        rot = [1.0, 0.0, 0.0,
               0.0, c,   -s, 
               0.0, s,    c]
    elif direction == 'y':
        # y direction
        rot = [c,   0.0, s,
               0.0, 1.0, 0.0,
               -s,  0.0, c]
    elif direction == 'z':
        # z direction
        rot = [c,   -s,  0.0,
               s,    c,  0.0,
               0.0, 0.0, 1.0]
    else:
        print 'invalid direction in rotate()'
        sys.exit(1)

    # apply the rotation
    n = len(x)/3
    z = []
    for i in range(n):
        ix = i*3
        iy = ix+1
        iz = ix+2
        z.append(  rot [0*3+0]*x[ix]\
                 + rot [0*3+1]*x[iy]\
                 + rot [0*3+2]*x[iz])
        z.append(  rot [1*3+0]*x[ix]\
                 + rot [1*3+1]*x[iy]\
                 + rot [1*3+2]*x[iz])
        z.append(  rot [2*3+0]*x[ix]\
                 + rot [2*3+1]*x[iy]\
                 + rot [2*3+2]*x[iz])

    return z


def translate(direction, dx, x):
    z = []
    if direction == 'x':
        # x direction
        z.append(x[0] + dx)
        z.append(x[1])
        z.append(x[2])
    elif direction == 'y':
        # y direction
        z.append(x[0])
        z.append(x[1] + dx)
        z.append(x[2])
    elif direction == 'z':
        # z direction
        z.append(x[0])
        z.append(x[1])
        z.append(x[2] + dx)

    return z



def print_config (x):
    # print arguments as a comment
    print '; configuration generated by gen-config.py from...'
    print ';  ',
    for i in range(len(sys.argv)):
        print '%s'%(sys.argv[i]),
    print ''
    print '(define x #('
    for i in range(len(x)/3):
        print '  %f %f %f ; %d'%(x[i*3],x[i*3+1],x[i*3+2],i)
    print '))'


def usage():
    print '$Id: gen-config.py,v 1.2 2007/11/18 05:20:20 kichiki Exp $'
    print 'USAGE:'
    print '\t-n : number of particles'
    print 'CONFIGURATION OPTIONS'
    print '\t-line d : make a line with d spacing in x direction'
    print '\t-arc d R t0 : make an arc\n'\
          '\t\t R is the radius of the arc\n'\
          '\t\t t0 is the init arg\n'\
          '\t\t d is the spacing in xz plane'
    print '\t-blob  phi a seed : make a spherical blob.\n'\
    '\t\t phi is the volume fraction\n'\
    '\t\t a is the radius of particles\n'\
    '\t\t seed is for the random number generator'
    print '\t-blob2 [xyz] phi a seed : make a circular blob in y=0 plane.\n'\
    '\t\t [xyz] : the plane on which particles are\n'\
    '\t\t phi is the volume fraction\n'\
    '\t\t a is the radius of particles\n'\
    '\t\t seed is for the random number generator'
    print '\t-torous R/r phi a seed : make a horizontal torous blob.\n'\
    '\t\t R is the radius of the torous center\n'\
    '\t\t r is the radius of the torous circle\n'\
    '\t\t phi is the volume fraction\n'\
    '\t\t a is the radius of particles\n'\
    '\t\t seed is for the random number generator'
    print 'OPERATION OPTIONS (you can give them as much as you want)'
    print '\t-rotate [xyz] D : rotate D degree in [xyz] direction'
    print '\t-translate [xyz] D : traslate D degree in [xyz] direction'
    sys.exit ()


def main():
    n = 0
    flag_conf = 0 # 0 == line
    r = 0.0
    op_type  = []
    op_dir   = []
    op_param = []
    i = 1
    while i < len(sys.argv):
        if sys.argv[i] == '-n':
            n = int(sys.argv[i+1])
            i += 2
        elif sys.argv[i] == '-line':
            flag_conf = 0 # line
            d = float(sys.argv[i+1])
            i += 2
        elif sys.argv[i] == '-arc':
            flag_conf = 1 # arc
            d = float(sys.argv[i+1])
            r  = float(sys.argv[i+2])
            t0 = float(sys.argv[i+3])
            i += 4
            # convert degree to radian
            t0 *= math.pi/180.0
        elif sys.argv[i] == '-blob':
            flag_conf = 2 # blob
            phi  = float(sys.argv[i+1])
            a    = float(sys.argv[i+2])
            seed = int(sys.argv[i+3])
            i += 4
        elif sys.argv[i] == '-blob2':
            flag_conf = 3 # blob2
            plxyz = sys.argv[i+1]
            phi   = float(sys.argv[i+2])
            a     = float(sys.argv[i+3])
            seed  = int(sys.argv[i+4])
            i += 5
        elif sys.argv[i] == '-torous':
            flag_conf = 4 # torous
            Rr    = float(sys.argv[i+1])
            phi   = float(sys.argv[i+2])
            a     = float(sys.argv[i+3])
            seed  = int(sys.argv[i+4])
            i += 5
        elif sys.argv[i] == '-rotate':
            op_type.append(0) # rotate
            op_dir.append(sys.argv[i+1])
            op_param.append(float(sys.argv[i+2]))
            i += 3
        elif sys.argv[i] == '-translate':
            op_type.append(1) # translate
            op_dir.append(sys.argv[i+1])
            op_param.append(float(sys.argv[i+2]))
            i += 3
        else:
            usage()

    # generate initial configuration
    if flag_conf == 0:
        # line
        x = []
        for i in range(n):
            x.append(float(i-(n/2))*d) # x component
            x.append(0.0)              # y component
            x.append(0.0)              # z component
    elif flag_conf == 1:
        # arc
        print 'R = %f, t0 = %f, d = %f'%(r, t0, d)
        x = []
        t = t0
        dt = 2.0 * math.asin (d/2.0/r)
        print 'dt = %f [radian] = %f [degree]'%(dt, dt*180.0/math.pi)
        for i in range(n):
            x.append(r*math.cos(t)) # x component
            x.append(0.0)           # y component
            x.append(r*math.sin(t)) # z component
            t += dt
    elif flag_conf == 2:
        # blob
        x = blob (n, phi, a)
    elif flag_conf == 3:
        # blob2
        x = blob2 (plxyz, n, phi, a)
    elif flag_conf == 4:
        # torous
        x = torous (n, Rr, phi, a)
    else:
        print 'invalid configuration frag'
        sys.exit(1)

    # rotation loop
    for i in range(len(op_type)):
        if op_type[i] == 0:
            # rotate
            x = rotate(op_dir[i], op_param[i], x)
        elif op_type == 1:
            # translate
            x = translate(op_dir[i], op_param[i], x)
        else:
            print 'invalid operation type'
            sys.exit(1)

    # print the configuration
    print_config (x)


if __name__ == "__main__":
    main()
