# dump stokes-netcdf
# Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
# $Id: stncdump.py,v 1.2 2007/04/20 02:23:12 kichiki Exp $
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


def usage():
    print '$Id: stncdump.py,v 1.2 2007/04/20 02:23:12 kichiki Exp $'
    print 'USAGE:'
    print '\t-f or --file : stokes-nc-file'
    print '\t-line        : all particles are in a single line for each time'\
          '(default: one particle per line)'
    sys.exit ()


def main():
    filename = ''
    flag_line = 0
    i = 1
    while i < len(sys.argv):
        if sys.argv[i] == '-f' or sys.argv[i] == '--file':
            filename = sys.argv[i+1]
            i += 2
        elif sys.argv[i] == '-line':
            flag_line = 1
            i += 1
        else:
            usage()
    if filename == '': usage()

    nc = stokes.stokes_nc_open (filename)
    stokes.stokes_nc_print_actives(nc, stokes.get_stdout())
    print ''
    
    lattice = stokes.darray(3)
    stokes.stokes_nc_get_l (nc, lattice)
    print 'lattice %f %f %f\n'%(lattice[0], lattice[3], lattice[2])

    pos = stokes.darray(nc.np  * nc.nvec)

    # print the infor about the output format
    if flag_line == 0:
        print '# t, i, x, y, z'
    else:
        print '# t, x, y, z (for particle 0), x, y, z (for 1), ... upto particle %d'%(nc.np)

    i = 0
    while 1:
        t = stokes.stokes_nc_get_time_step (nc, i)

        stokes.stokes_nc_get_data (nc, "x", i, pos)
        if flag_line == 0:
            for j in range(nc.np):
                x = pos[j*3]
                y = pos[j*3+1]
                z = pos[j*3+2]
                print '%f %d %f %f %f'%(t, j, x, y, z)
        else:
            print t,
            for j in range(nc.np):
                x = pos[j*3]
                y = pos[j*3+1]
                z = pos[j*3+2]
                print x, y, z,
            print ''

        i += 1

if __name__ == "__main__":
    main()