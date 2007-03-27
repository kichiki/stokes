# dump stokes-netcdf
# Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
# $Id: stncdump.py,v 1.1 2007/03/27 07:03:26 kichiki Exp $
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
import sys
#sys.path.append('/somewhere/ryuon/stokes/python')
import stokes


def usage():
    print '$Id: stncdump.py,v 1.1 2007/03/27 07:03:26 kichiki Exp $'
    print 'USAGE:'
    print '\t-f or --file : stokes-nc-file'
    sys.exit ()


def main():
    filename = ''
    i = 1
    while i < len(sys.argv):
        if sys.argv[i] == '-f' or sys.argv[i] == '--file':
            filename = sys.argv[i+1]
            i += 2
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

    i = 0
    while 1:
        t = stokes.stokes_nc_get_time_step (nc, i)

        stokes.stokes_nc_get_data (nc, "x", i, pos)
        for j in range(nc.np):
            x = pos[j*3]
            y = pos[j*3+1]
            z = pos[j*3+2]
            print '%f %d %f %f %f'%(t, j, x, y, z)

        i += 1

if __name__ == "__main__":
    main()
