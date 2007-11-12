# dump stokes-netcdf
# Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
# $Id: stncdump.py,v 1.7 2007/11/12 03:34:53 kichiki Exp $
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
    print '$Id: stncdump.py,v 1.7 2007/11/12 03:34:53 kichiki Exp $'
    print 'USAGE:'
    print '\t-f or --file : stokes-nc-file'
    print '\t-line        : all particles are in a single line for each time\n'\
          '\t\t (default: one particle per line)'
    print '\t-step n      : write the config at step n\n'\
          '\t\t step 1000 is the last step for 1000 run.\n'\
          '\t\t (step 0, the initial config, does not exist in the nc file.)\n'
    print '\t-com0        : set COM to 0 with -step option.\n'
    sys.exit ()



def main():
    filename = ''
    flag_line = 0
    step = -1
    flag_com0 = 0
    i = 1
    while i < len(sys.argv):
        if sys.argv[i] == '-f' or sys.argv[i] == '--file':
            filename = sys.argv[i+1]
            i += 2
        elif sys.argv[i] == '-line':
            flag_line = 1
            i += 1
        elif sys.argv[i] == '-step':
            step = int(sys.argv[i+1])
            step -= 1
            i += 2
        elif sys.argv[i] == '-com0':
            flag_com0 = 1
            i += 1
        else:
            usage()
    if filename == '': usage()

    nc = stokes.stokes_nc_open (filename)

    # pos[] : center of particles
    pos = stokes.darray(nc.np  * nc.nvec)

    # q[] : quaternion
    if nc.flag_q != 0:
        q  = stokes.darray(nc.np  * nc.nquat)
    else:
        q = []

    # lattice
    lattice = stokes.darray(3)
    stokes.stokes_nc_get_array1d (nc, 'l', lattice)

    # extract the config at the step
    if step >= 0:
        if step > nc.ntime:
            print 'out of the range %d <= %d'%(step, nc.ntime)
            sys.exit(1)

        # read the config at the step
        stokes.stokes_nc_get_data (nc, "x", step, pos)
        if nc.flag_q != 0:
            stokes.stokes_nc_get_data (nc, "q", step, q)

        comx = 0.0
        comy = 0.0
        comz = 0.0
        if flag_com0 != 0:
            for i in range(nc.np):
                comx += pos[i*3]
                comy += pos[i*3+1]
                comz += pos[i*3+2]
            comx /= float(nc.np)
            comy /= float(nc.np)
            comz /= float(nc.np)

        # print the config
        # print arguments as a comment
        print '; config at %d step of %s, generated by'%(step+1, filename)
        print ';  ',
        for i in range(len(sys.argv)):
            print '%s'%(sys.argv[i]),
        print ''
        print '(define x #('
        for i in range(nc.np):
            print '  %f %f %f ; %d'%(pos[i*3]   - comx,
                                     pos[i*3+1] - comy,
                                     pos[i*3+2] - comz,
                                     i)
        print '))'

        if nc.flag_q != 0:
            stokes.stokes_nc_get_data (nc, "q", step, q)
            print '(define q #('
            for i in range(nc.np):
                print '  %f %f %f %f ; %d'\
                      %(q[i*4],q[i*4+1],q[i*4+2],q[i*4+3],i)
            print '))'
        sys.exit(1)
        # done!

    if flag_line == 0:
        # print some general informations
        stokes.stokes_nc_print_actives(nc, stokes.get_stdout())
        print ''
        print 'ntime = %d'%(nc.ntime)
        print ''

        # imposed flows
        if nc.flag_ui0 == 1:
            ui0 = stokes.darray(3)
            stokes.stokes_nc_get_array1d (nc, "Ui0", ui0)
            print 'ui0 = %f %f %f'%(ui0[0], ui0[1], ui0[2])
        if nc.flag_oi0 == 1:
            oi0 = stokes.darray(3)
            stokes.stokes_nc_get_array1d (nc, "Oi0", oi0)
            print 'oi0 = %f %f %f'%(oi0[0], oi0[1], oi0[2])
        if nc.flag_ei0 == 1:
            ei0 = stokes.darray(5)
            stokes.stokes_nc_get_array1d (nc, "Ei0", ei0)
            print 'ei0 = %f %f %f %f %f'\
                  %(ei0[0], ei0[1], ei0[2], ei0[3], ei0[4])
        print ''

        print 'lattice %f %f %f\n'%(lattice[0], lattice[3], lattice[2])

    # print the infos about the output format
    if flag_line == 0:
        if nc.flag_q != 0:
            print '# t, i, x, y, z, q1, q2, q3, q4'
        else:
            print '# t, i, x, y, z'
    else:
        if nc.flag_q != 0:
            print '# t, x, y, z, q1, q2, q3, q4 (for particle 0),'\
                  ' ... upto particle %d'%(nc.np)
        else:
            print '# t, x, y, z (for particle 0), x, y, z (for 1),'\
                  ' ... upto particle %d'%(nc.np)

    for i in range(nc.ntime):
        t = stokes.stokes_nc_get_time_step (nc, i)

        stokes.stokes_nc_get_data (nc, "x", i, pos)
        if nc.flag_q != 0:
            stokes.stokes_nc_get_data (nc, "q", i, q)

        if flag_line == 0:
            for j in range(nc.np):
                x = pos[j*3]
                y = pos[j*3+1]
                z = pos[j*3+2]
                if nc.flag_q != 0:
                    print '%f %d %f %f %f %f %f %f %f'\
                          %(t, j, x, y, z, q[0], q[1], q[2], q[3])
                else:
                    print '%f %d %f %f %f'%(t, j, x, y, z)
        else:
            print t,
            for j in range(nc.np):
                x = pos[j*3]
                y = pos[j*3+1]
                z = pos[j*3+2]
                print x, y, z,
                if nc.flag_q != 0:
                    print q[0], q[1], q[2], q[3],
            print ''


if __name__ == "__main__":
    main()
