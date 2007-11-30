# visualization program for stokes-nc file by VTK
# Copyright (C) 2006-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
# $Id: stvis-vtk.py,v 1.13 2007/11/30 06:34:04 kichiki Exp $
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

#sys.path.append('/somewhere/vtk-python-binding/lib64/python2.4/site-packages')
import vtk
#from vtk.util.colors import peacock, tomato
from vtk.util.colors import *

#sys.path.append('/somewhere/ryuon/stokes/python')
import stokes

import ryuon_vtk


def usage():
    print '$Id: stvis-vtk.py,v 1.13 2007/11/30 06:34:04 kichiki Exp $'
    print 'USAGE:'
    print '\t-f or --file : stokes-nc-file'
    print '\t-step n      : draw every n steps (default: 1, all frames)'
    print '\t-b or --bond : draw bond connecting particles'
    print '\t-s           : print step instead of time'
    print '\t-ts          : print time with step'
    print '\t-bottom      : follow the camera at the bottom of the blob'
    print '\t             : (default) follow the center of mass'
    print '\t             : this option is valid only for non-periodic system'
    print '\t-p or --pic  : generate PNG images'
    sys.exit ()


def main():
    filename = ''
    nstep = 1
    flag_bond = 0
    flag_step = 0
    flag_bottom = 0
    flag_pic = 0
    i = 1
    while i < len(sys.argv):
        if sys.argv[i] == '-f' or sys.argv[i] == '--file':
            filename = sys.argv[i+1]
            i += 2
        elif sys.argv[i] == '-step':
            nstep = int(sys.argv[i+1])
            i += 2
        elif sys.argv[i] == '-b' or sys.argv[i] == '--bond':
            flag_bond = 1
            i += 1
        elif sys.argv[i] == '-s':
            flag_step = 1
            i += 1
        elif sys.argv[i] == '-ts':
            flag_step = 2
            i += 1
        elif sys.argv[i] == '-bottom':
            flag_bottom = 1
            i += 1
        elif sys.argv[i] == '-p' or sys.argv[i] == '--pic':
            flag_pic = 1
            i += 1
        else:
            usage()
    if filename == '' : usage()

    nc = stokes.stokes_nc_open (filename)
    #stokes.stokes_nc_print_actives(nc, stokes.get_stdout())
    print 'ntime = ', nc.ntime
    lattice = stokes.darray(3)
    stokes.stokes_nc_get_array1d (nc, 'l', lattice)
    print 'lattice = ', lattice[0], lattice[1], lattice[2]

    # x[] : center of particles
    x  = stokes.darray(nc.np  * nc.nvec)
    # q[] : quaternion
    if nc.flag_q != 0:
        q  = stokes.darray(nc.np  * nc.nquat)
    else:
        q = []
    # a[] : radius of mobile particles
    if nc.flag_a != 0:
        a = stokes.darray(nc.np)
        stokes.stokes_nc_get_array1d (nc, "a", a)
    else:
        a = []
    # af[] : radius of fixed particles
    if nc.flag_af != 0:
        af = stokes.darray(nc.npf)
        stokes.stokes_nc_get_array1d (nc, "af", af)
    else:
        af = []
    
    # fixed particles
    if nc.npf > 0:
        xf0  = stokes.darray(nc.npf * nc.nvec)
        stokes.stokes_nc_get_data0 (nc, "xf0", xf0)
        (pfActor,pfGlyph) = ryuon_vtk.make_pActor ((1,0,0), 1, 20)
        pfData = ryuon_vtk.make_pData(nc.npf, xf0, af)
        pfGlyph.SetInput(pfData)

    # then, make Actor
    if nc.flag_q != 0:
        # with quaternion
        pActors = []
        for i in range (nc.np):
            pActors.append(ryuon_vtk.make_pActor_with_quaternion())
    else:
        # no quaternion in the result
        (pActor,pGlyph) = ryuon_vtk.make_pActor((1,1,1), 1, 8)
        # pGlyph.SetInput(pData)

    # bond Actor
    if flag_bond != 0:
        (bondActor,bond) = ryuon_vtk.make_bondActor (khaki)


    ## prepare renderer
    ren = vtk.vtkRenderer()
    win = vtk.vtkRenderWindow()
    win.AddRenderer(ren)

    ren.SetBackground (0.1, 0.2, 0.4)
    win.SetSize(480,360)
    # iren = vtk.vtkRenderWindowInteractor ()
    # iren.SetRenderWindow(win)

    if nc.flag_q != 0:
        for i in range(nc.np):
            ren.AddActor(pActors[i])
    else:
        ren.AddActor(pActor)
    if nc.npf > 0:
        ren.AddActor(pfActor)
    if flag_bond != 0:
        ren.AddActor(bondActor)

    # periodic box
    if lattice[0] != 0.0 or lattice[1] != 0.0 or lattice[2] != 0.0:
        CubeEdges = ryuon_vtk.make_cubeActor (lattice)
        ren.AddActor(CubeEdges)

    # axes
    #axesActor = ryuon_vtk.make_axesActor()
    #ren.AddActor(axesActor)

    # text
    textActor = ryuon_vtk.make_textActor()
    ren.AddActor(textActor)


    aCamera = vtk.vtkCamera()
    ren.SetActiveCamera (aCamera)

    if lattice[0] != 0.0 or lattice[1] != 0.0 or lattice[2] != 0.0:
        # periodic boundary
        if lattice[0] > lattice[2]:
            l = lattice[0]
        else:                  
            l = lattice[2]

        aCamera.SetFocalPoint (0.5*lattice[0], 0,      0.4*lattice[2])
        aCamera.SetPosition   (0.5*lattice[0], -1.6*l, 0.5*lattice[2])

    # loop 
    while 1:
        i = 0
        while i < nc.ntime:
            t = stokes.stokes_nc_get_time_step (nc, i)
            stokes.stokes_nc_get_data (nc, "x", i, x)
            if nc.flag_q != 0:
                # with quaternion
                for j in range (nc.np):
                    stokes.stokes_nc_get_data (nc, "q", i, q)
                    m = ryuon_vtk.Q2M (q[j*4+0],q[j*4+1],q[j*4+2],q[j*4+3])
                    trans = vtk.vtkTransform()
                    trans.SetMatrix((m[0],m[3],m[6], x[j*3+0],
                                     m[1],m[4],m[7], x[j*3+1],
                                     m[2],m[5],m[8], x[j*3+2],
                                     0.0, 0.0, 0.0,  1.0));
                    pActors[j].SetUserTransform(trans)
            else:
                # no quaternion in the result
                pData = ryuon_vtk.make_pData(nc.np, x, a)
                pGlyph.SetInput(pData)
    
            # update bond Actor
            if flag_bond != 0:
                bData = ryuon_vtk.make_bData(nc.np, x)
                bond.SetInput(bData)

            if lattice[0] == 0.0 and lattice[1] == 0.0 and lattice[2] == 0.0:
                # non-periodic boundary
                #(cx,cy,cz, lx,ly,lz) = bounding_box (nc.np, x)
                (cx,cy,cz, lx0,ly0,lz0, lx1, ly1, lz1)\
                           = ryuon_vtk.bounding_box_ (nc.np, x)
                lx = lx1 - lx0
                ly = ly1 - ly0
                lz = lz1 - lz0
                # centered by COM
                if lx > lz:
                    l = lx
                else:                  
                    l = lz
                # prevent to go far away
                if l > 50: l = 50
    
                if flag_bottom == 0:
                    aCamera.SetFocalPoint (cx, cy,       cz)
                    aCamera.SetPosition   (cx, cy-3.0*l, cz+0.1*l)
                else:
                    # bottom align
                    aCamera.SetFocalPoint (cx, cy,       lz0)
                    aCamera.SetPosition   (cx, cy-3.0*l, lz0+0.1*l)

            if flag_step == 0:
                textActor.SetInput ('time: %f'%(t))
            elif flag_step == 1:
                textActor.SetInput ('step: %d'%(i))
            else:
                textActor.SetInput ('time: %f (%d step)'%(t,i))
            win.Render()

            if flag_pic == 1:
                w2if = vtk.vtkWindowToImageFilter()
                w2if.SetInput(win)
                wr = vtk.vtkPNGWriter()
                wr.SetInput(w2if.GetOutput())
                wr.SetFileName('%s-%05d.png'%(filename,i))
                wr.Write()

            i += nstep
        # end of while i < nc.ntime:
    # end of while 1:


if __name__ == "__main__":
    main()
