# visualization program for stokes-nc file by VTK
# Copyright (C) 2006-2008 Kengo Ichiki <kichiki@users.sourceforge.net>
# $Id: stvis-vtk.py,v 1.15 2008/06/03 03:00:39 kichiki Exp $
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


# INPUT
#  r     : radius of the cylinder
#  x0[3] : one end of the cylinder
#  x1[3] : the other end of the cylinder
#  rgb : color
#  opacity :
def make_cylinderActor (r, x0, x1, rgb, opacity):
    points = vtk.vtkPoints()
    lines  = vtk.vtkCellArray()
    lines.InsertNextCell(2)
    # point 0
    points.InsertNextPoint(x0[0], x0[1], x0[2])
    lines.InsertCellPoint(0)
    # point 1
    points.InsertNextPoint(x1[0], x1[1], x1[2])
    lines.InsertCellPoint(1)

    cData = vtk.vtkPolyData()
    cData.SetPoints(points)
    cData.SetLines(lines)

    c = vtk.vtkTubeFilter()
    c.SetNumberOfSides(8)
    c.SetInput(cData)
    c.SetRadius(r)

    cMapper = vtk.vtkPolyDataMapper()
    cMapper.SetInput(c.GetOutput())

    cActor = vtk.vtkActor()
    cActor.SetMapper(cMapper)
    cActor.GetProperty().SetColor(rgb[0], rgb[1], rgb[2])
    cActor.GetProperty().SetOpacity(opacity)

    return cActor

# INPUT
#  x[3]  : center of the sphere
#  r     : radius of the sphere
#  rgb : color
#  opacity :
def make_sphereActor (x, r, rgb, opacity):
    points = vtk.vtkPoints()
    points.InsertNextPoint(x[0], x[1], x[2])

    diameter = vtk.vtkDoubleArray()
    diameter.SetNumberOfComponents(1)
    diameter.InsertNextTuple1(2.0*r)

    pData = vtk.vtkPolyData()
    pData.SetPoints(points)
    pData.GetPointData().SetScalars(diameter)

    pSource = vtk.vtkSphereSource()
    pSource.SetPhiResolution(16)
    pSource.SetThetaResolution(16)

    pGlyph = vtk.vtkGlyph3D()
    pGlyph.SetSource(pSource.GetOutput())
    pGlyph.SetInput(pData)
    pGlyph.ScalingOn()
    pGlyph.SetScaleModeToScaleByScalar()

    pMapper = vtk.vtkPolyDataMapper()
    pMapper.ScalarVisibilityOff()
    pMapper.SetInput(pGlyph.GetOutput())

    pActor = vtk.vtkActor()
    pActor.SetMapper(pMapper)
    pActor.GetProperty().SetColor(rgb[0], rgb[1], rgb[2])
    pActor.GetProperty().SetOpacity(opacity)

    return pActor


def usage():
    print '$Id: stvis-vtk.py,v 1.15 2008/06/03 03:00:39 kichiki Exp $'
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
    print '\t-c or --cylinder : give radius to show cylinder'
    print '\t--sphere         : give radius to show radius'
    print '\t-d or --dumbbell : give 4 parameters (R1, R2, L, r), where'
    print '\t                   R1 : left cavity radius'
    print '\t                   R2 : right cavity radius'
    print '\t                   L  : length of the cylinder'
    print '\t                   r  : radius of the cylinder'
    print '\t--hex2d          : give 3 parameters (R, r, L), where'
    print '\t                   R : cavity radius'
    print '\t                   r : cylinder radius'
    print '\t                   L : lattice spacing'
    sys.exit ()


def main():
    filename = ''
    nstep = 1
    flag_bond = 0
    flag_step = 0
    flag_bottom = 0
    flag_pic = 0
    flag_cylinder = 0
    flag_sphere = 0
    flag_dumbbell = 0
    flag_hex2d = 0
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
        elif sys.argv[i] == '-c' or sys.argv[i] == '--cylinder':
            flag_cylinder = 1
            cylinder_radius = float(sys.argv[i+1])
            i += 2
        elif sys.argv[i] == '--sphere':
            flag_sphere = 1
            cavity_radius = float(sys.argv[i+1])
            i += 2
        elif sys.argv[i] == '-d' or sys.argv[i] == '--dumbbell':
            flag_dumbbell = 1
            cavity1_radius = float(sys.argv[i+1])
            cavity2_radius = float(sys.argv[i+2])
            cylinder_length = float(sys.argv[i+3])
            cylinder_radius = float(sys.argv[i+4])
            i += 5
        elif sys.argv[i] == '--hex2d':
            flag_hex2d = 1
            cavity_radius = float(sys.argv[i+1])
            cylinder_radius = float(sys.argv[i+2])
            lattice_spacing = float(sys.argv[i+3])
            i += 4
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

        if lattice[0] == 0.0 and lattice[1] == 0.0 and lattice[2] == 0.0:
            # non-periodic boundary
            pfData = ryuon_vtk.make_pData(nc.npf, xf0, af)
            (pfActor,pfGlyph) = ryuon_vtk.make_pActor ((1,0,0), 1, 20)
        else:
            # periodic boundary
            pfData = ryuon_vtk.make_pData_periodic(nc.npf, xf0, af,
                                               lattice)
            (pfActor,pfGlyph) = ryuon_vtk.make_pActor_periodic ((1,0,0),
                                                                1, 20,
                                                                lattice)
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

    # cylinder
    if flag_cylinder != 0:
        cActor = make_cylinderActor (r=cylinder_radius,
                                     x0=(0.0, 0.0, 0.0),
                                     x1=(200.0, 0.0, 0.0),
                                     rgb=(0.0, 0.5, 0.0),
                                     opacity=0.3)
        ren.AddActor(cActor)

    # sphere
    if flag_sphere != 0:
        cavActor = make_sphereActor (x=(0.0, 0.0, 0.0),
                                      r=cavity_radius,
                                      rgb=(0.0, 0.5, 0.0),
                                      opacity=0.3)
        ren.AddActor(cavActor)

    # dumbbell
    if flag_dumbbell != 0:
        cActor = make_cylinderActor (r=cylinder_radius,
                                     x0=(-0.5 * cylinder_length, 0.0, 0.0),
                                     x1=(+0.5 * cylinder_length, 0.0, 0.0),
                                     rgb=(0.0, 0.5, 0.0),
                                     opacity=0.3)
        ren.AddActor(cActor)

        theta1 = math.asin(cylinder_radius / cavity1_radius)
        cav1Actor = make_sphereActor (x=(-0.5 * cylinder_length\
                                         -cavity1_radius * math.cos(theta1),
                                         0.0, 0.0),
                                      r=cavity1_radius,
                                      rgb=(0.0, 0.5, 0.0),
                                      opacity=0.3)
        ren.AddActor(cav1Actor)
        
        theta2 = math.asin(cylinder_radius / cavity2_radius)
        cav2Actor = make_sphereActor (x=(+0.5 * cylinder_length\
                                         +cavity2_radius * math.cos(theta2),
                                         0.0, 0.0),
                                      r=cavity2_radius,
                                      rgb=(0.0, 0.5, 0.0),
                                      opacity=0.3)
        ren.AddActor(cav2Actor)
        
    # hex2d
    if flag_hex2d != 0:
        cav1Actor = make_sphereActor (x=(0.0, 0.0, 0.0),
                                      r=cavity_radius,
                                      rgb=(0.0, 0.5, 0.0),
                                      opacity=0.3)
        ren.AddActor(cav1Actor)
        
        cav2Actor = make_sphereActor (x=(lattice_spacing, 0.0, 0.0),
                                      r=cavity_radius,
                                      rgb=(0.0, 0.5, 0.0),
                                      opacity=0.3)
        ren.AddActor(cav2Actor)
        
        cav3Actor = make_sphereActor (x=(0.5*lattice_spacing, math.sqrt(3.0) * 0.5 *lattice_spacing, 0.0),
                                      r=cavity_radius,
                                      rgb=(0.0, 0.5, 0.0),
                                      opacity=0.3)
        ren.AddActor(cav3Actor)
        
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
                #if l > 50: l = 50
                if l > 100: l = 100
    
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
