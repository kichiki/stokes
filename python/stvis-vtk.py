# visualization program for stokes-nc file by VTK
# Copyright (C) 2006-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
# $Id: stvis-vtk.py,v 1.7 2007/05/13 23:13:47 kichiki Exp $
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
import vtk
from vtk.util.colors import peacock, tomato
import sys
import math
#sys.path.append('/somewhere/ryuon/stokes/python')
import stokes


# make an Actor for mobile particles at (0,0,0) with radius 1
# this is for the result with quaternion
def make_pActor ():
    pSource = vtk.vtkSphereSource()
    pSource.SetPhiResolution(20)
    pSource.SetThetaResolution(20)

    pGlyph = vtk.vtkGlyph3D()
    pGlyph.SetSource(pSource.GetOutput())

    #pGlyph.ScalingOn()
    #pGlyph.SetScaleModeToScaleByScalar()

    pos = vtk.vtkPoints()
    pos.InsertNextPoint(0, 0, 0)

    pData = vtk.vtkPolyData()
    pData.SetPoints(pos)

    diameter = vtk.vtkDoubleArray()
    diameter.SetNumberOfComponents(1)
    diameter.InsertNextTuple1(2.0)

    pData.GetPointData().SetScalars(diameter)
    pGlyph.SetInput(pData)

    pNormals = vtk.vtkPolyDataNormals()
    pNormals.SetInput(pGlyph.GetOutput())

    ###########################################################################
    # clip plane
    plane = vtk.vtkPlane()
    plane.SetOrigin(0, 0, 0)
    plane.SetNormal(0, 1, 0) # upper side

    clipper = vtk.vtkClipPolyData()
    clipper.SetInput(pNormals.GetOutput())
    clipper.SetClipFunction(plane)
    clipper.GenerateClipScalarsOn()
    clipper.GenerateClippedOutputOn()
    clipper.SetValue(0)

    clipMapper = vtk.vtkPolyDataMapper()
    clipMapper.SetInput(clipper.GetOutput())
    clipMapper.ScalarVisibilityOff()

    clipActor = vtk.vtkActor()
    clipActor.SetMapper(clipMapper)
    clipActor.GetProperty().SetColor(peacock)

    ###########################################################################
    # clipped part
    restMapper = vtk.vtkPolyDataMapper()
    restMapper.SetInput(clipper.GetClippedOutput())
    restMapper.ScalarVisibilityOff()

    restActor = vtk.vtkActor()
    restActor.SetMapper(restMapper)
    restActor.GetProperty().SetColor(tomato)

    ###########################################################################
    # assembly
    assembly = vtk.vtkAssembly()
    assembly.AddPart(clipActor)
    assembly.AddPart(restActor)
    assembly.RotateX(90)

    return assembly


# make transform matrix (3x3) by quaternion
def Q2M (q1,q2,q3,q4):
    m = []

    m.append(2.0*(q1*q1 + q4*q4 - .5))
    m.append(2.0*(q1*q2 + q3*q4))
    m.append(2.0*(q1*q3 - q2*q4))

    m.append(2.0*(q1*q2 - q3*q4))
    m.append(2.0*(q2*q2 + q4*q4 - .5))
    m.append(2.0*(q2*q3 + q1*q4))

    m.append(2.0*(q1*q3 + q2*q4))
    m.append(2.0*(q2*q3 - q1*q4))
    m.append(2.0*(q3*q3 + q4*q4 - .5))

    return m


# make pData for np particles
def make_pData (np, x, a):
    pos = vtk.vtkPoints()
    diameter = vtk.vtkDoubleArray()
    diameter.SetNumberOfComponents(1)
    for i in range(np):
        pos.InsertNextPoint(x[i*3], x[i*3+1], x[i*3+2])
        if a != []:
            diameter.InsertNextTuple1(2.0*a[i])
        else:
            diameter.InsertNextTuple1(2.0)

    # first make pData containing particle coordinates
    pData = vtk.vtkPolyData()
    pData.SetPoints(pos)
    pData.GetPointData().SetScalars(diameter)

    return pData


def usage():
    print '$Id: stvis-vtk.py,v 1.7 2007/05/13 23:13:47 kichiki Exp $'
    print 'USAGE:'
    print '\t-f or --file : stokes-nc-file'
    sys.exit ()


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
    #stokes.stokes_nc_print_actives(nc, stokes.get_stdout())
    print 'ntime = ', nc.ntime
    lattice = stokes.darray(3)
    stokes.stokes_nc_get_l (nc, lattice)
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
        stokes.stokes_nc_get_data0 (nc, "a", a)
    else:
        a = []
    # af[] : radius of fixed particles
    if nc.flag_af != 0:
        af = stokes.darray(nc.npf)
        stokes.stokes_nc_get_data0 (nc, "af", af)
    else:
        af = []
    

    # fixed particles
    if nc.npf > 0:
        xf0  = stokes.darray(nc.npf * nc.nvec)
        stokes.stokes_nc_get_data0 (nc, "xf0", xf0)
        pfData = make_pData(nc.npf, xf0, af)

        pfSource = vtk.vtkSphereSource()
        pfSource.SetPhiResolution(20)
        pfSource.SetThetaResolution(20)

        pfGlyph = vtk.vtkGlyph3D()
        pfGlyph.SetSource(pfSource.GetOutput())
        pfGlyph.SetInput(pfData)
        pfGlyph.ScalingOn()
        pfGlyph.SetScaleModeToScaleByScalar()

        pfMapper = vtk.vtkPolyDataMapper()
        pfMapper.ScalarVisibilityOff()
        pfMapper.SetInput(pfGlyph.GetOutput())

        pfActor = vtk.vtkActor()
        pfActor.SetMapper(pfMapper)
        # put Red color
        pfActor.GetProperty().SetColor(1,0,0)
        pfActor.GetProperty().SetOpacity(1)

    
    # then, make Actor
    if nc.flag_q != 0:
        # with quaternion
        pActors = []
        for i in range (nc.np):
            pActors.append(make_pActor())
    else:
        # no quaternion in the result
        pSource = vtk.vtkSphereSource()
        #pSource.SetPhiResolution(20)
        #pSource.SetThetaResolution(20)

        pGlyph = vtk.vtkGlyph3D()
        pGlyph.SetSource(pSource.GetOutput())
        # pGlyph.SetInput(pData)
        pGlyph.ScalingOn()
        pGlyph.SetScaleModeToScaleByScalar()

        pMapper = vtk.vtkPolyDataMapper()
        pMapper.ScalarVisibilityOff()
        pMapper.SetInput(pGlyph.GetOutput())

        pActor = vtk.vtkActor()
        pActor.SetMapper(pMapper)

        pActor.GetProperty().SetColor(1,1,1)
        pActor.GetProperty().SetOpacity(1)


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
    if nc.npf > 0: ren.AddActor(pfActor)

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
        for i in range(nc.ntime):
            stokes.stokes_nc_get_data (nc, "x", i, x)
            if nc.flag_q != 0:
                # with quaternion
                for j in range (nc.np):
                    stokes.stokes_nc_get_data (nc, "q", i, q)
                    m = Q2M (q[j*4+0],q[j*4+1],q[j*4+2],q[j*4+3])
                    t = vtk.vtkTransform()
                    t.SetMatrix((m[0],m[3],m[6], x[j*3+0],
                                 m[1],m[4],m[7], x[j*3+1],
                                 m[2],m[5],m[8], x[j*3+2],
                                 0.0, 0.0, 0.0,  1.0));
                    pActors[j].SetUserTransform(t)
            else:
                # no quaternion in the result
                pData = make_pData(nc.np, x, a)
                pGlyph.SetInput(pData)
    
            if lattice[0] == 0.0 and lattice[1] == 0.0 and lattice[2] == 0.0:
                # non-periodic boundary
                (cx,cy,cz, lx,ly,lz) = bounding_box (nc.np, x)
                if lx > lz:
                    l = lx
                else:                  
                    l = lz
                # prevent to go far away
                if l > 50: l = 50
    
                aCamera.SetFocalPoint (cx, cy,       cz)
                aCamera.SetPosition   (cx, cy-3.0*l, cz+0.1*l)
    
            win.Render()


if __name__ == "__main__":
    main()
