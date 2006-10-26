# visualization program for stokes-nc file by VTK
# Copyright (C) 2006 Kengo Ichiki <kichiki@users.sourceforge.net>
# $Id: stvis-vtk.py,v 1.2 2006/10/26 02:42:10 kichiki Exp $
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
import time
import vtk
import sys
#sys.path.append('/somewhere/ryuon/stokes/python')
import stokes

def pData_read(np, x):
    pos = vtk.vtkPoints()
    diameter = vtk.vtkDoubleArray()
    diameter.SetNumberOfComponents(1)
    for i in range(np):
        pos.InsertNextPoint(x[i*3], x[i*3+1], x[i*3+2])
        diameter.InsertNextTuple1(2.0)

    # first make pData containing particle coordinates
    pData = vtk.vtkPolyData()
    pData.SetPoints(pos)
    pData.GetPointData().SetScalars(diameter)

    return pData


def usage():
    print '$Id: stvis-vtk.py,v 1.2 2006/10/26 02:42:10 kichiki Exp $'
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
    #stokes.stokes_nc_print_actives(nc, stokes.get_stdout())
    print 'ntime = ', nc.ntime
    lattice = stokes.darray(3)
    stokes.stokes_nc_get_l (nc, lattice)
    print 'lattice = ', lattice[0], lattice[1], lattice[2]
    if lattice[0] > lattice[2]:
        l = lattice[0]
    else:                  
        l = lattice[2]

    x  = stokes.darray(nc.np  * nc.nvec)


    # then, make Actor
    pSource = vtk.vtkSphereSource()
    # pSource.SetPhiResolution(20)
    # pSource.SetThetaResolution(20)

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

    # fixed particles
    if nc.npf > 0:
        print 'npf = ', nc.npf
        xf0  = stokes.darray(nc.npf * nc.nvec)
        stokes.stokes_nc_get_data0 (nc, "xf0", xf0)
        pfData = pData_read(nc.npf, xf0)

        pfGlyph = vtk.vtkGlyph3D()
        pfGlyph.SetSource(pSource.GetOutput()) # use the same pSource
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

    
    ## prepare renderer
    ren = vtk.vtkRenderer()
    win = vtk.vtkRenderWindow()
    win.AddRenderer(ren)
    # iren = vtk.vtkRenderWindowInteractor ()
    # iren.SetRenderWindow(win)

    ren.AddActor(pActor)
    if nc.npf > 0: ren.AddActor(pfActor)
    ren.SetBackground (0.1, 0.2, 0.4)
    win.SetSize(480,360)

    aCamera = vtk.vtkCamera()
    ren.SetActiveCamera (aCamera)
    aCamera.SetFocalPoint (0.5*lattice[0], 0, 0.2*lattice[2])
    aCamera.SetPosition (0.5*lattice[0], 2.0*l, 0.1*lattice[2])

    # loop
    while 1:
        for i in range(nc.ntime):
            stokes.stokes_nc_get_data (nc, "x", i, x)

            pData = pData_read(nc.np, x)
            pGlyph.SetInput(pData)
            win.Render()


if __name__ == "__main__":
    main()


  
