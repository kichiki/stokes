# utility routines of VTK
# Copyright (C) 2006-2008 Kengo Ichiki <kichiki@users.sourceforge.net>
# $Id: ryuon_vtk.py,v 1.2 2008/05/08 03:08:43 kichiki Exp $
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

#import sys
#sys.path.append('/somewhere/vtk-python-binding/lib64/python2.4/site-packages')
import vtk
#from vtk.util.colors import peacock, tomato
from vtk.util.colors import *


# make an Actor for mobile particles at (0,0,0) with radius 1
# this is for the result with quaternion
def make_pActor_with_quaternion():
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


def make_cubeActor (lattice, radius=0.1):
    CubeModel = vtk.vtkCubeSource()
    CubeModel.SetXLength(lattice[0])
    CubeModel.SetYLength(lattice[1])
    CubeModel.SetZLength(lattice[2])
    CubeModel.SetCenter(.5*lattice[0], .5*lattice[1], .5*lattice[2])
    Edges = vtk.vtkExtractEdges()
    Edges.SetInput(CubeModel.GetOutput())
    Tubes = vtk.vtkTubeFilter()
    Tubes.SetInput(Edges.GetOutput())
    Tubes.SetRadius(radius)
    Tubes.SetNumberOfSides(6)
    Tubes.UseDefaultNormalOn()
    Tubes.SetDefaultNormal(.577, .577, .577)
    # Create the mapper and actor to display the cube edges.
    TubeMapper = vtk.vtkPolyDataMapper()
    TubeMapper.SetInput(Tubes.GetOutput())
    CubeEdges = vtk.vtkActor()
    CubeEdges.SetMapper(TubeMapper)
    CubeEdges.GetProperty().SetDiffuseColor(khaki)
    CubeEdges.GetProperty().SetSpecular(.4)
    CubeEdges.GetProperty().SetSpecularPower(10)
    #ren.AddActor(CubeEdges)
    return(CubeEdges)


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


# make pData for np particles under periodic boundary condition
def make_pData_periodic (np, x, a, lattice):
    pos = vtk.vtkPoints()
    diameter = vtk.vtkDoubleArray()
    diameter.SetNumberOfComponents(1)
    # primary cell
    for i in range(np):
        pos.InsertNextPoint(x[i*3], x[i*3+1], x[i*3+2])
        if a != []:
            diameter.InsertNextTuple1(2.0*a[i])
        else:
            diameter.InsertNextTuple1(2.0)
    # image cells
    for ix in range(-1,2):
        for iy in range(-1,2):
            for iz in range(-1,2):
                if ix == 0 and iy == 0 and iz == 0: continue
                for i in range(np):
                    pos.InsertNextPoint(x[i*3  ]+float(ix)*lattice[0],
                                        x[i*3+1]+float(iy)*lattice[1],
                                        x[i*3+2]+float(iz)*lattice[2])
                    if a != []:
                        diameter.InsertNextTuple1(2.0*a[i])
                    else:
                        diameter.InsertNextTuple1(2.0)

    # first make pData containing particle coordinates
    pData = vtk.vtkPolyData()
    pData.SetPoints(pos)
    pData.GetPointData().SetScalars(diameter)

    return pData


# make bData for bond
def make_bData (np, x):
    points = vtk.vtkPoints()
    lines = vtk.vtkCellArray()
    lines.InsertNextCell(np)
    for i in range(np):
        points.InsertNextPoint(x[i*3], x[i*3+1], x[i*3+2])
        lines.InsertCellPoint(i)

    # first make pData containing particle coordinates
    bData = vtk.vtkPolyData()
    bData.SetPoints(points)
    bData.SetLines(lines)

    return bData

# INPUT
#  rgb = (r,g,b)
#  opacity
#  res : phi- and theta-resolutions for vtkSphereSource. (default: 8)
def make_pActor (rgb, opacity, res):
    pSource = vtk.vtkSphereSource()
    pSource.SetPhiResolution(res)
    pSource.SetThetaResolution(res)

    pGlyph = vtk.vtkGlyph3D()
    pGlyph.SetSource(pSource.GetOutput())
    #pGlyph.SetInput(pData)
    pGlyph.ScalingOn()
    pGlyph.SetScaleModeToScaleByScalar()

    pMapper = vtk.vtkPolyDataMapper()
    pMapper.ScalarVisibilityOff()
    pMapper.SetInput(pGlyph.GetOutput())

    pActor = vtk.vtkActor()
    pActor.SetMapper(pMapper)
    # put Red color
    pActor.GetProperty().SetColor(rgb[0], rgb[1], rgb[2])
    pActor.GetProperty().SetOpacity(opacity)

    return (pActor,pGlyph)


# INPUT
#  rgb = (r,g,b)
#  opacity
#  res : phi- and theta-resolutions for vtkSphereSource. (default: 8)
def make_pActor_periodic (rgb, opacity, res, lattice):
    pSource = vtk.vtkSphereSource()
    pSource.SetPhiResolution(res)
    pSource.SetThetaResolution(res)

    pGlyph = vtk.vtkGlyph3D()
    pGlyph.SetSource(pSource.GetOutput())
    #pGlyph.SetInput(pData)
    pGlyph.ScalingOn()
    pGlyph.SetScaleModeToScaleByScalar()

    pMapper = vtk.vtkPolyDataMapper()
    pMapper.ScalarVisibilityOff()
    pMapper.SetInput(pGlyph.GetOutput())

    # periodic boundary
    pNormals = vtk.vtkPolyDataNormals()
    pNormals.SetInput(pGlyph.GetOutput())

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

    plane2 = vtk.vtkPlane()
    plane2.SetOrigin(0, 0, 0)
    plane2.SetNormal(1, 0, 0)
    clipper2 = vtk.vtkClipPolyData()
    clipper2.SetInput(clipper.GetOutput())
    clipper2.SetClipFunction(plane2)
    clipper2.GenerateClipScalarsOn()
    clipper2.GenerateClippedOutputOn()
    clipper2.SetValue(0)

    plane3 = vtk.vtkPlane()
    plane3.SetOrigin(0, 0, 0)
    plane3.SetNormal(0, 0, 1)
    clipper3 = vtk.vtkClipPolyData()
    clipper3.SetInput(clipper2.GetOutput())
    clipper3.SetClipFunction(plane3)
    clipper3.GenerateClipScalarsOn()
    clipper3.GenerateClippedOutputOn()
    clipper3.SetValue(0)

    plane4 = vtk.vtkPlane()
    plane4.SetOrigin(0, lattice[1], 0)
    plane4.SetNormal(0, -1, 0)
    clipper4 = vtk.vtkClipPolyData()
    clipper4.SetInput(clipper3.GetOutput())
    clipper4.SetClipFunction(plane4)
    clipper4.GenerateClipScalarsOn()
    clipper4.GenerateClippedOutputOn()
    clipper4.SetValue(0)

    plane5 = vtk.vtkPlane()
    plane5.SetOrigin(lattice[0], 0, 0)
    plane5.SetNormal(-1, 0, 0)
    clipper5 = vtk.vtkClipPolyData()
    clipper5.SetInput(clipper4.GetOutput())
    clipper5.SetClipFunction(plane5)
    clipper5.GenerateClipScalarsOn()
    clipper5.GenerateClippedOutputOn()
    clipper5.SetValue(0)

    plane6 = vtk.vtkPlane()
    plane6.SetOrigin(0, 0, lattice[2])
    plane6.SetNormal(0, 0, -1)
    clipper6 = vtk.vtkClipPolyData()
    clipper6.SetInput(clipper5.GetOutput())
    clipper6.SetClipFunction(plane6)
    clipper6.GenerateClipScalarsOn()
    clipper6.GenerateClippedOutputOn()
    clipper6.SetValue(0)

    clipMapper = vtk.vtkPolyDataMapper()
    #clipMapper.SetInput(clipper.GetOutput())
    #clipMapper.SetInput(clipper2.GetOutput())
    #clipMapper.SetInput(clipper3.GetOutput())
    clipMapper.SetInput(clipper6.GetOutput())
    clipMapper.ScalarVisibilityOff()

    pActor = vtk.vtkActor()
    pActor.SetMapper(clipMapper)

    pActor.GetProperty().SetColor(rgb[0], rgb[1], rgb[2])
    pActor.GetProperty().SetOpacity(opacity)

    return (pActor,pGlyph)


def make_bondActor (color):
    bond = vtk.vtkTubeFilter()
    bond.SetNumberOfSides(8)
    # bond.SetInput(bData)
    bond.SetRadius(0.2)

    bondMapper = vtk.vtkPolyDataMapper()
    bondMapper.SetInput(bond.GetOutput())

    bondActor = vtk.vtkActor()
    bondActor.SetMapper(bondMapper)
    bondActor.GetProperty().SetDiffuseColor(khaki)
    bondActor.GetProperty().SetSpecular(.4)
    bondActor.GetProperty().SetSpecularPower(10)

    return (bondActor,bond)


def make_axesActor():
    axes = vtk.vtkAxes()
    axes.SetOrigin(0,0,0)
    axes.SetScaleFactor(4)
    axesTubes = vtk.vtkTubeFilter()
    axesTubes.SetInput(axes.GetOutput())
    axesTubes.SetRadius(axes.GetScaleFactor()/25.0)
    axesTubes.SetNumberOfSides(6)
    axesMapper = vtk.vtkPolyDataMapper()
    axesMapper.SetInput(axesTubes.GetOutput())
    axesActor = vtk.vtkActor()
    axesActor.SetMapper(axesMapper)
    
    #ren.AddActor(axesActor)
    return (axesActor)


def make_textActor ():
    textActor = vtk.vtkTextActor()
    #textActor.SetInput (text)
    textActor.ScaledTextOn()
    textActor.SetDisplayPosition(5, 5)
    # Set coordinates to match the old vtk.vtkScaledTextActor default value
    textActor.GetPosition2Coordinate().SetCoordinateSystemToNormalizedViewport()
    textActor.GetPosition2Coordinate().SetValue(0.6, 0.1)
    tprop = textActor.GetTextProperty()
    tprop.SetFontSize(18)
    tprop.SetFontFamilyToArial()
    tprop.SetJustificationToCentered()
    tprop.BoldOn()
    tprop.ItalicOn()
    tprop.ShadowOn()
    #tprop.SetColor(0, 0, 1)
    return (textActor)
    

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


def bounding_box_ (np, x):
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

    #lx = lx1 - lx0
    #ly = ly1 - ly0
    #lz = lz1 - lz0

    #return (cx,cy,cz, lx,ly,lz)
    return (cx,cy,cz, lx0,ly0,lz0, lx1,ly1,lz1)
