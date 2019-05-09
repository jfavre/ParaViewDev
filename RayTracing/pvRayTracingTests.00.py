# Paraview master
#
# Testing raytracing and pathtracing
# https://www.physicsclassroom.com/class/refrn/Lesson-1/Refraction-and-Sight
#
# Jean M. Favre, Swiss National Supercomputing Centre
# Thu May  9 09:37:03 CEST 2019


#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

mf1 = "/users/jfavre/Projects/ParaView/ospray_mats.json"
materialLibrary1 = GetMaterialLibrary()
print("using materials: {:}".format(mf1))
materialLibrary1.LoadMaterials = mf1

Batch = True
if Batch:
  ProgressivePasses = 0
else:
  ProgressivePasses = 5
  
SamplesPerPixel = 20

# Create a new 'Light'
light1 = CreateLight()
light1.Position = [50.0, 50.0, -150.0]

# get the material library
materialLibrary1 = GetMaterialLibrary()

# create light
# Create a new 'Render View'
renderView1 = GetRenderView()
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [0.0, 0.4856886910475378, 0.0]
renderView1.UseLight = 0
renderView1.StereoType = 0
renderView1.CameraPosition = [-2.0845170131514843, 1.251255184574042, 1.885297086960957]
renderView1.CameraFocalPoint = [38.85293334447471, -12.375838003738075, -31.700665545414243]
renderView1.CameraViewUp = [0.20406248231693128, 0.9682879470153982, -0.14414213462883688]
renderView1.CameraParallelScale = 14.151484590759953
renderView1.Background = [0.1803921568627451, 0.20392156862745098, 0.21176470588235294]
renderView1.EnableRayTracing = 0 # Classic OpenGL
renderView1.AdditionalLights = light1
renderView1.OSPRayMaterialLibrary = materialLibrary1

imageRes = [1024, 1024]
renderView1.ViewSize = imageRes
SetActiveView(renderView1)

# create a new 'Programmable Source'
programmableSource1 = ProgrammableSource(guiName="Mirrors")
programmableSource1.OutputDataSetType = 'vtkUnstructuredGrid'
programmableSource1.Script = """
import vtk
from vtk import VTK_WEDGE
import numpy as np
from vtk.numpy_interface import dataset_adapter as dsa
executive = self.GetExecutive()
outInfo = executive.GetOutputInformation(0)
XYZ = np.array([ 1.0, 0.0, -1.0,
               3.0, 0.0, -5.0,
               2.0, 0.0,  5.0,
               1.0,  1.0, -1.0,
               3.0,  1.0, -5.0,
               2.0,  1.0,  5.0 ])
nnodes = XYZ.shape[0]//3
CONNECTIVITY = np.array([6, 0,1,2,3,4,5])
nelts = 1
CELL_TYPES = np.full((nelts), VTK_WEDGE, np.ubyte)
CELL_OFFSETS = np.arange(nelts)
CELL_OFFSETS = 0 * CELL_OFFSETS
output.SetCells(CELL_TYPES, CELL_OFFSETS, CONNECTIVITY)
output.Points = XYZ.reshape((nnodes,3))
"""
programmableSource1.ScriptRequestInformation = ''
programmableSource1.PythonPath = ''

# create a new 'Cylinder'
cylinder2 = Cylinder()
cylinder2.Resolution = 64
cylinder2.Height = 0.75
cylinder2.Radius = 0.015
cylinder2.Center = [0.0, 0.35, 0.0]

# create a new 'Cylinder'
cylinder1 = Cylinder()
cylinder1.Resolution = 64
cylinder1.Height = 0.5
cylinder1.Radius = 0.25
cylinder1.Center = [0.0, 0.25, 0.0]

# create a new 'Plane'
plane1 = Plane()
plane1.Origin = [-10.0, 0.0, -10.0]
plane1.Point1 = [-10.0, 0.0, 10.0]
plane1.Point2 = [10.0, 0.0, -10.0]
plane1.XResolution = 200
plane1.YResolution = 200

# show data from plane1
plane1Display = Show(plane1, renderView1)

# trace defaults for the display properties.
plane1Display.Representation = 'Surface With Edges'
plane1Display.ColorArrayName = ['POINTS', '']
plane1Display.LineWidth = 0.1
plane1Display.Ambient = 0.2
plane1Display.Diffuse = 0.8

# show data from cylinder1
cylinder1Display = Show(cylinder1, renderView1)

# trace defaults for the display properties.
cylinder1Display.Representation = 'Surface'
cylinder1Display.AmbientColor = [0.0, 0.0, 0.0]
cylinder1Display.ColorArrayName = ['POINTS', '']
cylinder1Display.Ambient = 1.0
cylinder1Display.Diffuse = 0.8

# show data from cylinder2
cylinder2Display = Show(cylinder2, renderView1)

# trace defaults for the display properties.
cylinder2Display.Representation = 'Surface'
cylinder2Display.ColorArrayName = [None, '']
cylinder2Display.Orientation = [-20.0, 0.0, 0.0]

# show data from programmableSource1
programmableSource1Display = Show(programmableSource1, renderView1)

# trace defaults for the display properties.
programmableSource1Display.Representation = 'Surface With Edges'
programmableSource1Display.ColorArrayName = [None, '']

# First do classic OpenGL rendering with Ambient Light only
programmableSource1Display.Ambient = 1.0
programmableSource1Display.Diffuse = 0.0
if Batch:
  SaveScreenshot("RayTracing-tutorial.00.png", ImageResolution=imageRes)

# Second do classic OpenGL rendering with Ambient and Diffuse light
programmableSource1Display.Ambient = 0.2
programmableSource1Display.Diffuse = 0.8
if Batch:
  SaveScreenshot("RayTracing-tutorial.01.png", ImageResolution=imageRes)


# Third do OSPRay rendering without shadows
renderView1.EnableRayTracing = 1
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.LightScale = 1.2
renderView1.Shadows = 0
renderView1.AmbientSamples = 1
renderView1.Denoise = 1
renderView1.ProgressivePasses = ProgressivePasses  # must launch with --enable-streaming
  
# try values between 10 and 200 to eliminate noise
renderView1.SamplesPerPixel = SamplesPerPixel

if Batch:
  SaveScreenshot("RayTracing-tutorial.02.png", ImageResolution=imageRes)

# Fourth do OSPRay rendering with shadows
renderView1.Shadows = 1
if Batch:
  SaveScreenshot("RayTracing-tutorial.03.png", ImageResolution=imageRes)

# Fifth do OSPRay rendering with *soft* shadows
renderView1.BackEnd = 'OSPRay pathtracer'
light1.Radius = 5
if Batch:
  SaveScreenshot("RayTracing-tutorial.04.png", ImageResolution=imageRes)

# Sixth do OSPRay rendering with shadows and a reflective material
programmableSource1Display.OSPRayMaterial = 'copper'
if Batch:
  SaveScreenshot("RayTracing-tutorial.05.png", ImageResolution=imageRes)

# Seventh do OSPRay rendering. A wood material
cylinder1Display.OSPRayMaterial = 'wood'
if Batch:
  SaveScreenshot("RayTracing-tutorial.06.png", ImageResolution=imageRes)

# Eigth do OSPRay rendering. A refractive material
cylinder1Display.OSPRayMaterial = 'water'
if Batch:
  SaveScreenshot("RayTracing-tutorial.07.png", ImageResolution=imageRes)

# Ninth do OptiX rendering with *soft* shadows
renderView1.BackEnd = 'OptiX pathtracer'
cylinder1Display.OSPRayMaterial = 'None'
programmableSource1Display.OSPRayMaterial = 'None'
if Batch:
  SaveScreenshot("RayTracing-tutorial.08.png", ImageResolution=imageRes)

# Tenth do OptiX rendering with shadows and a reflective material
programmableSource1Display.OSPRayMaterial = 'copper'
if Batch:
  SaveScreenshot("RayTracing-tutorial.09.png", ImageResolution=imageRes)

# Eleventh OptiX rendering. A wood material
cylinder1Display.OSPRayMaterial = 'wood'
if Batch:
  SaveScreenshot("RayTracing-tutorial.10.png", ImageResolution=imageRes)

# 12-th do Optix rendering. A refractive material
cylinder1Display.OSPRayMaterial = 'water'
if Batch:
  SaveScreenshot("RayTracing-tutorial.11.png", ImageResolution=imageRes)





