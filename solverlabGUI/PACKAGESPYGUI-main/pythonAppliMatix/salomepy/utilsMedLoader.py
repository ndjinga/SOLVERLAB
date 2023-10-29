#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


__doc__ = """\
utilities to medloader
thanks to Christophe Bourcier
"""

import os
import sys
import pprint as PP #pretty print
import MEDLoader as ML

import xyzpy.loggingXyz as LOG

logger = LOG.getLogger()

########################################################################
def getMeshFromFile(fileName):
  """get mesh from med file"""
  realPath = os.path.realpath(os.path.expandvars(fileName))
  if not os.path.isfile(realPath):
    logger.warning("Inexisting file '%s' as '%s'" % (fileName, realPath))
    return None
  try:
    medmesh = ML.MEDFileUMesh.New(realPath)
    return medmesh
  except:
    logger.warning("Problem reading med file '%s'" % realPath)
    return None

########################################################################
def getGroupsFromMesh(medmesh):
  """
  get mesh groups by dimension
  return a list of tuples = [ (levelName, groupsForThisLevel), ... ]
  with levelsNames = [ Nodes, Beams, Shells, Volumes ]
  """
  groups = medmesh.getGroupsNames()
  meshDim = medmesh.getMeshDimension()
  logger.debug("Mesh dimension: %s Groups:\n%s" % (meshDim, PP.pformat(groups))
  #"noeuds", "poutres", "coques", "volumes"
  levels = [("Nodes", 0), ("Beams", 1), ("Shells", 2), ("Volumes", 3)]
  res = []
  for level, dim in levels[0:meshDim+1]:
    if dim == 0:
      # groupes de noeuds
      meshDimRelToMaxExt = 1
    else:
      # groupes de mailles: 0 (dimension max du maillage), -1 (dimension inférieure), -2 (encore inférieure)
      meshDimRelToMaxExt = dim - meshDim
    groups_dim = medmesh.getGroupsOnSpecifiedLev(meshDimRelToMaxExt) 
    logger.debug("Groups of %s: %s" % (level, PP.pformat(groups_dim)))
    res.append((level, groups_dim))
  return res


########################################################################
def testColors():
  """
  Visualization of named colors.
  Simple plot example with the named colors and its visual representation.
  display geometric groups and mesh groups with random color.
  """
  import six
  import numpy as np
  import matplotlib.pyplot as plt
  from matplotlib import colors

  colors_ = list(six.iteritems(colors.cnames))

  # Add the single letter colors.
  for name, rgb in six.iteritems(colors.ColorConverter.colors):
      hex_ = colors.rgb2hex(rgb)
      colors_.append((name, hex_))

  # Transform to hex color values.
  hex_ = [color[1] for color in colors_]
  # Get the rgb equivalent.
  rgb = [colors.hex2color(color) for color in hex_]
  return


########################################################################
def testCreateAndDisplayBox():
  import time
  import math
  import salome
  import random

  salome.salome_init()
  theStudy = salome.myStudy

  import GEOM
  from salome.geom import geomBuilder
  import SALOMEDS

  geompy = geomBuilder.New()


  import SalomePyQt
  sgPyQt=SalomePyQt.SalomePyQt()

  ## Refresh the GUI in 7.5.1 (synchronous python console)
  def processEvents():
    try:
      sgPyQt.processEvents()
    except:
      # Not in 7.5.1: nothing to do
      pass

  # Create OCC view
  the_view = sgPyQt.getActiveView()
  if the_view == -1 or sgPyQt.getViewType(the_view) != "OCCViewer":
      sgPyQt.createView("OCCViewer")
      processEvents()


  nb_spheres = 50
  r_min = 1
  r_max = 5

  x_max = 100
  y_max = 100
  z_max = 100


  geom_gui = salome.ImportComponentGUI("GEOM")

  # Create a box
  box = geompy.MakeBoxDXDYDZ(x_max, y_max, z_max)
  entry = geompy.addToStudy(box, "box")

  # Display the box in wireframe
  geom_gui.setDisplayMode(entry, 0)
  geom_gui.createAndDisplayGO(entry)
  processEvents()

  random.seed()

  # Place some spheres randomly inside the box
  spheres = []
  colors = []
  for i in range(nb_spheres):
      x = random.uniform(0, x_max)
      y = random.uniform(0, y_max)
      z = random.uniform(0, z_max)
      r = random.uniform(r_min, r_max)
      sphere = geompy.MakeSphere(x, y, z, r)
      entry = geompy.addToStudy(sphere, "sphere_%i"%(i+1))

      color = [random.randint(0, 255), random.randint(0, 255), random.randint(0, 255)]
      colors.append(color)

      # Display sphere in random color, shading mode
      geom_gui.setDisplayMode(entry, 1)
      geom_gui.setColor(entry, *color)
      geom_gui.createAndDisplayGO(entry)
      processEvents()
      
      spheres.append(sphere)

  logger.info("MakePartition...")
  processEvents()
  time0 = time.time()

  box_part = geompy.MakePartition([box], spheres)
  time1 = time.time()
  entry = geompy.addToStudy(box_part, "box_part")

  logger.info("done in %.3f s." % (time1-time0))
  processEvents()

  geom_gui.setDisplayMode(entry, 0)
  geom_gui.createAndDisplayGO(entry)

  logger.info("Getting sub-shapes...")
  processEvents()
  time0 = time.time()

  sub_shapes = []
  for sphere in spheres:
    name = sphere.GetName()
    sub_sphere = geompy.GetInPlace(box_part, sphere, 1)
    geompy.addToStudyInFather(box_part, sub_sphere, name)
    sub_shapes.append(sub_sphere)
    
  time1 = time.time()
  logger.info("done in %.3f s."%(time1-time0)) 
  processEvents()

  # update the object browser
  salome.sg.updateObjBrowser(0)
  processEvents()

  import  SMESH, SALOMEDS
  from salome.smesh import smeshBuilder
  smesh_gui = salome.ImportComponentGUI('SMESH')

  logger.info("Building the mesh...")
  processEvents()
  smesh = smeshBuilder.New()

  Mesh_1 = smesh.Mesh(box_part)
  algo2d = Mesh_1.Triangle(algo=smeshBuilder.MG_CADSurf)
  params = algo2d.Parameters()
  params.SetGeometricMesh( 1 )
  params.SetAngleMesh( 8 )
  Mesh_1.Tetrahedron(algo=smeshBuilder.MG_Tetra)
  processEvents()
  time0 = time.time()
  isDone = Mesh_1.Compute()
  time1 = time.time()

  logger.info("Mesh generated in %.3f s." % (time1-time0))
  processEvents()

  # Create a VTK view
  the_view = sgPyQt.getActiveView()
  if the_view == -1 or sgPyQt.getViewType(the_view) != "VTKViewer":
      sgPyQt.createView("VTKViewer")
      processEvents()

  # Create the groups of mesh for each sphere and display them
  for i, shape in enumerate(sub_shapes):
    gr = Mesh_1.Group(shape)
    color = colors[i]
    # for SMESH, color must be in [0, 1]
    color = [c/255. for c in color]
    gr.SetColor( SALOMEDS.Color(*color))
    # Display group of sphere
    entry = salome.ObjectToID(gr)
    smesh_gui.CreateAndDisplayActor(entry)
    salome.sg.FitAll()
    processEvents()


  salome.sg.updateObjBrowser(0)


