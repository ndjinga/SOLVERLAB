#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END

"""\
example module managing paraview display on miscellanous types of files
have to be launch in SALOME (or OMNIORB_CONFIG set, see paraviewDisplayFromStandaloneExample.py)
"""

import os
import pprint as PP
import pvsimple as PV
import paraviewpy.paraviewDisplay as PVD

verbose = True
  
if verbose: print("__name__", __name__)

testDir = os.path.join(os.path.split(os.path.realpath(PVD.__file__))[0])
aVtkFile = os.path.join(testDir, "test", "voronoi_10grains_voxelized.vtk")
   
def test(aVtkFile):
  if verbose: print("PVD.paraviewDisplayMicrostructureVtk(%s)" % aVtkFile)
  PVD.paraviewDisplayMicrostructureVtk(aVtkFile)


test(aVtkFile)

if verbose: 
  print("dir():\n", PP.pformat(dir()))
  print("OMNIORB_CONFIG='%s'" % os.getenv("OMNIORB_CONFIG"))


    
