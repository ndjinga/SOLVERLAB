#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


__doc__ = """\
utilities to test, try and get elementary SALOME API
"""

import sys
import os

def hasSalomeEnv():
  """test if a salome env is set"""
  kernel = os.getenv("KERNEL_ROOT_DIR") #why not?
  if kernel == None: return False
  if os.path.isdir(kernel): return True #bad env
  return False
  
def getStudy():
  if not hasSalomeEnv(): return None
  try:
    # imports pour Salome
    import salome
    return salome.myStudy
  except:
    return None
    
def getGeompy():
  if not hasSalomeEnv(): return None
  try:
    # imports pour Salome
    import salome
    salome.salome_init()
    from salome.geom import geomBuilder
    geompy = geomBuilder.New()
    return geompy
  except:
    return None


  
