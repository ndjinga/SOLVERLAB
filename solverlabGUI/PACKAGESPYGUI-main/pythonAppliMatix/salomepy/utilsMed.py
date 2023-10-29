#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


__doc__ = """\
utilities to med conversions
"""

import os

import xyzpy.loggingXyz as LOG

logger = LOG.getLogger()

def isSameFile(file1, file2):
  """follow relative paths and symbolic links"""
  rfile1=os.path.realpath(file1)
  rfile2=os.path.realpath(file2)
  return (rfile1==rfile2)

def med2sauv(file_in, file_out=None):
  import medutilities as MU #needs import _MEDLoader salome libmedloader.so so not in main
  if file_out==None:
    dire, base=os.path.split(file_in) 
    file_out = os.path.join(dire, os.path.splitext(base)[0] + ".sauv")
  if not isSameFile(file_in, file_out):
    MU.convert(file_in, "MED", "GIBI", file_out=file_out)
    return file_out
  else:
    logger.error("med2sauv: not done: files are same: "+file_in+" -> "+file_out)
    return None

def sauv2med(file_in, file_out=None):
  import medutilities as MU #needs import _MEDLoader salome libmedloader.so so not in main
  if file_out==None:
    dire, base=os.path.split(file_in) 
    file_out = os.path.join(dire, os.path.splitext(base)[0] + ".med")
  if not isSameFile(file_in, file_out):
    MU.convert(file_in, "GIBI", "MED", file_out=file_out)
    return file_out
  else:
    logger.error("sauv2med: not done: files are same: "+file_in+" -> "+file_out)
    return None

