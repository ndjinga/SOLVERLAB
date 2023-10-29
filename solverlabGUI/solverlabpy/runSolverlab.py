#!/usr/bin/env python

import sys
import os

def modify(data):
  # insert code here to modify the model
  return data

def runSolverlab(data):
  import solverlabpy.diffusionEq as DEQ
  import solverlabpy.LaunchDiffusionEq as LDEQ

  n = data.Analysis.caseSolverlab.Equation
  equation = data.Model[int(n)]
  if equation.__class__ == DEQ.DiffusionEq: # check which script to launch
    # for i in range(1):
    #   modify(data)
    LDEQ.launchDEQ(data)


def loadXML(file):
  import xyzpy.utilsXyz as UXYZ
  import solverlabpy.modelSvl

  return UXYZ.fromFileXml(file)



if __name__ == "__main__":
  # os.environ["PACKAGESPY_ROOT_DIR"] = "/volatile/catB/ym268439/packagespy"
  packagespydir = os.getenv("PACKAGESPY_ROOT_DIR", None) # check all environement variable to use the model from packagespy
  solverlabdir = os.getenv("SOLVERLABGUI_ROOT_DIR", None)
  if packagespydir is None or solverlabdir is None:
    raise ValueError("PACKAGESPY_ROOT_DIR or SOLVERLABGUI_ROOT_DIR environ variable not set, fix it.")
  sys.path.insert(0, os.path.join(packagespydir, "packagespy"))
  sys.path.insert(0, solverlabdir)

  if len(sys.argv) == 2:
    data = loadXML(sys.argv[1]) # import model data from .xml file in arg
    runSolverlab(data)
  else:
    raise ValueError("missing argument")