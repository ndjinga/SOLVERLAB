import os
import sys

solverlabdir, _ = os.path.split(__file__)
solverlabdir, _ = os.path.split(solverlabdir)
os.environ["SOLVERLABGUI_ROOT_DIR"] = solverlabdir
sys.path.insert(0, solverlabdir)

PACKAGEPY_ROOT_DIR = os.getenv("PACKAGESPY_ROOT_DIR")
if os.getenv("USER") == "ym268439":
  os.environ["PACKAGESPY_ROOT_DIR"] = "/volatile/catB/ym268439/packagespy" # for yacine temporary TODO remove that line
  PACKAGEPY_ROOT_DIR = os.getenv("PACKAGESPY_ROOT_DIR")

if PACKAGEPY_ROOT_DIR:
  print("packagepy find")
  sys.path.insert(0, os.path.join(PACKAGEPY_ROOT_DIR, "packagespy"))
else:
  print("packagepy not find notebook may not work")



