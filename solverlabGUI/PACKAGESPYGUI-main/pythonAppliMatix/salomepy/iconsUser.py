#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


"""
IconsUser assume indirection from a simple name to local ./image/namefile.png 
{key=simplename : value=namefile} (without .png extension)
to have global generic icons image files (as files .png)

others local packages python may have local icons image files in
otherpackagepy/resources/*.png

simply loaded throught $PYTHONPATH with getIconFromName()
example: 
  #from file ./image/paraview64x64.png
  aIcon = getIconFromName("paraview") #from ./image/paraview64x64.png
  #from file .../otherpackagepy/resources/paraview.png
  aIcon = getIconFromName("otherpackagepy.resources.paraview")

see http://doc-snapshot.qt-project.org/4.8/widgets-icons.html
"""

import os
import sys
from PyQt5 import QtGui, QtWidgets
import xyzpy.loggingXyz as LOG
import fnmatch

logger = LOG.getLogger()

verbose = False

IconsUser={
  "paraviewOfficial" : "paraviewOfficial64x64", #not fun
  "paraview" : "paraview64x64",
  "paraviewPatternFile" : "paraview64x64",
  "castem" : "castemMod64x64", #"earthquake100x100",
  "castemDgibiFile" : "castemMod64x64", #"earthquake100x100",
  "castemSauvFile" : "castemMod64x64",
  "medFile" : "salomeLampe180x180",
  "csvFile" : "butblue",
  "inpFile" : "butgreen",
  "pngFile" : "butpurple",
  "test" : "working64x64",
  "tests" : "working64x64",
  "noIcon" : "noIcon64x64",
  "help" : "helpred",
  "editor": "kate",
  "edit": "kate",
  "clear": "clearModel",
  "run" : "run",
  "runstandalone" : "runstandalone",
  "browsefile": "browsefile",
  "webbrowser": "webbrowser",
  "root": "root-cyan-128",
  "open" : "open",
  "save" : "save",
  "openxml" : "openxml",
  "savexml" : "savexml",
  "opencsv" : "opencsv",
  "savecsv" : "savecsv",
  "showHideColumn" : "showHideColumn",
  "listresult": "listresult",
  "TODO": "working64x64",
  "datafromfilemed": "datafromfilemed",
  
  "refresh" : "refresh",
  "clearModel" : "clearModel",
  "datainformation" : "infoblue",
  "dialogwidget" : "widget",
  "deleteitem" : "trash-empty",
  "deletechilditem" : "trash-empty",
  "insertitemup" : "arrowub",
  "insertitembelow" : "arrowdb",
  "deleteitemup" : "arrowur",
  "deleteitembelow" : "arrowdr",
  "resetinexistingchilditem" : "refresh",
  "expand" : "expand",
  "collapse" : "collapse",
  "sablier": "sablier",
  "cpp": "cpp",
  "search": "search",
  "python": "python",
  "tgz": "tgz",
  "user": "user",
  "doc": "doc",
  #"" : "",
  #cnf crescendo
  "opencnf" : "opencnf",
  "savecnf" : "savecnf",
  "crescendo" : "crescendo",
  "testCrescendoCode" : "crescendoworking64x64",
  #eki ekinox
  "ekinox" : "EKINOX",
  "openeki" : "openeki",
  "saveeki" : "saveeki",
  "testEkinoxCode" : "ekinoxworking",
  "testNumodisCode" : "numodisworking",
  #nmd numodis
  "numodis" : "NUMODIS",
  "opennmdscratch" : "opennmdscratch",
  "opennmd" : "opennmd",
  "savenmd" : "savenmd",
  #oly olympe
  "irfu" : "irfu_icon",
  "irfu_icon" : "irfu_icon",
  "openolyscratch" : "irfu_icon",
  #"openoly" : "opennmd",
  #"saveoly" : "savenmd",
  #mcg microgen
  "microgen" : "MICROGEN",
  "openmcgscratch" : "openmcgscratch",
  #"savemcg" : "savemcg",
  "testmicrogencode" : "microgenworking",
  "voronoimcg" : "voronoimcg",
  "combsmcg" : "combsmcg",
  "combsadvancedmcg" : "advanced",
  "voronoigeometrymcg": "cube",
  "voronoiinputpointsmcg": "points",
  "voronoimeshmcg": "ExecMESH",
  "voronoistatmcg": "ExecSTAT",
  "voronoiadvancedmcg": "advanced",
  "combsmeshmcg": "isocahedre",
  "combsassembleoptionsmcg": "glue",
  "combsmeshoptionsmcg": "isocahedre",
  "combscropoptionsmcg": "scissors",
  "combsstatisticoptionsmcg": "ExecSTAT",
  "combsspheremcg": "combsspheremcg",
  "combsmultiLayerspheremcg": "combsmultiLayerspheremcg",
  "geometry": "geometry",
  "ebsdmcg": "voronoi2D",
  "ebsdcleaningmcg": "ExecNETTOYAGE",
  "ebsdmicro2Dmcg": "ExecEBSD",
  "ebsdvoxelizemcg": "ExecVoxelize",
  "ebsdsurf2medmcg": "disketteRouge",
  #"": "",
  #amt amitex
  "amitex" : "AMITEX",
  "openamt" : "openamt",
  "openamtscratch" : "openamtscratch",
  "saveamt" : "saveamt",
  "testAmitexCode" : "amitexworking",
  "microstructureamt": "microstructureamt",
  "listmicrostructuresamt": "microstructureamt",
  "materiauxamt": "materiauxamt",
  "listmateriauxamt": "materiauxamt",
  "chargementamt": "haltere",
  "listchargementsamt": "haltere",
  "algorithmeamt": "algorithmeamt",
  "listalgorithmesamt": "algorithmeamt",
  #iradina
  "simulationira": "run",
  "ionbeamira": "projectiledrt",
  "targetira": "target",
  "structureira": "voronoi2D",
  "listmaterialira": "materiauxamt",
  "materialira": "materiauxamt",
  "caseira": "run",
  "browseelement": "periodictable",
  #dart
  "dart" : "DART",
  "simulationdrt": "run",
  "projectiledrt": "projectiledrt",
  "targetdrt": "target",
  "listtargetcomponentsdrt": "target",
  "targetcomponentdrt": "target",
  "algorithmedrt": "algorithmeamt",
  "listresultdrt": "listresult",
  "opendrt" : "opendrt",
  "opendrtscratch" : "opendrtscratch",
  "savedrt" : "drivesavedrt",
  "testDartCode" : "dartworking64x64",
  "isotopedrt": "isotope",
  #cmdc
  "cmdc" : "CMDC",
  "opencmdc" : "opencmdc",
  "opencmdcscratch" : "noIcon64x64", #"opensmdcscratch",
  "savecmdc" : "noIcon64x64", #"savecmdc",
  "testCmdcCode" : "noIcon64x64", #"Cmdcworking",
  #plot
  "plot" : "plot",
  "plot2" : "plot2",
  "plotn" : "plotn",
  "plotd" : "plotd",
  "ploth" : "ploth",
  "plotEllipse" : "working64x64",
  "addColumn" : "addColumn",
  #numodis
  "controlnmd" : "ExecControl",
  "listcontrolnmd" : "ExecControl",
  "dislocationnmd" : "ExecDislocation",
  "listdislocationnmd" : "ExecDislocation",
  "grainnmd" : "ExecGrain",
  "listgrainnmd" : "ExecGrain",
  "topologynmd" : "ExecTopology",
  "listtopologynmd" : "ExecTopology",
  "visualizenmd" : "ExecVisualize",
  "importnmd" : "ExecImport",
  "precipitatenmd" : "ExecPrecipitate",
  "listprecipitatenmd" : "ExecPrecipitate",
  "newnmd" : "ExecNew",
  "exportnmd" : "ExecExport",
  #iradinagui
  "iradinagui" : "ion-cannon-blast",
  "helpIraCode" : "helpgreen",
  "helpIraGui" : "helpblue",
  "helpCode" : "helpgreen",
  "helpGui" : "helpblue",
  }


IconsUser_fnpatterns=[
  ("*.sauv", "castemSauvFile"),
  ("*.dgibi", "castemDgibiFile"),
  ("*PV_*.py", "paraviewPatternFile"),
  ("*.med", "medFile"),
  ("*.csv", "csvFile"),
  ("*.inp", "inpFile"),
  ("*.png", "pngFile"),
  ("*.test", "test"),
  #("", ),
  ]

IconsUserloaded = {} #load only one time, in memory
IconPath = os.path.join(os.path.dirname(__file__),"images")

logger.debug( "iconUser path '%s'" % IconPath )

def _getIcon(name):
  #find as name=path/name(+".png" or ".gif")
  res = None
  for ext in [".png", ".gif", ".xpm"]:
    aname = name + ext
    if os.path.isfile(aname):
      if aname in list(IconsUserloaded.keys()):
       return IconsUserloaded[aname]
      try:
        res = QtGui.QIcon(aname)
        IconsUserloaded[aname] = res
        return res
      except:
        res=None
        #could be an icon but empty if file not found
    else:
      if verbose: logger.warning("No file icon for '%s'" % aname)
    
  if res==None:
    if verbose: logger.warning("No icon for '%s.*'" % name)
    aname = os.path.join(os.path.dirname(__file__), "images", "noIcon64x64.png")
    if aname in IconsUserloaded:
      res = IconsUserloaded[aname]
    else:
      res = QtGui.QIcon(aname)
      IconsUserloaded[aname] = res
  return res

def findInPythonPath(text):
  """
  from text as import syntax (like "otherpackagepy.resources.paraview")
  for file /.../otherpackagepy/resources/paraview.png
  try import relative directory (throught $PYTHONPATH)
  for validate resources directory (as usage here)
  """

  """
  example:
    export PYTHONPATH=/export/home/wambeke/SALOME-7.8.0_CV/SOURCES/CASSIS/src
    python
    import cassispy.resources 
  """

  import os
  import traceback
  if verbose: logger.info("findInPythonPath text '%s'" % text)
  s = text.split(".")            #"cassispy.resources.cassis".split(".")
  if len(s) <= 1: return None
  aPath = os.path.join(*s[1:-1]) #resources
  aImport = ".".join(s[0:-1])    #cassispy.resources
  aName = s[-1]
  try:
    aModule = __import__(aImport, globals(), locals(), []) 
    tmp = os.path.realpath(aModule.__file__)  #/.../cassispy/__init__.pyc
    tmp = os.path.dirname(tmp)                #/.../cassispy
    localIconFileName = os.path.join(tmp, aPath, aName)
    if verbose: logger.warning("aImport %s as %s" % (aImport, localIconFileName))
    return localIconFileName
  except Exception as e:
    #trace = traceback.format_exc() #better explicit verbose problem
    #raise Exception("%s\n%s\ncan't import file: '%s.png'" % (trace, e, text))
    if verbose: logger.warning("findInPythonPath: No icon for '%s.png'" % text)
    return None


def getIconFromName(text):
  """load only one time, store in IconsUserloaded, and resend same instance (same id)"""
  path = IconPath
  if text in list(IconsUser.keys()): #something like aIcon = getIconFromName("paraview")
    afile = os.path.join(path,IconsUser[text])
    return _getIcon(afile)

  afile = findInPythonPath(text) #something like aIcon = getIconFromName("otherpackagepy.resources.paraview")
  if afile != None:
    return _getIcon(afile)

  afile = os.path.join(path,IconsUser["noIcon"])
  if verbose: logger.warning("No icon for '%s'" % text)
  return _getIcon(afile)

  
def getIconFileName(text):
  if text in list(IconsUser.keys()):
    res = IconsUser[text] + ".png"
  else:
    res = IconsUser["noIcon"] + ".png"
  return res

def getIconFromPattern(text):
  name = None
  for pat, nam in IconsUser_fnpatterns:
    if fnmatch.fnmatch(text, pat):
      name = nam
      break
  if name == None:
    return None
  else:
    return getIconFromName(name)

"""
#sert a rien... resoudre pb par
#from salomepy.onceQApplication import OnceQApplication
#app = OnceQApplication()

#sert a rien...
try:
  print "test entry Problem Must construct a QApplication before a QPaintDevice"
  getIconFromName("noIcon")
  print "test end ok Problem Must construct a QApplication before a QPaintDevice"
except:
  #erreur non catchee core dumped
  print "Detect Problem Must construct a QApplication before a QPaintDevice"
"""
