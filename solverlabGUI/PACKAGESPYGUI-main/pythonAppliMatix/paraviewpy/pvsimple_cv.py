#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END

#
# Author: Adrien Bruneton (CEA)

r""" This module is a direct forward to the initial 'simple' module of ParaView.
On top of that it also establishes a connection to a valid PVServer whose address
is provided by the PVSERVER engine.
"""

"""
find . -name "pvsimple_cv.py"
cp ./SOURCES/PACKAGESPY/pythonAppliMatix/paraviewpy/pvsimple_cv.py ./INSTALL/PACKAGESPY/pythonAppliMatix/paraviewpy/pvsimple_cv.py
python ./SOURCES/PACKAGESPY/pythonAppliMatix/paraviewpy/test/paraviewDisplayFromStandaloneExample.py
"""

__DEBUG = 1   # increase if you want more verbosity

print("!!!!!pvsimple_cv!!!!!")

def __my_log(msg):
    if __DEBUG:
      print("[PARAVIS] %s" % msg)

def __getFromGUI():
    """ Identify if we are running inside SALOME's embedded interpreter.
    @return a value strictly greater than 0 if we are in SALOME's embedded interpreter
    @return 2 if we are in Salome embedded Python console.
    """
    import salome_iapp
    ret = 0
    if salome_iapp.IN_SALOME_GUI:
      ret += 1
    try:
      if __IN_SALOME_GUI_CONSOLE:  # only defined if we are in SALOME's embedded console (not only GUI)
        ret += 1
    except NameError:
      pass
    return ret

def ShowParaviewView():
    """
    If the import is made from SALOME embedded console, the ParaView application needs to
    be instanciated to avoid a future crash.
    """
    if __getFromGUI():
      __my_log("Initializing ParaView main elements, please be patient ...")
      import SalomePyQt
      sgPyQt = SalomePyQt.SalomePyQt()
      viewIds = sgPyQt.findViews("ParaView")
      if len(viewIds):
        sgPyQt.setViewVisible(viewIds[0], True)
        sgPyQt.activateView(viewIds[0])
      else:
        sgPyQt.createView("ParaView")
      # Now let the GUI main loop process the initialization event posted above
      sgPyQt.processEvents()
      __my_log("ParaView initialized.")

## The below has to called BEFORE importing paraview!!! This is crazy, but it has to be.
#ShowParaviewView()

import paraview
import pvserver
from paraview import simple

def SalomeConnectToPVServer():
    """
    Automatically connect to the right PVServer when not ("inside SALOME GUI" and "already connected").
    """
    __my_log("Connecting to PVServer ...")
    server_url = ""
    if True: #try:
        isGUIConnected = pvserver.myPVServerService.GetGUIConnected()
        if isGUIConnected and __getFromGUI():
            __my_log("Importing pvsimple from GUI and already connected. Won't reconnect.")
            return
        server_url = pvserver.myPVServerService.FindOrStartPVServer(0)
        # Extract host and port from URL:
        a = server_url.split(':')
        b = a[1].split('//')
        host, port = b[-1], int(a[-1])
        simple.Connect(host, port)
        __my_log("Connected to %s!" % server_url)
        if __getFromGUI():
            pvserver.myPVServerService.SetGUIConnected(True)
    else: #except Exception as e:
        __my_log("*******************************************")
        __my_log("** Could not connect to a running PVServer!")
        __my_log("*******************************************")
        raise e
    pass

if __getFromGUI() < 1:
    # Only if not in GUI (otherwise the createView will do the connection)
    SalomeConnectToPVServer()
#del SalomeConnectToPVServer

# Forward namespace of simple into current pvsimple:
for name in dir(simple):
  if not name.startswith("__"):
    globals()[name] = getattr(simple, name)
del simple
