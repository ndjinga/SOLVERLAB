#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import os
from PyQt5 import QtGui, QtCore, QtWidgets
import xyzpy.loggingXyz as LOG
import salomepy.iconsUser as IUSR
import traceback

#import traceback

logger = LOG.getLogger()

verbose = False

"""
define and store commons action(s) and slot(s) and event(s)
for multiple views or widgets or other instances.
"""

__commonActions__ = {}
__doneOnce__ = [False]
_messDone = []
_inUnittest = True

class QActionXyz(QtWidgets.QAction):
  #http://pyqt.sourceforge.net/Docs/PyQt5/new_style_signals_slots.html
  #http://qt-project.org/doc/qt-5.1/qtwidgets/qaction.html ... Detailed Description
  defaultSignal = QtCore.pyqtSignal(dict) #for debug test #like **kwargs , for useful
  GeneralHelpSignal = QtCore.pyqtSignal(object)
  AllTestsSignal = QtCore.pyqtSignal(object)

  def __init__(self, name=None, parent=None, text=None):
    super(self.__class__, self).__init__(parent)
    if name != None:
      self.name = name
      self.setText(name)
    else:
      self.name = None
      self.setText("?")
    if text != None:
      self.setText(text)
    self.slot = None
    self.signal = None
    self.lastDefaultSlotArgs = [] #all last Args
    self.isConnected = False
    self.initialAction = None
    self.iconFile = ""

  def setEnable(self, status=True):
    if status:
      self.setEnabled()
    else:
      self.setDisabled()

  def setAction(self, slot=None, signal=None, shortcut=None, tooltip=None, icon=None, verbose=True):
    if shortcut != None: self.setShortcut(shortcut)
    if tooltip != None: self.setToolTip(tooltip)
    if tooltip != None: self.setStatusTip(tooltip)
    logger.debug("action '%s' with tooltip '%s'" % (self.name, tooltip))
    if icon!=None:
      self.setIcon(IUSR.getIconFromName(icon))
    self.iconFile = IUSR.getIconFileName(icon)
    if self.slot != None:
      if verbose: logger.warning("setAction : slot of common signal ''%s' is affected yet" % self.name)
      return False
    #TODO explicitely stupid if affect again same slot to same signal
    #could be protected if tests... __future__
    if slot != None: #menu and buttons are triggered (without args)
      self.triggered.connect(self.internalSlot)
      pass
    if signal != None:
      #print "dir(signal)",dir(signal),signal
      signal.connect(slot)
      self.isConnected = True
      pass
    self.slot = slot
    self.signal = signal
    return True

  def internalSlot(self, checkable):
    logger.debug("internalSlot of common signal '%s' from triggered." % self.name)
    if not self.isConnected:
      self.signal.connect(self.slot)
      self.isConnected = True
    self.signal.emit(None)
    pass

  def defaultSlot(self, args=None, verbose=False):
    #it's work: QtCore.QObject object but oriented-object
    #http://http://qt-project.org/doc/qt-5.0/qtcore/qobject.html#sender
    #Warning: This function violates the object-oriented principle of modularity.
    #However, getting access to the sender might be useful when many signals are connected to a single slot.
    #sender = self.sender()
    if verbose: logger.info("default slot of action of common name '%s', args : %s" % (self.name, args))
    self.lastDefaultSlotArgs.append(args)

  def createNewAction(self, name=None, slot=None, text=None, toolTip=None, verbose=False):
    """create (duplicate self) new QActionXyz instance, with/for another name and new slot... or not"""
    action = QActionXyz(self.name, self.parent())
    if verbose: logger.info("createNewAction from action '%s'" % self.name)
    if name != None:
      action.name = name
      action.setText(name)
    if text != None:
      action.setText(text)
    action.slot = slot
    action.signal = self.signal
    action.setShortcut(self.shortcut())
    if toolTip==None:
      action.setToolTip(self.toolTip())
    else:
      action.setToolTip(toolTip)
    action.setStatusTip(self.statusTip())
    action.setIcon(self.icon())
    action.initialAction = self
    action.triggered.connect(action.slot)
    return action

#def getCommonActions(): #for __future__
#  return __commonActions__

#def verifyCommonActions(): #for  __future__
#  """test if same signal for same slot in __commonActions__ (multiples calls to same slot for same signal)"""
#  return True

def getCommonActionByName(name):
  #global __doneOnce__
  logger.debug("for action %s", name)
  if __doneOnce__[0] == False:
    createGeneralActions() #because Must construct a QApplication before a QPaintDevice
    __doneOnce__[0] = True
  try:
    return __commonActions__[name]
  except:
    return None

def addInCommonActions(action, verbose=True):
  if action.name == None:
    if verbose: logger.error("actionFactoryXyz : add action without name : %s\n%s" % (action.name, "".join(traceback.format_stack()[-2:])))
    return False
  try:
    if type(action.name) != str:
      if verbose: logger.error("actionFactoryXyz : add action without name : %s\n%s" % (action.name, "".join(traceback.format_stack()[-2:])))
      return False
    if action.name in list(__commonActions__.keys()):
      #print __commonActions__
      if verbose:
        #mess = "existing yet action '%s'" % action.name
        mess = "actionFactoryXyz : existing yet action with name '%s' (to check if not unittest)" % action.name
        if mess not in _messDone and not _inUnittest: #avoid message warning in unittest
          logger.warning(mess)
          _messDone.append(mess)
      return False
    else:
       logger.debug("actionFactoryXyz : create action '%s'" % action.name)
       __commonActions__[action.name] = action
       return True
  except Exception as e:
    logger.error("actionFactoryXyz : exception on action : %s\n  %s" % (action.name, e))
    return False

def createGeneralActions():
  #warning Must construct a QApplication before a QPaintDevice (for creating icons)
  logger.debug("create general common actions")
  createActionGeneralHelp()

def createActionGeneralHelp():
  action = QActionXyz(name="GeneralHelp", text="Help")
  ok = action.setAction( slot=GeneralHelpAction, signal=action.GeneralHelpSignal,
                         shortcut=None, tooltip="general help", icon="help" )
  if ok:
    addInCommonActions(action)

def GeneralHelpAction(args=None):
  """
  args is a map. like {"nameUrlHelp": anUrlString}
  affich a args["nameUrlHelp"] salome browser for doc
  for example. nameUrlHelp=$MODULE_ROOT_DIR/share/doc/salome/gui/MODULE/index.html
  use xdg-open, not SalomePyQt.SalomePyQt().helpContext(maDoc,"")
  you can use this pattern to create your action everywhere in your code,
  and put it in factory with addInCommonActions
  """
  import xyzpy.utilsXyz as UXYZ
  nameBrowser = UXYZ.getBrowser() # "xdg-open" or firefox or chrome...
  try:
    nameUrlHelp = args["nameUrlHelp"]
  except:
    mess = "no nameUrlHelp in args: '%s'" % str(args)
    logger.error(mess)
    return
  if os.path.exists(nameUrlHelp):
    cmd = nameBrowser + " " + nameUrlHelp + " &"
    proc = SP.Popen(str(cmd), shell=True)
    """
    #not used because sometimes independent (popened) from salome
    import SalomePyQt
    sgPyQt = SalomePyQt.SalomePyQt()
    sgPyQt.helpContext(maDoc,"")
    """
  else:
    logger.error("inexisting nameUrlHelp: '%s'" % nameUrlHelp )
