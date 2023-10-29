#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import unittest
import sys
from salomepy.onceQApplication import OnceQApplication
import xyzpy.actionsFactoryXyz as ACFX
import xyzpy.loggingXyz as LOG

logger = LOG.getLogger()
verbose = False
otherSlotArgs = []

def myOtherSlot(args):
  if verbose: print("myOtherSlot",args)
  otherSlotArgs.append(args)
  return

class TestCase(unittest.TestCase):
  def test_000(self):
    LOG.pushLevel("CRITICAL")
    app = OnceQApplication()
    self.assertNotEqual(app, None)
    a1 = ACFX.QActionXyz()
    self.assertEqual(ACFX.addInCommonActions(a1, verbose), False)
    a1.name = "essai"
    self.assertEqual(ACFX.addInCommonActions(a1), True)
    self.assertEqual(ACFX.getCommonActionByName("essai"), a1)
    self.assertEqual(ACFX.getCommonActionByName("truc"), None)
    self.assertEqual(ACFX.addInCommonActions(a1), False) #done yet
    LOG.popLevel()

  def test_010(self):
    LOG.pushLevel("CRITICAL")
    a1 = ACFX.QActionXyz(name='essai010')
    self.assertEqual(ACFX.addInCommonActions(a1), True)
    res = a1.setAction(slot=a1.defaultSlot, signal=a1.defaultSignal, shortcut=None, tooltip=None, icon=None)
    self.assertEqual(res, True)
    self.assertEqual(a1.lastDefaultSlotArgs, [])
    anArg = {"kwarg1": "hello"} #like **kwargs , for useful, or something else..;
    a1.signal.emit(anArg)
    """
    #curious: http://qt-project.org/doc/qt-5.1/qtwidgets/qaction.html#trigger
    a1.setChecked(False)
    a1.trigger()
    a1.triggered.emit(True)
    a1.triggered.emit(False)
    a1.triggered.emit(True)
    a1.trigger()
    a1.setChecked(True)
    a1.trigger()
    """
    self.assertEqual(a1.lastDefaultSlotArgs, [anArg])

    res = a1.setAction(slot=a1.defaultSlot, signal=a1.defaultSignal, shortcut=None, tooltip=None, icon=None, verbose=verbose)
    self.assertEqual(res, False)
    a1.signal.emit(anArg)
    self.assertEqual(a1.lastDefaultSlotArgs, [anArg, anArg])

    a2 = ACFX.QActionXyz(name='essai010bis')
    #explicitely stupid affect again same slot to same signal
    #TODO could be protected if tests... __future__
    res = a2.setAction(slot=a1.defaultSlot, signal=a1.defaultSignal, shortcut=None, tooltip=None, icon=None)
    self.assertEqual(res, True)
    a1.lastDefaultSlotArgs = []
    self.assertEqual(a1.lastDefaultSlotArgs, [])
    anArg2 = {"kwarg1": "hellobis"}
    a2.signal.emit(anArg2) #only one emit for two calls slot
    #stupid: slot call twice because connect twice
    self.assertEqual(a1.lastDefaultSlotArgs, [anArg2, anArg2])
    LOG.popLevel()

  def test_011(self):
    LOG.pushLevel("CRITICAL")
    a1 = ACFX.getCommonActionByName('essai010')
    self.assertNotEqual(a1, None)
    a2 = ACFX.getCommonActionByName('essai010bis')
    self.assertEqual(a2, None)
    a1.lastDefaultSlotArgs = []
    anArg3 = {"kwarg1": "helloter"}
    a1.signal.emit(anArg3)
    #stupid continue: slot call twice even a garbage collecting?
    self.assertEqual(a1.lastDefaultSlotArgs, [anArg3, anArg3])
    LOG.popLevel()

  def test_100(self):
    LOG.pushLevel("CRITICAL")
    a1 = ACFX.QActionXyz(name='essai0100')
    res = a1.setAction(slot=a1.defaultSlot, signal=a1.defaultSignal, shortcut=None, tooltip=None, icon=None)
    a1.signal.emit({"myname": a1.name})
    self.assertEqual(a1.lastDefaultSlotArgs, [{'myname': 'essai0100'}])

    a2 = a1.createNewAction('essai100bis', myOtherSlot, verbose=verbose)
    self.assertEqual(a2.name, 'essai100bis')
    self.assertEqual(a1.lastDefaultSlotArgs,  [{'myname': 'essai0100'}])
    self.assertEqual(a2.lastDefaultSlotArgs, [])

    self.assertEqual(otherSlotArgs, [])
    a2.trigger()
    self.assertEqual(otherSlotArgs, [False])

    a2.triggered.emit(True)
    self.assertEqual(otherSlotArgs, [False, True])
    a2.signal.emit({"myname": a2.name})
    #here we know write a local def mySlot ending by a2.signal.emit(args)
    #to reach a1.slot with myArgs
    self.assertEqual(a1.lastDefaultSlotArgs, [{'myname': 'essai0100'}, {'myname': 'essai100bis'}])
    LOG.popLevel()

if __name__ == '__main__':
  verbose = False #True
  unittest.main()
  pass
