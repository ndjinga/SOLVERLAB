#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END

__all__=[]

import xyzpy.utilsXyz as UXYZ
import xml.etree.ElementTree as ET

verbose = False

########################################################################################
class ModelSimpleTest(object):
  """
  a simple model, ONLY for example test.
  inherits object, without nothing else to answer request
  user HAVE TO USE BaseFreeXyz _XyzConstrainBase
  """
  index = [0] # unambigous objectName

  def __init__(self, nameObject=""):
    self._lastReceiveRequest = None # firstly only for unittest
    if nameObject != "":
      self._objectName = nameObject
    else:
      self._objectName = "ModelSimpleTest[%s]" % self.index
      self.index[0] += 1              # unambigous objectName
    self._controller = None

  def setController(self, controller):
    """really could be None if no use in view without MVC pattern"""
    if self._controller == None:
      self._controller = controller
      return
    raise Exception("ModelSimpleTest.setController done yet for %s as %s" % (self.objectName(), self.getController().objectName()))

  def getController(self):
    return self._controller

  def objectName(self):
    """as Qt objectName"""
    return self._objectName

  def toXml(self, **kwargs):
    """kwarg are for optional future option of added details in xml tree"""
    #return UXYZ.toXml(self, **kwargs)
    res = ET.Element("ModelSimpleTest")
    res.text = "...TODO..."
    return res

  def lastReceiveRequest(self):
    """only for unittest"""
    return self._lastReceiveRequest

  def receiveRequest(self, aRequest):
    """
    asynchronous treatment of a request through a controller owner of a
    requestToModelSignal (QtCore.pyqtSignal)

    in controller
    - create signal and connect it
      >> aController.requestToModelSignal.connect(aModel.receiveRequest)

    in controller
    - have sendRequestToModel method where
      >> aController.requestToModelSignal.emit(aRequest)

    in view (for example)
    - have getController, AND DO NOT directly modify the model
      >> aController = aView.getController()
      >> aRequest = controller.getRequest() # for example
      >> # now define what request details exactly (modify model or else)
      >> aController.sendRequestToModel(aRequest)

    obviously request type typed yet
    it is simple test
    """
    if verbose:
      print("ModelSimpleTest %s receiveRequest:\n%s" % (self.objectName(), aRequest))
    self._lastReceiveRequest = aRequest # simple only for unittest
    res = True
    aRequest = UXYZ.fromXml(aRequest) # recreates a request through xml
    if "TestRequestFromView" in str(aRequest): # return the request to views for test
      res = self.getController().sendRequestToViews(aRequest)
      if not res: print("ModelSimpleTest.receiveRequest: Problem with request TestRequestFromView")
    return res
