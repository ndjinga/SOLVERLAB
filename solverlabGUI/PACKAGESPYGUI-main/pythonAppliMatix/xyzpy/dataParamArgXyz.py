#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END

import string
import time
import os
import sys
import subprocess
import json
import pprint #pretty print

from xyzpy.baseXyz import BaseFreeXyz, BaseFreeAllXyz

import xyzpy.loggingXyz as LOG
logger = LOG.getLogger()

verbose = False

########################################################################################
# json utils
########################################################################################
def dumper(obj):
    """goal is to json explore subclass object as dict"""
    return obj.__dict__

def dumperType(obj):
    """goal is to get a "_type" trace json subclass object, but ignore attributes begining with '_'"""
    typeatt = "_type"
    if hasattr(obj, "_helpDict"):  #get help tooltip
      aDict = {}
      for k in list(obj.__dict__.keys()):
        if k[0] != "_" or k == typeatt:
          v = obj.__dict__[k]
          try:
            h = obj._helpDict[k][0]
          except:
            h = ""
          aDict[k] = v
          aDict["%s__help" % k] = h
          aDict["%s__treepath" % k] = v.getTreePyName()
    else:
      aDict = dict((k,v) for k, v in list(obj.__dict__.items()) if k[0] != "_" or k == typeatt)
    
    if typeatt not in aDict: aDict[typeatt] = obj.__class__.__name__
    return aDict

def jsonDumps(obj):
    """to get direct default jsonDumps method"""
    return json.dumps(obj, default=dumperType, sort_keys=True, indent=2)


########################################################################################
class DataParamArgXyz(object):
  """
  base class for user inherited another one:
  
  - manages input data life and correspondences
    between data model format Xyz
    and its file.xml representation, and other formats
  - formats are conf.py or conf.json and ?TODO? line argparse
  - input data are parameters of an user inherited class 
    whih have fonctionalities of a class 'engine'
  - give __repr__(self) human-eyes readable.
  - TODO give saveState(self) human-eyes readable.
  - TODO give loadState(self)
  """
  def __init__(self):
    self._model_initial = None #set only once, no modifications#to be modified
    self.internals = None #to be modified
    self._indent = 4
    
  def setModelInitial( self, aModel ):
    if self._model_initial !=  None:
      msg = "model initial set yet."
      logger.error(msg)
      raise Exception(msg)
    if not issubclass(aModel.__class__, BaseFreeXyz):
      msg = "model initial class '%s' have to be subclass of 'BaseFreeXyz'" % aModel.__class__.__name__
      logger.error(msg)
      raise Exception(msg)
    """
    here we store all on-choice user initial and internal engine attributes, (explicitly programming)
    saved/restored through methods saveState(self), loadState(self)
    it could be done by pickle, but unreadable..., so xml ot json possibilities are provided
    """
    self._model_initial = aModel.duplicate() #no modifications
    self.internals = BaseFreeXyz()
    self.internals.model_current = aModel.duplicate() #with modifications... later in instance life
    self.internals.model_addings = self.createModelAddings() #an BaseFreeXyz instance, for default
    if verbose:
      logger.warning("model initial:\n%s" % self._model_initial.jsonDumps())
      #logger.warning("model initial: \n%s" % self.jsonDumps(self._model_initial)) #same at this time
      logger.warning("model_addings:\n%s" % self.jsonDumps(self.internals.model_addings))
      #logger.warning("model current %s" % self._model_current.jsonDumps())
      logger.warning("DataParamArgXyz with info help and type from model:\n%s" % self)
      #logger.warning("DataParamArgXyz with info help and type from model and pretty print:\n%s" % self.pformat())

  def jsonDumps(self, obj=None):
    if obj != None: 
      return jsonDumps(obj)
    return jsonDumps(self)
    
  def jsonloads(self, aStr):
    return json.loads(aStr)
  
  def jsonDumpsFile(self, obj, aFile):
    res = jsonDumps(obj)
    with open(aFile, "w") as f: f.write( res )
    return

  def jsonloadsFile(self, aFile):
    with open(aFile, "w") as f: aStr = f.read()
    return json.loads(aStr)

  def createModelAddings(self):
    """user have define his method createModelAddings as his Xyz class, or not on choice"""
    res = BaseFreeAllXyz() #just permissive... risky
    """
    if verbose: #only for example True:
      #user can do all that on BaseFreeAllXyz class, get xml json methods etc...
      res.essaiParam1 = 123
      res.essaiParam2 = [1,2,3]
      res.essaiParam3 = "123"
      res.essaiParam4 = {123: 456}
      res.essaiParam5 = None
      res.essaiParam6 = False
      #print "ressss",res
    """
    logger.warning("DataParamArgXyz.createModelAddings virtual as 'BaseFreeAllXyz', permissive")
    return res
    
  def __repr__(self):
    """focus self.internals, are only interesting"""
    return jsonDumps(self.internals)

  def pprint(self):
    """pretty print with json"""
    data = json.loads(jsonDumps(self.internals))
    pprint.pprint(data, indent=self._indent)
    
  def pformat(self):
    """pretty format with json"""
    data = json.loads(jsonDumps(self.internals))
    return pprint.pformat(data, indent=self._indent)

  def toStrFileConf(self):
    """TODO begin give saveState with json"""
    data = json.loads(jsonDumps(self.internals))
    print("toStrFileConf TODO with data:\n",pprint.pformat(data, indent=self._indent))
    return
    
