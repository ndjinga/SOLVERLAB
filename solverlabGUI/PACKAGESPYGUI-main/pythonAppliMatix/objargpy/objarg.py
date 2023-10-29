#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
make class to implement argparse

provides usage as argparse namespaces
 self is a current instance derived object, could be modified
 self.getInitialNamespace(), a copy of the initial namespace

memorize initial values allowing reset modified values in object
 self.resetCurrentNamespace()

provides str(self) for print simple current namespace
self.reprAll() for print recursive tree of current namespace

see objargTest.py
"""

import os
import sys
import argparse as AP
import pprint as PP #pretty print

verbose = False

#########################################
def exec2or3(aString):
  """python 2-3 exec() compatibility"""
  if verbose: print("exec %s"% aString)
  namespace = {}
  try:
    exec(aString, namespace)
  except:
    msg = "error in '%r'" % aString
    raise Exception(msg)
  return namespace

#########################################
class ObjArg(object):

  def __init__(self, args=[]):
    super(ObjArg, self).__init__()
    self._initialArgs = args
    self.setCurrentNamespace(args)

  def getHelp(self):
    return self.getParser().format_help()

  def printHelp(self):
    self.getParser().print_help()

  def getInitialNamespace(self):
    return self.getParser().parse_args(self._initialArgs) #a copy

  def setCurrentNamespace(self, args):
    res = self.getParser().parse_args(args)
    self._initAttributes(res)

  def resetCurrentNamespace(self):
    res = self.getParser().parse_args(self._initialArgs)
    self._initAttributes(res)

  def _initAttributes(self, aNamespace):
    for nameAttr, value in aNamespace._get_kwargs():
      setattr(self, nameAttr, value)

  def _get_kwargs(self):
    """
    as _get_kwargs from namespace argparse
    returns list of attributes, avoid name attributes bebining by '_'
    """
    return [(name, self.__dict__[name]) for name in sorted(self.__dict__.keys()) if name[0] != '_']

  def errorParser(self, msg):
    """
    avoid
    python2.7/argparse.py", line 2362, and exit _sys.exit(status)
    raise Exception instead
    """
    mymsg = "ObjArg parser error: %s" % msg
    if verbose: print(mymsg)
    raise Exception(mymsg)


  #################################################
  #to string method for print
  ################################################# 

  def __str__(self):
    """simple representation and not recursive, is only for first depth info, if tree"""
    return "%s(\n%s\n)" % (self.__class__.__name__, PP.pformat(self._get_kwargs()))

  def reprAll(self, indent=2):
    """
    recursive indented representation, if tree
    user have to avoid loop, no test
    """
    res = self.__class__.__name__ + "([\n"
    for name, val in self._get_kwargs():
      if issubclass(val.__class__, ObjArg):
        res += " "*indent + "('" + name + "', " +val.reprAll(indent+2) + ",\n"
      else: 
        res += " "*indent + str((name, val)) + ",\n"
    res += " "*(indent-2) + "])" 
    return res

  #################################################
  #parser filter methods
  #################################################

  def filter_list(self, string):
    """
    parser filter from string 'xx,yy,zz,...'
    returns list (if not error with python exec(value=[xx,yy,zz,...]))
    """
    res = exec2or3("value = [%s]" % string)
    return res["value"]

  def filter_list_float(self, string):
    """
    parser filter from string 'xx,yy,zz,...'
    returns list (if not error with python exec(value=[xx,yy,zz,...]))
    """
    value = self.filter_list(self, string)
    value = [float(v) for v in value]
    return value

  def filter_list_int(self, string):
    """
    parser filter from string 'xx,yy,zz,...'
    returns list (if not error with python exec(value=[xx,yy,zz,...]))
    """
    value = self.filter_list(self, string)
    value = [int(v) for v in value]
    return value

  def filter_range(self, string):
    """
    parser filter from string 'vmin,vmax'
    returns list (if not error with python exec(value=[xx,yy]))
    """
    try:
      value = self.filter_list_float(self, string)
      if len(value) != 2:
        msg = "%r is not a range of 2 float" % string
        raise Exception(msg)
    except:
      msg = "%r is not a list of 2 float" % string
      raise Exception(msg)
    return value

  def filter_range_float(self, string):
    """
    parser filter from string 'vmin,vmax'
    returns list (if not error with python exec(value=[xx,yy]))
    """
    return self.filter_range(self, string)

  def filter_range_int(self, string):
    """
    parser filter from string 'vmin,vmax'
    returns list (if not error with python exec(value=[xx,yy]))
    """
    try:
      value = self.filter_list_int(self, string)
      if len(value) != 2:
        msg = "%r is not a range of 2 int" % string
        raise Exception(msg)
    except:
      msg = "%r is not a list of 2 int" % string
      raise Exception(msg)
    return value

  def filter_square(self, string):
    value = int(string)
    sqrtv = value**(.5)
    if sqrtv != int(sqrtv):
      msg = "%r is not a perfect square" % string
      raise Exception(msg)
    return value

  def filter_int_positive(self, string):
    try:
      value = int(string)
      if value < 0:
        msg = "%r is not a positive integer number" % string
        raise Exception(msg)
    except:
      msg = "%r is not a correct positive integer number" % string
      raise Exception(msg)
    return value

  def filter_float_positive(self, string):
    try:
      value = float(string)
      if value < 0:
        msg = "%r is not a positive float number" % string
        raise Exception(msg)
    except:
      msg = "%r is not a correct positive float number" % string
      raise Exception(msg)
    return value

  def filter_existing_file(self, string):
    try:
      ok = os.path.isfile(string)
      if not ok:
        msg = "%r is not an existing file" % string
        raise Exception(msg)
    except:
      msg = "%r is not a correct existing file name" % string
      raise Exception(msg)
    return value

  #################################################
  #virtual parser for example
  ################################################# 

  def getParser(self):
    """
    as virtual example getParser of user needs inerited classes
    user can really create all filter_xxx method(s) to smart parsing
    arguments as '_toto' are forbidden, not keyworded parser arg, reserved for user needs 
    """
    parser = AP.ArgumentParser(
      prog='ObjArg',
      description='ObjArg use argparse to init/control initial scalar values of predeterminated attributes. This is virtual example implementation: you have to create your inherited class where you define your method getParser()',
      #usage="Test design of ObjArgObjArg",
      argument_default=None)
    parser.add_argument(
      '-v', '--verbose', 
      help='set verbose, for debug',
      action='store_true',
    )
    parser.add_argument(
      '-f', '--firstAttInt', 
      help='my first attribute integer',
      type=int,
      default=1,
      choices=[1, 5, 7]
    )
    parser.add_argument(
      '-s', '--secondAttFloat', 
      help='my second attribute float',
      type=float,
      default=2.2
    )
    parser.add_argument(
      '-t', '--thirdAttList', 
      help='my third attribute list comma separated',
      type=self.filter_list,
      default="123,456.7,8910"
    )
    parser.add_argument(
      '-ps', '--perfectSquare', 
      help='my attribute perfect square',
      type=self.filter_square,
      default=9
    )

    parser.error = self.errorParser
    return parser


