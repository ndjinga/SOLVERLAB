#!/usr/bin/env python
#-*- coding:utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


"""
define class as a simple dictionary with keys
with pretty print __str__ and __repr__ (indented as recursive)
and jsonDumps()

| Usage:
| >> import catchAll as CAA
| >> a = CAA.CatchAll()
| >> a.tintin = "reporter"
| >> a.milou = "dog"
| >> print("a=%s" % a)
| >> print("tintin: %s" % a.tintin)
"""

import pprint as PP
import json

########################################################################################
# json utils
########################################################################################
def dumper(obj):
    """to json explore subclass object as dict"""
    return obj.__dict__

def dumperType(obj):
    """\
    to get a "_type" to trace json subclass object,
    but ignore all attributes begining with '_'
    """
    typeatt = "_type"
    aDict = dict((k,v) for k, v in obj.__dict__.items() if k[0] != "_" or k == typeatt)
    if typeatt not in aDict: aDict[typeatt] = obj.__class__.__name__
    return aDict

def jsonDumps(obj):
    """to get direct default jsonDumps method"""
    return json.dumps(obj, default=dumperType, sort_keys=True, indent=2)


########################################################################################
class CatchAll(object):
  """
  class as simple dynamic dictionary 
  with predefined keys as properties in
  inherited classes through __init__ method. Or NOT. 
  with pretty print __str__ and __repr__ (indented as recursive)
  with jsonDumps()
  
  | Usage:
  | >> import catchAll as CAA
  | >> a = CAA.CatchAll()
  | >> a.tintin = "reporter"
  | >> a.milou = "dog"
  | >> print("a=%s" % a)
  | >> print("tintin: %s" % a.tintin)
  | 
  | as
  
  | >> a = {}
  | >> a["tintin"] = "reporter"
  | >> a["milou"] = "dog"
  | >> print("tintin: %s" % a["tintin"]
  """
  
  def __repr__asList(self):
    """\
    goal is to be unambiguous
    an ordered list representation is better for test (and visualize) (in)equality
    """
    aList = []
    for k in sorted(self.__dict__.keys()):
      if k[0] != '_':
        aList.append( [k, self.__dict__[k]]  )
    return self.__class__.__name__ + " = " + aList.__repr__()

  def __repr__(self):
    """goal is to be unambiguous, easy human readeable"""
    return self._reprIndent()
  
  def _reprIndent(self, indent=0):
    res = ""
    newIndent = indent + 2
    for k in sorted(self.__dict__.keys()):
      if k[0] != '_':
        kk = self.__dict__[k]
        if issubclass(CatchAll, kk.__class__):
          res += "\n" + " "*newIndent + "%s: %s" % (k, kk._reprIndent(newIndent))
        else:
          skk = self._indent(PP.pformat(kk), newIndent)
          res += "\n" + " "*newIndent + "%s: %s" % (k, skk)
    return self.__class__.__name__ + "(" + res + ")"
  
  def _indent(self, txt, indent):
    txts = txt.split("\n")
    if len(txt) > 1:
      return ("\n" + " "*indent).join(txts)
    else:
      return txt

  def jsonDumps(self):
    return jsonDumps(self)

