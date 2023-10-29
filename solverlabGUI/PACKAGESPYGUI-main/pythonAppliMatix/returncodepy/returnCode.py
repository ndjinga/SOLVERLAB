#!/usr/bin/env python
#-*- coding:utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


"""
This file contains ReturnCode class

| Usage:
| >> import returnCode as RCO
"""

import pprint as PP

# global module variable
_OK_STATUS = "OK"
_KO_STATUS = "KO"
_NA_STATUS = "NA" # not applicable
_UNKNOWN_STATUS = "ND" # not defined
_KNOWNFAILURE_STATUS = "KF"
_TIMEOUT_STATUS = "TIMEOUT"

#####################################################
class ReturnCode(object):
  """
  assume simple return code for methods, with explanation as 'why'.
  Obviously why is 'why it is not OK', 
  but also why is 'why it is OK' (if you want). 
  Optionaly contains a return value as self.getValue()
  
  | Usage:
  | >> import returnCode as RCO
  | 
  | >> aValue = doSomethingToReturn()
  | >> return RCO.ReturnCode("OK", "there is no problem here", aValue)
  | >> return RCO.ReturnCode("KO", "there is a problem here because etc", None)
  | >> return RCO.ReturnCode("TIMEOUT_STATUS", "too long here because etc")
  | >> return RCO.ReturnCode("NA", "not applicable here because etc")
  | 
  | >> rc = doSomething()
  | >> print("short returnCode string", str(rc))
  | >> print("long returnCode string with value", repr(rc))
  | 
  | >> rc1 = RCO.ReturnCode("OK", ...)
  | >> rc2 = RCO.ReturnCode("KO", ...)
  | >> rcFinal = rc1 + rc2
  | >> print("long returnCode string with value", repr(rcFinal)) # KO!
  | 
  | >> rc = doSomething()
  | >> if rc.isOk(): doSomethingAsOK()
  | >> if not rc.isOk(): doSomethingAsKO()
  | 
  | >> rc = doSomething().raiseIfKo() # raise Exception if KO
  | >> doSomethingWithValue(rc.getValue()) # here i am sure that is OK
  """

  # redunctant but useful class variables
  OK_STATUS = _OK_STATUS
  KO_STATUS = _KO_STATUS
  NA_STATUS = _NA_STATUS # not applicable
  UNKNOWN_STATUS = _UNKNOWN_STATUS # not defined
  KNOWNFAILURE_STATUS = _KNOWNFAILURE_STATUS
  TIMEOUT_STATUS = _TIMEOUT_STATUS

  # an integer for sys.exit(anInteger)
  # OKSYS and KOSYS seems equal on linux or windows
  OKSYS = 0  # OK 
  KOSYS = 1  # KO
  NASYS = 2  # KO not applicable return code
  NDSYS = 3  # KO not defined return code
  KFSYS = 4  # KO known failure return code
  TOSYS = 5  # KO time out
  
  _TOSYS = { 
    OK_STATUS: OKSYS,
    KO_STATUS: KOSYS,
    NA_STATUS: NASYS,
    UNKNOWN_STATUS: NDSYS,
    KNOWNFAILURE_STATUS: KFSYS,
    TIMEOUT_STATUS: TOSYS, 
  }
  _DEFAULT_WHY = "No given explanation"
  _DEFAULT_VALUE = None

  def __init__(self, status=None, why=None, value=None):
    self._why = self._DEFAULT_WHY 
    self._value = self._DEFAULT_VALUE
    if status is None:
      self._status = self.UNKNOWN_STATUS
    else:
      self.setStatus(status, why, value)
    
  def __repr__(self):
    """complete with value, 'ok, why, value' message"""
    res = '%s: %s --value: %s' % (self._status, self._why, PP.pformat(self._value))
    return res
  
  def __str__(self):
    """without value, only simple 'ok, why' message"""
    res = '%s: %s' % (self._status, self._why)
    return res

  def indent(self, text, amount=5, ch=' '):
    """indent multi lines message"""
    padding = amount * ch
    res = ''.join(padding + line for line in text.splitlines(True))
    return res[amount:]

  def __add__(self, rc2):
    """allows expression 'returnCode1 + returnCode2 + ...' """
    isOk = self.isOk() and rc2.isOk()
    newWhy = self._toList(self.getWhy()) + self._toList(rc2.getWhy())
    newValue = self._toList(self.getValue()) + self._toList(rc2.getValue())    
    if isOk: 
      return ReturnCode("OK", newWhy, newValue)
    else:
      return ReturnCode("KO", newWhy, newValue)
    
  def __radd__(self, other):
    # see http://www.marinamele.com/2014/04/modifying-add-method-of-python-class.html
    if other == 0:
      return self
    else:
      return self.__add__(other) 
    
  def _toList(self, strOrList):
    """internal use"""
    if type(strOrList) is not list: 
      return [strOrList]
    else:
      return strOrList

  def toSys(self):
    """return system return code as bash or bat"""
    try:
      return self._TOSYS[self._status]
    except:
      return self._TOSYS[self.NA_STATUS]

  def toXmlPassed(self):
    """return xml  return code as '0' (passed) or '1' (not passed)"""
    if self.isOk(): 
      return "0"
    else:
      return "1"
    
  def getWhy(self):
    """return why as str or list if sum or some ReturnCode"""
    return self._why
    
  def setWhy(self, why):
    self._why = why
    
  def getValue(self):
    return self._value
    
  def setValue(self, value):
    """choice as not deep copying if mutables value"""
    # TODO deepcopy maybe for value, not yet
    self._value = value
    
  def setStatus(self, status, why=None, value=None):
    if why is None: 
      aWhy = self._DEFAULT_WHY
    else:
      aWhy = why
      
    if status in list(self._TOSYS.keys()):
      self._status = status
      self._why = aWhy
    else:
      self._status = self.NA_STATUS
      self._why = "Error status '%s' for '%s'" % (status, aWhy)
      
    if value is not None:
      # TODO deepcopy maybe for value, not yet
      self._value = value
    else:
      self._value = self._DEFAULT_VALUE
      
  def getStatus(self):
    return self._status

  def isOk(self):
    """
    return True if ok.
    inexisting method isKo(), use more explicit/readability 'if not res.isOk()'
    """
    return (self._status == self.OK_STATUS)
  
  def raiseIfKo(self):
    """
    raise an exception with message why if not ok, else return self.
    This trick is to write usage
    
    | Usage:
    | >> rc = doSomething().raiseIfKo() # raise Exception if KO
    | >> doSomethingWithValue(rc.getValue()) # here i am sure that is OK
    """
    if self.isOk(): 
      return self
    else:
      raise Exception(self.getWhy())

def ReturnCodeFromList(aListOfReturnCodes):
  """
  Create ReturnCode from list of ReturnCode
  
  convenience over "+" operand
  """
  res = "OK"
  whyes = []
  for rc in aListOfReturnCodes:
    if not rc.isOk():
      res = "KO"
    whyes.append(str(rc))
  reswhy = "\n  ".join(whyes)
  return ReturnCode(res, "\n  " + reswhy)
    
    