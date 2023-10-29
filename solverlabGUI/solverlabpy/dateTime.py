#!/usr/bin/env python
#-*- coding:utf-8 -*-

#  Copyright (C) 2018-20xx  CEA/DEN
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA

"""
This file contains DateTime and DeltaTime class

| Usage:
| >> import dateTime as DATT
| >> ini = DATT.DateTime("now")
| >> # some stuff
| >> fin = DATT.DateTime("now")
| >> duration = DATT.DeltaTime(ini, fin)
"""

import datetime as DT
import time as TI

# global module variable
verbose = False

#####################################################
class DateTime(object):
  """
  assume storing a date and hour, and conversions
  
  | Usage:
  | >> import dateTime as DATT
  | >> now = DATT.DateTime("now")
  | >> print("now is %s" % now)
  """

  FORMAT_HUMAN = '%Y-%m-%d %H:%M:%S' # human readable
  FORMAT_FILE = '%Y%m%d_%H%M%S' # for file name
  
  FORMAT_HOUR_CONFIG = '%H%M%S' # for config pyconf
  FORMAT_DATE_CONFIG = '%Y%m%d' # for config pyconf
  FORMAT_DATEHOUR_CONFIG = '%Y%m%d_%H%M%S' # for config as FORMAT_FILE
  
  FORMAT_PACKAGE = '%Y-%m-%d %H:%M' # for another
  FORMAT_XML = '%Y/%m/%d %Hh%Mm%Ss' # for log file xml
  
  MSG_UNDEFINED = "UndefinedTime"

  def __init__(self, when=None):
    self._time = None # set "UndefinedTime", else is a float
    if verbose: print("when", when)
    if type(when) == str:
      if when == "now":
        self._time = TI.time() # is a float
      else:
        raise Exception("DateTime: unknown when '%s'" % when)
    elif type(when) == self.__class__:
      self._time = when.getValue()
    elif type(when) == DT.datetime:
      # convert from datetime to time
      self._time = TI.mktime(when.timetuple())
    elif (type(when) == float) and (when > 1e9): # 1526469510 is may 2018
      self._time = when
    else:
      # UndefinedTime
      if verbose:# for debug 
        msg = "DateTime: unknown when %s '%s' implies 'UndefinedTime'" % (type(when), when)
        #raise Exception(msg)
        print(msg)
      pass
    
  def __add__(self, seconds):
    """add seconds"""
    return DateTime(self._time + seconds)
    
  def __repr__(self):
    """complete with type class as 'DateTime(2018-05-07 12:30:55)'"""
    res = "DateTime(%s)" % self
    return res
  
  def __str__(self):
    """human readable, sortable, as '2018-05-07 12:30:55'"""
    if self.isOk():
      res = TI.strftime(self.FORMAT_HUMAN, self.localTime())
    else:
      res = self.MSG_UNDEFINED
    return res
  
  def __eq__(self, other):
    return self._time == other._time

  def __gt__(self, other):
    return self._time > other._time

  def __ge__(self, other):
    return self._time >= other._time
  
  def addSeconds(self, secs):
    """add seconds at time"""
    self.raiseIfKo()
    self._time += secs
  
  def localTime(self):
    if self.isOk():
      return TI.localtime(self._time)
    else:
      return None
    
  def toStrFile(self):
    """use self.FORMAT_FILE, sortable, 2018-05-07... as '20180507_235958'"""
    if self.isOk():
      res = TI.strftime(self.FORMAT_FILE, self.localTime())
    else:
      res = self.MSG_UNDEFINED
    return res

  def toStrHuman(self):
    """use self.FORMAT_HUMAN"""
    if self.isOk():
      res = TI.strftime(self.FORMAT_HUMAN, self.localTime())
    else:
      res = self.MSG_UNDEFINED
    return res
  
  def toStrHourConfig(self):
    if self.isOk():
      res = TI.strftime(self.FORMAT_HOUR_CONFIG, self.localTime())
    else:
      res = self.MSG_UNDEFINED
    return res
  
  def toStrDateConfig(self):
    if self.isOk():
      res = TI.strftime(self.FORMAT_DATE_CONFIG, self.localTime())
    else:
      res = self.MSG_UNDEFINED
    return res
  
  def toStrDateHourConfig(self):
    if self.isOk():
      res = TI.strftime(self.FORMAT_DATEHOUR_CONFIG, self.localTime())
    else:
      res = self.MSG_UNDEFINED
    return res
    
  def toStrPackage(self):
    if self.isOk():
      res = TI.strftime(self.FORMAT_PACKAGE, self.localTime())
    else:
      res = self.MSG_UNDEFINED
    return res
    
  def toStrXml(self):
    if self.isOk():
      res = TI.strftime(self.FORMAT_XML, self.localTime())
    else:
      res = self.MSG_UNDEFINED
    return res
    
  def getValue(self):
    return self._time
    
  def toSeconds(self):
    return self._time
    
  def setValue(self, time):
    """choice as not deep copying if mutables value"""
    # TODO deepcopy maybe for value, not yet
    self._time = time
    
  def getSecondsToNow(self):
    delta = TI.time() - self._time
    return delta

  def isOk(self):
    """return True if ok"""
    return self._time is not None
  
  def raiseIfKo(self):
    """
    raise an exception with message why if not ok, else return self.
    This trick is to write usage
    
    | Usage:
    | >> aTimeOk = aTime.raiseIfKo() # raise Exception if KO
    | >> doSomethingWithaTimeOk(aTimeOk) # here i am sure that is OK
    """
    if self.isOk(): 
      return self
    else:
      raise Exception("DateTime not initialized")

#####################################################
class DeltaTime(object):
  """
  assume storing a duration, delta between two DateTime, and conversions
  
  | Usage:
  | >> import dateTime as DATT
  | >> t1 = DATT.DateTime("now")
  | >> time.sleep(3)
  | >> t2 = DATT.DateTime("now")
  | >> delta = DATT.DeltaTime(t1, t2)
  | >> print("delta time is %s" % delta)
  """

  MSG_UNDEFINED = "UndefinedDeltaTime"

  def __init__(self, t1=None, t2=None):
    try: 
      self._t1 = DateTime(t1)
    except:
      self._t1 = DateTime()
    try:
      self._t2 = DateTime(t2)
    except:
      self._t2 = DateTime()
    if type(t1) == str:
      if t1 == "now":
        self._t1 = DateTime(t1)
      else:
        raise Exception("DeltaTime: unknown t1 '%s'" % t1)
    
    if type(t2) == str:
      if t1 == "now":
        self._t2 = DateTime(t2)
      else:
        raise Exception("DeltaTime: unknown t2 '%s'" % t1)
        
  def __repr__(self):
    """complete with type class as 'DeltaTime(345.67)'"""
    res = "DeltaTime(t1=%s, t2=%s)" % (self._t1, self._t2)
    return res
  
  def __str__(self):
    """human readable, seconds, sortable, as '345.67'"""
    if self.isOk():
      res = "%s" % self.toSeconds()
    else:
      res = self.MSG_UNDEFINED
    return res
  
  def toSeconds(self):
    if self.isOk():
      res = self._t2.toSeconds() - self._t1.toSeconds()
    else:
      res = self.MSG_UNDEFINED
    return res
  
  def toMinutes(self):
    if self.isOk():
      res = (self._t2.toSeconds() - self._t1.toSeconds()) / 60
    else:
      res = self.MSG_UNDEFINED
    return res
  
  def toStrHuman(self):
    """automatic best unity, hours or minutes or seconds"""
    if self.isOk():
      res = self._t2.toSeconds() - self._t1.toSeconds()
      if res < 0: 
        sign = "-"
        res = abs(res)
      else:
        sign = ""
      if res < 0: return sign + "%.3fs" % res
      if res < 10: return sign + "%.3fs" % res
      if res < 60: return sign + "%is" % int(res)
      if res < 3600: return sign + "%im%is" % (int(res/60), int(res%60))
      return self.toStrHms()
    else:
      res = self.MSG_UNDEFINED
    return res
  
  def toStrHms(self):
    """all unities, hours and minutes and seconds as '2h34m56s'"""
    if self.isOk():
      res = self._t2.toSeconds() - self._t1.toSeconds()
      if res < 0: 
        sign = "-"
        res = abs(res)
      else:
        sign = ""
      hh = int(res/3600)
      mm = int(res%3600)/60
      ss = int(res%60)
      return sign + "%ih%im%is" % (hh, mm, ss)
    else:
      res = self.MSG_UNDEFINED
    return res
  
  def setT1(self, t):
    self._t1 = DateTime(t)
  
  def setT2(self, t):
    self._t2 = DateTime(t)
  
  def getT1(self, t):
    return DateTime(self._t1)
  
  def getT2(self, t):
    return DateTime(self._t2)
  
  def getValue(self):
    """idem toSeconds()"""
    return self.toSeconds()
  
  def isOk(self):
    """return True if ok"""
    return self._t1.isOk() and self._t2.isOk()
  
  def raiseIfKo(self):
    """
    raise an exception with message why if not ok, else return self.
    This trick is to write usage
    
    | Usage:
    | >> aDeltaTimeOk = adeltaTime.raiseIfKo() # raise Exception if KO
    | >> doSomethingWithaDeltaTimeOk(aDeltaTimeOk) # here i am sure that is OK
    """
    if not self._t1.isOk() and not self._t2.isOk(): 
      raise Exception("DeltaTime t1 and t2 are not initialized")
    if not self._t1.isOk():
      raise Exception("DeltaTime t1 not initialized")
    if not self._t2.isOk():
      raise Exception("DeltaTime t2 not initialized")
    return self # is ok
 
    
##############################################################################
# date utilities
##############################################################################
def sleep(seconds):
    """as time.sleep(seconds)"""
    TI.sleep(seconds)
  
def getWeekDayNow():
    """Returns monday as 0, tuesday as 1 etc."""
    return DT.date.weekday(DT.date.today())

def fromTimeStamp(val):
    """Returns datetime.datetime"""
    return DT.datetime.fromtimestamp(val)
    
def fromDateHourConfig(datehour):
    """
    datehour as pyconf config.VARS.datehour 'YYYYMMDD_HHMMSS'.
    Returns datetime.datetime
    """
    Y, m, dd, H, M, S = date_to_datetime(datehour)
    t0 = DT.datetime(int(Y), int(m), int(dd), int(H), int(M), int(S))
    return t0
    
def parse_date(date):
    """
    Transform as pyconf config.VARS.datehour 'YYYYMMDD_HHMMSS'
    to 'YYYY-MM-DD hh:mm:ss'.
    
    :param date: (str) The date to transform
    :return: (str) The date in the new format
    """
    if len(date) != 15:
        return date
    res = "%s-%s-%s %s:%s:%s" % (date[0:4],
                                 date[4:6],
                                 date[6:8],
                                 date[9:11],
                                 date[11:13],
                                 date[13:15])
    return res

def date_to_datetime(date):
    """
    From a string date as pyconf config.VARS.datehour 'YYYYMMDD_HHMMSS'
    returns [year, month, day, hour, minutes, seconds]
    
    :param date: (str) The date in format YYYYMMDD_HHMMSS
    :return: (tuple) as (str,str,str,str,str,str)
      The same date and time in separate variables.
    """
    Y = date[:4]
    m = date[4:6]
    dd = date[6:8]
    H = date[9:11]
    M = date[11:13]
    S = date[13:15]
    # print "date_to_datetime", date, [Y, m, dd, H, M, S]
    return Y, m, dd, H, M, S

def timedelta_total_seconds(timedelta):
    """
    Replace total_seconds from datetime module 
    in order to be compatible with old python versions
    
    :param timedelta: (datetime.timedelta) 
      The delta between two dates
    :return: (float) 
      The number of seconds corresponding to timedelta.
    """
    return (
        timedelta.microseconds + 0.0 +
        (timedelta.seconds + timedelta.days * 24 * 3600) * 10 ** 6) / 10 ** 6

