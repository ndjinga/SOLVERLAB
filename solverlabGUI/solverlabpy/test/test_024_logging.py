#!/usr/bin/env python
#-*- coding:utf-8 -*-

#  Copyright (C) 2010-2023  CEA
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
test logger for solverlabGui
"""


"""
| see: http://sametmax.com/ecrire-des-logs-en-python
|
|  # creation d'un formateur qui va ajouter le temps, le niveau
|  # de chaque message quand on ecrira un message dans le log
|  formatter = logging.Formatter('%(asctime)s :: %(levelname)s :: %(message)s')
|
|  # creation d'un handler qui va rediriger une ecriture du log vers
|  # un fichier en mode 'append', avec 1 backup et une taille max de 1Mo
|  file_handler = RotatingFileHandler('activity.log', 'a', 1000000, 1)
|
|  # on lui met le niveau sur DEBUG, on lui dit qu'il doit utiliser le formateur
|  # cree precedement et on ajoute ce handler au logger
|  file_handler.setLevel(logging.DEBUG)
|  file_handler.setFormatter(formatter)
|  logger.addHandler(file_handler)
|   
|  # creation d'un second handler qui va rediriger chaque ecriture de log
|  # sur la console
|  stream_handler = logging.StreamHandler()
|  stream_handler.setLevel(logging.DEBUG)
|  logger.addHandler(stream_handler)
|
|  # Après 3 heures, on peut enfin logguer
|  # Il est temps de spammer votre code avec des logs partout :
|  logger.info('Hello')
|  logger.warning('Testing %s', 'foo')
"""

import os
import sys
import unittest
import pprint as PP
import logging as LOGI
from logging.handlers import BufferingHandler

import debogpy.debug as DBG

verbose = False #True

_TRACE = LOGI.INFO - 2 # just below info

class LoggerTest(LOGI.Logger):
  """
  Elementary prototype for logger
  add a level TRACE as log.trace(msg) 
  below log.info(msg)
  above log.debug(msg)
  to assume store long log asci in files txt under/outside files xml
  
  see: /usr/lib64/python2.7/logging/xxx__init__.py etc.
  """
  
  def __init__(self, name, level=LOGI.INFO):
    """
    Initialize the logger with a name and an optional level.
    """
    super(LoggerTest, self).__init__(name, level)
    LOGI.addLevelName(_TRACE, "TRACE")
    # LOGI.TRACE = _TRACE # only for coherency,
    
  def trace(self, msg, *args, **kwargs):
    """
    Log 'msg % args' with severity '_TRACE'.

    To pass exception information, use the keyword argument exc_info with
    a true value, e.g.

    logger.trace("Houston, we have a %s", "long trace to follow")
    """
    if self.isEnabledFor(_TRACE):
        self._log(_TRACE, msg, args, **kwargs)

class TestCase(unittest.TestCase):
  "Test the debug.py"""
  
  initialLoggerClass = [] # to keep clean module logging
  
  def test_000(self):
    # one shot setUp() for this TestCase
    self.initialLoggerClass.append(LOGI._loggerClass)
    LOGI.setLoggerClass(LoggerTest)
    if verbose:
      DBG.push_debug(True)
      # DBG.write("assert unittest", [a for a in dir(self) if "assert" in a])
    pass
  
  def test_999(self):
    # one shot tearDown() for this TestCase
    if verbose:
      DBG.pop_debug()
    LOGI.setLoggerClass(self.initialLoggerClass[0])
    return
  
  def test_010(self):
    # LOGI.setLoggerClass(LoggerTest) # done once in test_000
    name = "testLogging"
    lgr = LOGI.getLogger(name) # create it
    lgr.setLevel("DEBUG")
    self.assertEqual(lgr.__class__, LoggerTest)
    self.assertEqual(lgr.name, name)
    self.assertIn("trace", dir(lgr))
    self.assertEqual(LOGI.getLevelName(LOGI.INFO), "INFO")
    self.assertEqual(LOGI.getLevelName(_TRACE), "TRACE")
    
    # creation d'un handler pour chaque log sur la console
    formatter = LOGI.Formatter('%(levelname)-8s :: %(message)s')
    # stream_handler = LOGI.handlers.StreamHandler() # log outputs in console
    stream_handler = LOGI.handlers.BufferingHandler(1000) # log outputs in memory
    stream_handler.setLevel(LOGI.DEBUG)
    stream_handler.setFormatter(formatter)
    lgr.addHandler(stream_handler)
    # print # skip one line if outputs in console
    lgr.warning("!!! test warning")
    lgr.info("!!! test info")
    lgr.trace("!!! test trace")
    lgr.debug("!!! test debug")
    self.assertEqual(len(stream_handler.buffer), 4)
    rec = stream_handler.buffer[-1]
    self.assertEqual(rec.levelname, "DEBUG")
    self.assertEqual(rec.msg, "!!! test debug")
    self.assertEqual(stream_handler.get_name(), None) # what to serve ?
    
  def test_020(self):
    # LOGI.setLoggerClass(LoggerTest)
    name = "testLogging"
    lgr = LOGI.getLogger(name) #  find it as created yet in test_010
    stream_handler = lgr.handlers[0]
    rec = stream_handler.buffer[-1]
    self.assertEqual(rec.levelname, "DEBUG")
    self.assertEqual(rec.msg, "!!! test debug")

      
          
    
if __name__ == '__main__':
    unittest.main(exit=False)
    pass

