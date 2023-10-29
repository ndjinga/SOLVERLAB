#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
http://sametmax.com/ecrire-des-logs-en-python/
https://docs.python.org/3/library/time.html#time.strftime

essai utilisation logger plusieurs handler format different

    /usr/lib/python2.7/logging/__init__.pyc

    init MyLogger, fmt='%(asctime)s :: %(levelname)-8s :: %(message)s', level='20'

    2018-03-11 18:51:21 :: INFO     :: test logger info
    2018-03-11 18:51:21 :: WARNING  :: test logger warning
    2018-03-11 18:51:21 :: ERROR    :: test logger error
    2018-03-11 18:51:21 :: CRITICAL :: test logger critical

    init MyLogger, fmt='None', level='10'

    2018-03-11 18:51:21 :: DEBUG    :: test logger debug
    test logger debug
    2018-03-11 18:51:21 :: INFO     :: test logger info
    test logger info
    2018-03-11 18:51:21 :: WARNING  :: test logger warning
    test logger warning
    2018-03-11 18:51:21 :: ERROR    :: test logger error
    test logger error
    2018-03-11 18:51:21 :: CRITICAL :: test logger critical
    test logger critical
"""

import os
import sys
import logging
import pprint as PP

print(logging.__file__)

def initMyLogger(fmt=None, level=None):
  # http://sametmax.com/ecrire-des-logs-en-python/
  # https://docs.python.org/3/library/time.html#time.strftime
  print("\ninit MyLogger, fmt='%s', level='%s'\n" % (fmt, level))
  logger = getMyLogger()
  handler = logging.StreamHandler() # Logging vers console
  if fmt is not None:
    formatter = logging.Formatter(fmt, "%Y-%m-%d %H:%M:%S")
    handler.setFormatter(formatter)
  logger.addHandler(handler)
  if level is not None:
    logger.setLevel(level)
  else:
    logger.setLevel(logger.INFO) # ou DEBUG
  # logger.info('\n' + PP.pformat(dir(logger)))
  
def getMyLogger():
  return logging.getLogger('MyLogger')

def testLogger1():
  logger = getMyLogger()
  logger.debug('test logger debug')
  logger.info('test logger info')
  logger.warning('test logger warning')
  logger.error('test logger error')
  logger.critical('test logger critical')
  
if __name__ == "__main__":
  initMyLogger('%(asctime)s :: %(levelname)-8s :: %(message)s', level=logging.INFO)
  testLogger1()
  initMyLogger(level=logging.DEBUG)
  testLogger1()
