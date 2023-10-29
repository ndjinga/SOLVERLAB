#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
http://sametmax.com/ecrire-des-logs-en-python/
https://docs.python.org/3/library/time.html#time.strftime

essai utilisation logger un handler format different 
sur info() pas de format et su other format

    /usr/lib/python2.7/logging/__init__.pyc

    init MyLogger, fmt='%(asctime)s :: %(levelname)-8s :: %(message)s', level='20'

    test logger info
    2018-03-11 18:51:51 :: WARNING  :: test logger warning
    2018-03-11 18:51:51 :: ERROR    :: test logger error
    2018-03-11 18:51:51 :: CRITICAL :: test logger critical

"""

import os
import sys
import logging
import pprint as PP

print(logging.__file__)

class MyFormatter(logging.Formatter):
  def format(self, record):
    # print "", record.levelname #type(record), dir(record)
    if record.levelname == "INFO": 
      return str(record.msg)
    else:
      return super(MyFormatter, self).format(record)

def initMyLogger(fmt=None, level=None):
  # http://sametmax.com/ecrire-des-logs-en-python/
  # https://docs.python.org/3/library/time.html#time.strftime
  print("\ninit MyLogger, fmt='%s', level='%s'\n" % (fmt, level))
  logger = getMyLogger()
  handler = logging.StreamHandler() # Logging vers console
  if fmt is not None:
    # formatter = logging.Formatter(fmt, "%Y-%m-%d %H:%M:%S")
    formatter =MyFormatter(fmt, "%Y-%m-%d %H:%M:%S")
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

