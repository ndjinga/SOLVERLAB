#!/usr/bin/env python
#-*- coding:utf-8 -*-

import os
import sys
import time

fifoname = "./pipefifo_001"

mess = "hello from write_pipefifo_001.py"
if len(sys.argv) == 2:
  mess = mess + " <%s>" % str(sys.argv[1])

#test presence of pipefifo_001
if not os.path.exists(fifoname):
  os.mkfifo(fifoname)   # create a named pipe file
  print("create fifo %s" % fifoname)

# two types of open of named pipe fifo, 
# blocking on open 'r' other side

for i in range(10):
  if False: #True:
    pipeout = os.open(fifoname, os.O_WRONLY)
    os.write(pipeout, mess + " %03d" % i)
    os.close(pipeout)

  else:
    with open(fifoname, "w") as f:
      f.write(mess + " %03d" % i)

  time.sleep(3)

print("end %s" % __file__)
