#!/usr/bin/env python
#-*- coding:utf-8 -*-

r"""
To distinguish messages better, we can mandate a separator character in the pipe. 
An end-of-line makes this easy, because we can wrap the pipe descriptor 
in a file object with os.fdopen, and rely on the file object's readline method 
to scan up through the next \n separator in the pipe. 
Example implements this scheme. 

wrap pipe input in stdio file object
to read by line, and close unused pipe fds in both processes

to launch:
  >> essai_names_pipe_2.py

see http://flylib.com/books/en/2.723.1.41/1/
"""

import os, time

def child(pipeout):
    print('begin child()')
    mess = 'Hello from child %03d\n'
    zzz = 7
    for i in range(zzz):                         # roll 0 to zzz
        time.sleep(1)                            # make parent wait
        os.write(pipeout, mess % i)              # send to parent
    print('end chid()')
    return
         
def parent():
    print('begin parent()')
    pipein, pipeout = os.pipe()                  # make 2-ended pipe
    if os.fork() == 0:                           # in child, write to pipe
        os.close(pipein)                         # close input side here
        child(pipeout)
        print('end of child in parent')
    else:                                        # in parent, listen to pipe
        os.close(pipeout)                        # close output side here
        pipein = os.fdopen(pipein)               # make stdio input object
        while True:
            line = pipein.readline()[:-1]        # blocks until data sent
            print('Parent %d got "%s" at %s' % (os.getpid(), line, time.strftime("%Hh %Mm %Ss")))
            if line == "": break
        print('end parent()')

parent()
print('end of all')
