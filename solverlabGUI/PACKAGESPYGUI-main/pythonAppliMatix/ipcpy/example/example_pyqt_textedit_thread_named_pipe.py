#!/usr/bin/env python
#-*- coding:utf-8 -*-

_help = r"""

see http://flylib.com/books/en/2.723.1.41/1/
https://stackoverflow.com/questions/37252756/simplest-way-for-pyqt-threading/37256736

file ./pipefifo
     ./pipefifo: fifo (named pipe)

to launch 2 independent processes with pipe communication:
  >> #first process write
  >> write_pipefifo_001.py "yiiii" &
  >> ls -alt  #created  prw-r--r--  1 wambeke lgls    0 30 oct.  14:30 pipefifo_001
  >> #second process read, GUI, and click 'launchParentFifo' before 30 sec
  >> example_pyqt_textedit_thread_named_pipe.py & 

"""


import os
import sys
import time
import inspect
import random
import select
import pprint as PP
import subprocess as SP
import threading as THRD 
from PyQt5 import QtGui, QtCore, QtWidgets

verbose = True
verboseEvent = False

fifoname = './pipefifo'                          # must open same name
fifoname001 = './pipefifo_001'                   # must open same name


#######################################################
#######################################################
def createLaunchDemo3():
  name = "LaunchDemo_3.bash"
  with open(name, 'w') as f: f.write(_launchDemo3)
  chmodax(name)

def chmodax(aFile):
  """chmod a+x path"""
  st = os.stat(aFile)
  os.chmod(aFile, st.st_mode | 0o111)


#######################################################
#######################################################
def child():
  print('begin child()')
  mess = 'Hello from child %03d zzz'
  pipeout = os.open(fifoname001, os.O_WRONLY)  # blocking open fifo pipe file
  print('Child os.open pipeout done')
  zzz = 10
  for i in range(zzz):                         # roll 0 to zzz
    time.sleep(1)
    toSend = mess % i
    #print('Child %d WILL send "%s" at %s' % (os.getpid(), toSend[:-1], time.strftime("%Hh %Mm %Ss")))
    os.write(pipeout, toSend)
    print('Child %d send "%s" at %s' % (os.getpid(), toSend, time.strftime("%Hh %Mm %Ss")))
  print('end chid()')
  return
         
def parent():
  print('begin parent()')
  pipein = open(fifoname, 'r')                 # open fifo as stdio object
  print('Parent open pipein done')
  while True:
    line = pipein.readline()            # blocks until data sent
    print('Parent %d got "%s" at %s' % (os.getpid(), line[:-1] , time.strftime("%Hh %Mm %Ss")))
    if line == "": break
  print('end parent()')


#######################################################
#######################################################
class ThreadTextedit(THRD.Thread):
  """Thread emit signal to append text in widget textedit"""
  noThread = [1]

  def __init__(self, textedit):
    super(ThreadTextedit, self).__init__()
    self.textedit = textedit
    self.myName = "ThreadTextedit_%03d" % self.noThread[0]
    self.noThread[0] += 1
    self.index = 1
    self.toStop = False

  def run(self):
    """Code during executing thread"""
    message = "append text %03d from thread %s"
    for i in range(20):
      if self.toStop: break
      self.textedit.appendTextSignal.emit(message % (i, self.myName))
      attente = 0.2
      attente += random.randint(1, 60) / 100
      time.sleep(attente)
    print('end %s.run()' % self.myName)

  def stop(self):
    self.toStop = True

  def __repr__(self):
    res = "ThreadTextedit(%s)" % self.myName
    return res

#######################################################
class ThreadTexteditFifo(THRD.Thread):
  """Thread emit signal to append text in widget textedit from fifo"""
  noThread = [1]

  def __init__(self, textedit):
    super(ThreadTexteditFifo, self).__init__()
    self.textedit = textedit
    self.myName = "ThreadTexteditFifo_%03d" % self.noThread[0]
    self.fifoname = fifoname+"_%03d" % self.noThread[0]
    self.noThread[0] += 1
    self.index = 1
    self.child = None            #gnome-terminal child not launched
    self.toStop = False
    self.nonBlock = True
    self.debug = True
    self.ii = 0
    if not os.path.exists(self.fifoname):
      os.mkfifo(self.fifoname)   # create a named pipe file
      print("create fifo %s" % self.fifoname)
      time.sleep(.5)

  def run(self):
    if self.nonBlock:
      self.run_nonblock_os()
    else:
      self.run_block()

  def run_block(self):
    """Code during executing thread"""
    print('begin %s.run_block()' % self.myName)
    pipein = open(self.fifoname, 'r')  #blocking if not open "w" at other pipe side
    print('Parent open pipein %s done' % self.fifoname)
 
    message = "append text %03d from thread %s: '%s'"
    i = 0
    while True:
      if self.toStop: break
      line = pipein.readline()   # blocks until data sent
      self.textedit.appendTextSignal.emit(message % (i, self.myName, line[:-1]))
      i += 1
      if line == "\n": break
    print('end %s.run_block()' % self.myName)

  def run_nonblock_os(self):
    """
    Code during executing thread
    From the Open Group page for the open() function 
    O_NONBLOCK
      When opening a FIFO with O_RDONLY or O_WRONLY set: 
      If O_NONBLOCK is set:
          An open() for reading only will return without delay. An open()
          for writing only will return an error if no process currently
          has the file open for reading.

      If O_NONBLOCK is clear:
          An open() for reading only will block the calling thread until a
          thread opens the file for writing. An open() for writing only
          will block the calling thread until a thread opens the file for
          reading.
    """
    import errno
    print('begin %s.run_nonblock()' % self.myName)

    pipein = os.open(self.fifoname, os.O_RDONLY | os.O_NONBLOCK) # open fifo pipe file as fd
    print('Parent open pipein NONBLOCK %s done %s' % (self.fifoname, pipein) )
 
    message = "append text %03d from thread %s: '%s'"
    i = 0
    
    while True:
      i += 1
      if self.toStop: 
       print("NONBLOCK os toStop %s" % self.toStop)
       break
      #line = os.read(pipein, 10) #pipein.readline()

      try:
        line = os.read(pipein, 1024)
      except OSError as exc:
        line = "OSError %s (errno.EAGAIN=%s)" % (exc.errno, errno.EAGAIN)
        print(line)
        #print dir(errno)
        pass

      if line == "":
        self.textedit.appendTextSignal.emit(message % (i, self.myName, "NOMESSAGE"))
      else:
        self.textedit.appendTextSignal.emit(message % (i, self.myName, line))
      if line == "QUIT": break #as toStop:
      time.sleep(1)      
    print('end %s.run_nonblock()' % self.myName)

  def run_nonblock_posix(self):
    """
    Code during executing thread
   
    A process can open a FIFO in nonblocking mode. 
    In this case, opening or read-only will succeed even if no-one 
    has opened on the write side yet, opening for write-only will fail
    with ENXIO (no such device or address) 
    unless the other end has already been opened.
    
    So the first solution is opening FIFO with O_NONBLOCK.
    In this case you can check errno: 
    if it is equal to ENXIO, then you can try opening FIFO later
    """
    import posix
    print('begin %s.run_nonblock_posix()' % self.myName)

    pipein = posix.open(self.fifoname, posix.O_RDONLY |  posix.O_NONBLOCK)
    print('Parent open posix pipein NONBLOCK %s done %s' % (self.fifoname, pipein) )
       
    message = "append text %03d from thread %s: '%s'"
    i = 0

    while True:
      if self.toStop: break
      i += 1
      print("NONBLOCK posix something here")
      line = pipein.readline()
      self.textedit.appendTextSignal.emit(message % (i, self.myName, line[:-1]))
      if line == "\n": break
      time.sleep(.5)      
    print('end %s.run_nonblock()' % self.myName)

  def runChild(self):
    if self.debug:
      self.runChild_debug_os()
    else:
      self.runChild_nodebug()

  def runChild_nodebug(self):
    if self.child != None:
      print("child launched yet %s" % self.myName)
      return
    #create bash script for gnome-terminal
    #launch bash script
    cmd = "echo 'runChild %s' ; example_pyqt_textedit_thread_named_pipe.py -child ; sleep 10" % self.myName
    proc = SP.Popen(cmd, shell=True)
    self.child = proc
    print('end %s.runChild() %s "%s"' % (self.myName, self.fifoname, cmd))

  def runChild_debug_os(self):
    mess = 'Hello from child_debug_os %03d xxx'
    pipeout = os.open(self.fifoname, os.O_WRONLY)     # open fifo pipe file
    print('child_debug open pipeout done')
    self.ii += 1
    toSend = mess % self.ii
    os.write(pipeout, toSend)
    os.close(pipeout)
    print('Child %d send "%s" at %s' % (os.getpid(), toSend, time.strftime("%Hh %Mm %Ss")))
    print('end child_debug()')

  def runChild_debug(self):
    mess = 'Hello from child_debug open %03d yyy'
    with open(self.fifoname, "w") as pipeout:    # open fifo pipe file
      print('child_debug open pipeout done')
      self.ii += 1
      toSend = mess % self.ii
      pipeout.write(toSend)
      print('Child %d send "%s" at %s' % (os.getpid(), toSend, time.strftime("%Hh %Mm %Ss")))
    print('end child_debug()')

  def stop(self):
    self.toStop = True
    #and remove file fifoname
    #pipeout = os.open(self.fifoname, os.O_WRONLY)     # open fifo pipe file as fd
    #os.write(pipeout, "")
    if os.path.exists(self.fifoname):
      print("NO remove fifo %s" % self.fifoname)
      #os.remove(self.fifoname)

  def __repr__(self):
    res = "ThreadTexteditFifo(%s, child=%s)" % (self.myName, self.child)
    return res


#######################################################
#######################################################
def printdir(obj):
  print("\n******** %s ***********" % obj.__repr__())
  for i in dir(obj): print(i)

def printSlot():
  if verboseEvent: print("Slot: %s" % inspect.stack()[1][3])



#######################################################
#######################################################
_launchChildFifo=r"""
#!/usr/bin/env bash

#this is demonstration using gnome-terminal and example_pyqt_textedit_thread_named_pipe.py

#create bash scripts for gnome-terminal --command [bash scripts]
cat > launchChildFifo.bash << EOL
#!/usr/bin/env bash
echo '***** begin launchChildFifo *****'
pwd
ls -alt
file ./pipefifo
echo 'example_pyqt_textedit_thread_named_pipe.py -child'
example_pyqt_textedit_thread_named_pipe.py -child
echo '***** end launchChildFifo.bash in 10 sec *****'
sleep 10

EOL
chmod a+x launchChildFifo.bash

#find . -name "launch*.bash" -exec cat -n {} \;
more l*.bash | cat

#launch demo
gnome-terminal --command launchChildFifo.bash

echo '***** end launchDemo *****'

"""

class Button(QtWidgets.QPushButton):

  def __init__(self, text, mainWin):
    super(Button, self).__init__(text)
    #if verbose: print("create Button %s" % text)
    self.mainWin = mainWin
    self.clicked.connect(self.slotButton)
    self.index = 0

  def slotButton(self):
    """
    select case ugly, but for test:
    avoid connect and create slot in QMainWindowEdit
    """
    name = self.text()
    mess = "Button %s" % name
    print(mess)
    self.index += 1

    if name == "quit":
      self.mainWin.close()
      return
    
    if name == "appendText":
      self.mainWin.centralEdit.append("append text %03d from %s" % (self.index, mess))
      return

    if name == "appendTextSignal":
      self.mainWin.centralEdit.appendTextSignal.emit("append text %03d from %s" % (self.index, mess))
      return

    if name == "launchThread":
      self.mainWin.centralEdit.runThread()
      return

    if name == "listOfThreadsAlive":
      res = PP.pformat(THRD.enumerate())
      self.mainWin.centralEdit.append("listOfThreadsAlive:\n%s" % res)
      return

    if name == "launchParentFifo":
      self.mainWin.centralEdit.runParentFifo()
      return

    if name == "launchChildFifo":
      self.mainWin.centralEdit.runChildFifo()
      return


#######################################################
#######################################################
class TextEdit(QtWidgets.QTextEdit):

  appendTextSignal = QtCore.pyqtSignal(str)
  
  def __init__(self, *args, **kwargs):
    super(TextEdit, self).__init__(*args)
    self.appendTextSignal.connect(self.slotAppend)
    self.threads = []

  def slotAppend(self, aText):
    #print("slotAppend '%s'" % aText)
    self.append(aText)
  
  def runThread(self):
    aThread = ThreadTextedit(self)
    self.threads.append(aThread)
    #print("runThread %s.start()" % aThread.myName)
    self.append("runThread %s.start()" % aThread.myName)
    aThread.start()
  
  def runParentFifo(self):
    aThread = ThreadTexteditFifo(self)
    self.threads.append(aThread)
    #print("runThreadFifo %s.start()" % aThread.myName)
    self.append("runParentFifo %s.start()" % aThread.myName)
    aThread.start()

  def closeAllThreads(self):
    for th in self.threads: #THRD.enumerate():
      print("close", th)
      th.stop()
  
  def runChildFifo(self):
    print("TODO runChildFifo as gnome-terminal")
    for th in self.threads: #THRD.enumerate():
      try:
        if th.child == None: #do only once on first child None
          th.runChild()
          return
      except:
        pass #no runChild() for this thread
    return 
     
    '''aThread = ThreadTexteditFifo(self)
    self.threads.append(aThread)
    #print("runThreadFifo %s.start()" % aThread.myName)
    self.append("runThreadFifo %s.start()" % aThread.myName)
    aThread.start()'''
  
 
class QMainWindowEdit(QtWidgets.QMainWindow):

  def __init__(self, *args, **kwargs):
    super(QMainWindowEdit, self).__init__(*args)
    self.setWindowModality(QtCore.Qt.NonModal)
    self.centralEdit = None
    self.docks = []
    self.buttons = {}
    self.__addCentral()
    self.__addDocks()
    self.resize(900, 500)

  #def close(self):
  #  if verbose: print("QMainWindowEdit.close")
  #  return super(QMainWindowEdit, self).close()

  def closeEvent(self, event):
    if verbose: print("QMainWindowEdit.closeEvent")
    self.centralEdit.closeAllThreads()
    time.sleep(3)   
    res = PP.pformat(THRD.enumerate())
    print("listOfThreadsAlive:\n%s" % res)
    return super(QMainWindowEdit, self).closeEvent(event)

  def event(self, event):
    if verboseEvent: print("QMainWindowEdit.event %s" % event)
    return super(QMainWindowEdit, self).event(event)

  def __addCentral(self):
    central = TextEdit() #QtWidgets.QTextEdit()
    self.setCentralWidget(central)
    self.centralEdit = central

  def __addDocks(self):
    self.docks = []
    dock = QtWidgets.QDockWidget("DockForButtons", self)
    self.leftWidget = QtWidgets.QWidget()
    layout = QtWidgets.QVBoxLayout()

    buttons = """createPipeFifo
deletePipeFifo
appendText
appendTextSignal
launchThread
launchParentFifo
launchChildFifo
listOfThreadsAlive
quit""".split()
    for b in buttons: 
      bt = Button(b, self)
      self.buttons[b] = bt
      layout.addWidget(bt)
    self.leftWidget.setLayout(layout)
    dock.setWidget(self.leftWidget)

    #slots done in Button
    #self.buttons[bxx].clicked.connect(self.sbxx) #slots done in Button

    pos = QtCore.Qt.LeftDockWidgetArea
    dock.setAllowedAreas(pos)
    self.docks.append(dock)
    for dock in self.docks:
      self.addDockWidget(pos, dock)

  def sbxx(self):
    printSlot()
    #TODO something
    return


#######################################################
#######################################################
def run_2():
  if sys.argv[1] == "-parent":
    parent()                               # run as parent process
  elif sys.argv[1] == "-child":
    child()                                # run as child process
  else:
    print(_help)
  print('end of run_2')

def run_1():
  app = QtWidgets.QApplication(sys.argv)
  fen = QMainWindowEdit()
  fen.show()
  app.exec_()
  #avoid Erreur de segmentation (core dumped)
  del(fen)
  del(app)
  print('end of GUI run_1')



if __name__ == '__main__':
  if len(sys.argv) == 1:
    run_1()
    sys.exit()
  if len(sys.argv) == 2:
    run_2()
    sys.exit()
  print(_help)
  sys.exit()



