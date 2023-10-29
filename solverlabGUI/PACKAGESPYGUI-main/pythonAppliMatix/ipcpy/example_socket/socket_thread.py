#!/usr/bin/env python
# coding: utf-8 

"""
this is demonstration using gnome-terminal and sockets

#launch demo execute:
  
  gnome-terminal --geometry 75x25+50+50 -x bash -c "pwd; read -p quit? -n1" &
  
  gnome-terminal --geometry 75x25+0+0 -x bash -c "socket_thread.py -server; read -p quit? -n1" &
  gnome-terminal --geometry 75x25+500+0 -x bash -c "socket_thread.py -client 1; read -p quit? -n1" &
  gnome-terminal --geometry 75x25+1000+0 -x bash -c "socket_thread.py -client 2; read -p quit? -n1" &
  gnome-terminal --geometry 75x25+1000+0 -x bash -c "socket_thread.py -client 3; read -p quit? -n1" &
  gnome-terminal --geometry 75x25+1000+0 -x bash -c "socket_thread.py -client 4; read -p quit? -n1" &
  gnome-terminal --geometry 75x25+1000+0 -x bash -c "socket_thread.py -client 5; read -p quit? -n1" &
  gnome-terminal --geometry 75x25+1000+0 -x bash -c "socket_thread.py -client 6; read -p quit? -n1" &
  cat /etc/services | grep 2222
  echo '***** end launch socket *****'

"""

import sys
import socket
import threading
import time

_ENDMESSAGE_ = "__endSocket__"
_PORT_ = 22222

# http://apprendre-python.com/page-reseaux-sockets-python-port

def getTime():
  return time.strftime("%Hh %Mm %Ss")

def getLocDist(aSocket):
  return "loc_addr %s dist_addr %s" % (aSocket.getsockname(), aSocket.getpeername())

def toHex(aStr):
  res = ""
  for i in aStr: res += hex(ord(i)) + " "
  return res

#########################################
# server
#########################################
class ClientThread(threading.Thread):
  def __init__(self, ip, port, clientsocket):
    threading.Thread.__init__(self)
    self.ip = ip
    self.port = port
    self.clientsocket = clientsocket
    print("Thread Server: [+] Nouveau thread pour %s" % getLocDist(self.clientsocket))

  def run(self):
    print("Thread Server: connection de %s %s" % (self.ip, self.port))
    while True: 
      receive = self.clientsocket.recv(512)
      print("Thread Server: receive from %i: '%s'" % (self.port, receive))
      if receive == _ENDMESSAGE_: break
      res = "reply from server with '%s' at %s" % (receive, getTime())
      self.clientsocket.send(res)
    print("Thread Server: end thread by client %i deconnected on %s" % (self.port, receive))

def runServer():
  tcpsock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
  tcpsock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
  tcpsock.bind(("", _PORT_))

  while True:
    tcpsock.listen(10)
    print("Server: listening...")
    (clientsocket, (ip, port)) = tcpsock.accept()
    newthread = ClientThread(ip, port, clientsocket)
    newthread.start()

#########################################
# client
#########################################
def runClientNo(no):
  runClient("Client%s" % no)

def runClient(name):
  s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
  while True:
    try:
      s.connect(("", _PORT_))
      break
    except:
      print("%s: can't connect, wait and retry"  % name)
      time.sleep(1)
    
  print("%s: connected to server %s" % (name, getLocDist(s)))

  for i in range(10):
    message = "%s YOO %0d" % (name, i)
    s.send(message)
    print("%s: send: '%s'" % (name, message))
    received = s.recv(9999999)
    print("%s: receive: '%s' as '%s'" % (name, received, toHex(received)))
    time.sleep(1)
  message = _ENDMESSAGE_
  print("%s: send: '%s'" % (name, message))
  s.send(message)
  print('%s: end of client' % name)

if __name__ == '__main__':
  if len(sys.argv) < 2:
    print(__doc__)
  else:
    if sys.argv[1] == "-server":
      runServer()
    elif sys.argv[1] == "-client":
      runClientNo(sys.argv[2])  
    else:
      print(_help)
    print('end of %s' % __file__)



