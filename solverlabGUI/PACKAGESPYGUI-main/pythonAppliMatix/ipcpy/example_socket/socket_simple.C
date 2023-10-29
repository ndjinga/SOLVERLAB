
/*
this is demonstration using gnome-terminal and ROOT sockets

to launch: root<->root
  gnome-terminal --geometry 75x25+0+0 -x bash -c "root -q -l socket_simple.C -server ; read -p quit? -n1" &
  sleep 1
  gnome-terminal --geometry 75x25+500+0 -x bash -c "root -q -l socket_simple.C -client ; read -p quit? -n1" &

or to launch: root<->python mix
  gnome-terminal --geometry 75x25+0+0 -x bash -c "root -q -l socket_simple.C -server ; read -p quit? -n1" &
  sleep 1
  gnome-terminal --geometry 75x25+500+0 -x bash -c "socket_thread.py -client 1; read -p quit? -n1" &
  # seems prependin header plus message recive in python

see:
  https://root.cern.ch/root/html534/guides/users-guide/Networking.html#a-server-with-multiple-sockets
  $ROOTSYS/tutorials/net/hserv.C 
  $ROOTSYS/tutorials/net/hclient.C.
*/


void runServer() {  
  /* 
  server
  */
  cout << "BEGIN runServer()" << endl;
  TServerSocket *ss = new TServerSocket(22222, kTRUE);
  TSocket *socket = ss->Accept();
  cout << "Server: accept client" << endl;
  /*
  TString message = "Hello from server\n";
  socket->Send(message);
  cout << "Server: send '" << message << "'" << endl;
  */
  char message2[64] = "AAAHello from server ROOT char!AAA";
  socket->Send(message2);
  cout << "Server: send '" << message2 << "'" << endl;
  cout << "socket->GetBytesSent() '" << socket->GetBytesSent() << "'" << endl;
  ss->Close();
}

void runClient() {
  /* 
  client
  On the clientside, we create a socket and ask the socket to receive input.
  */
  cout << "BEGIN runClient()" << endl;
  TSocket *socket = new TSocket("localhost", 22222);
  TMessage *message;

  socket->Recv(message);
  if (message->What() == kMESS_STRING) {
    char str[64];
    cout << "Client: as 'kMESS_STRING'" << endl;
    cout << "Client: receive '" << message->ReadString(str, 64) << "'" << endl;
  }
  else {
    cout << "Client: receive what '" << message->->What() << "' unknown TODO " << endl;
  }

}

void socket_simple() {
  TString nameFile = gSystem->BaseName(__FILE__);
  cout << "BEGIN " << nameFile << endl;
  int argc = gApplication->Argc();
  for( int ai=0; ai<argc; ai++ ) printf("  arg[%d]='%s'\n", ai, gApplication->Argv(ai));
  TString arg = gApplication->Argv(argc - 1);
  if (arg == "-server") runServer(); 
  if (arg == "-client") runClient();
  cout << "END " << nameFile << endl;
}
