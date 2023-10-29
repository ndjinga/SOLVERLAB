#!/usr/bin/env python
#-*- coding:utf-8 -*-

_help = r"""
Named Pipes (Fifos)

On some platforms, it is also possible to create a pipe that exists as a file. 
Such files are called "named pipes" (or sometimes, "fifos"), 
because they behave just like the pipes created within the previous programs, 
but are associated with a real file somewhere on your computer, 
external to any particular program. 
Once a named pipe file is created, processes read and write it using normal file operations. 
Fifos are unidirectional streams, but a set of two fifos can be used to implement bidirectional communication just as we did for anonymous pipes in the prior section.

Because fifos are files, they are longer-lived than in-process pipes 
and can be accessed by programs started independently. 
The unnamed, in-process pipe examples thus far depend on the fact 
that file descriptors (including pipes) are copied to child processes. 
With fifos, pipes are accessed instead by a filename visible to all programs 
regardless of any parent/child process relationships. 
Because of that, they are better suited as IPC mechanisms 
for independent client and server programs; 
for instance, a perpetually running server program may create and listen 
for requests on a fifo, that can be accessed later by arbitrary clients 
not forked by the server.

In Python, named pipe files are created with the os.mkfifo call, 
available today on Unix-like platforms and Windows NT (but not on Windows 95/98). 
This only creates the external file, though; 
to send and receive data through a fifo, it must be opened 
and processed as if it were a standard file. 
Example 3-19 is a derivation of the pipe2.py script listed earlier, 
written to use fifos instead of anonymous pipes.
Example 3-19. PP2E\System\Processes\pipefifo.py

To distinguish messages better, we can mandate a separator character in the pipe. 
An end-of-line makes this easy, because we can wrap the pipe descriptor 
in a file object with os.fdopen, and rely on the file object's readline method 
to scan up through the next \n separator in the pipe. 
Example implements this scheme. 

Because the fifo exists independently of both parent and child, 
there's no reason to fork here -- 
the child may be started independently of the parent, 
as long as it opens a fifo file by the same name. 
Here, for instance, on Linux the parent is started in one xterm window, 
and then the child is started in another. 
Messages start appearing in the parent window only after the child is started

see http://flylib.com/books/en/2.723.1.41/1/

file ./pipefifo
     ./pipefifo: fifo (named pipe)

to launch:
  >> essai_names_pipe_3.py -parent  #in one console
  >> ls -alt
  >> essai_names_pipe_3.py -child   #in another console

  or...
  >>> LaunchDemo_3.bash
"""

_launchDemo3=r"""
#!/usr/bin/env bash

#this is demonstration using gnome-terminal and essai_names_pipe_3.py

#create bash scripts for gnome-terminal --command [bash scripts]
cat > launchParent.bash << EOL
#!/usr/bin/env bash
echo '***** begin launchParent *****'
pwd
ls -alt
file ./pipefifo
echo 'essai_names_pipe_3.py -parent'
essai_names_pipe_3.py -parent
echo '***** end launchParent.bash in 10 sec *****'
sleep 10

EOL
chmod a+x launchParent.bash

cat > launchChild.bash << EOL
#!/usr/bin/env bash
echo '***** begin launchChild *****'
pwd
ls -alt
file ./pipefifo
echo 'essai_names_pipe_3.py -child'
essai_names_pipe_3.py -child
echo '***** end launchChild.bash in 10 sec *****'
sleep 10

EOL
chmod a+x launchChild.bash

#find . -name "launch*.bash" -exec cat -n {} \;
more l*.bash | cat

#launch demo
gnome-terminal --command launchParent.bash && gnome-terminal --command launchChild.bash

echo '***** end launchDemo *****'

"""

#########################################################
# named pipes; os.mkfifo not avaiable on Windows 95/98;
# no reason to fork here, since fifo file pipes are 
# external to processes--shared fds are irrelevent;
#########################################################

import os, time, sys
fifoname = './pipefifo'                          # must open same name

def child():
    print('begin child()')
    mess = 'Hello from child %03d\n'
    pipeout = os.open(fifoname, os.O_WRONLY)     # open fifo pipe file as fd
    print('Child open pipeout done')
    zzz = 100 
    for i in range(zzz):                         # roll 0 to zzz
        time.sleep(1)
        toSend = mess % i
        #print 'Child %d WILL send "%s" at %s' % (os.getpid(), toSend[:-1], time.strftime("%Hh %Mm %Ss"))
        os.write(pipeout, toSend)
        #toSend[:-1] as remove '\n'
        print('Child %d send "%s" at %s' % (os.getpid(), toSend[:-1], time.strftime("%Hh %Mm %Ss")))
    print('end chid()')
    return
         
def parent():
    print('begin parent()')
    pipein = open(fifoname, 'r')                 # open fifo as stdio object
    print('Parent open pipein done')
    while True:
        line = pipein.readline()[:-1]            # blocks until data sent
        print('Parent %d got "%s" at %s' % (os.getpid(), line, time.strftime("%Hh %Mm %Ss")))
        if line == "": break
    print('end parent()')

def createLaunchDemo3():
    name = "LaunchDemo_3.bash"
    with open(name, 'w') as f: f.write(_launchDemo3)
    chmodax(name)

def chmodax(aFile):
    """chmod a+x path"""
    st = os.stat(aFile)
    os.chmod(aFile, st.st_mode | 0o111)

if __name__ == '__main__':
    if len(sys.argv) != 2:
      print(_help)
      createLaunchDemo3()
    else:
      if not os.path.exists(fifoname):
          os.mkfifo(fifoname)                    # create a named pipe file
      if sys.argv[1] == "-parent":
          parent()                               # run as parent process
      elif sys.argv[1] == "-child":
          child()                                # run as child process
      else:
          print(_help)
      print('end of all')


#for info
r"""
>> gnome-terminal --help-all
  Utilisation :
    gnome-terminal [OPTION...]

  Options de l'aide :
    -h, --help                         Affiche les options de l'aide
    --help-all                         Affiche toutes les options de l'aide
    --help-gtk                         Affiche les options GTK+
    --help-terminal                    Afficher les options du terminal
    --help-window-options              Afficher les options par fenêtre
    --help-terminal-options            Afficher les options par terminal

  Options GTK+
    --class=CLASSE                     Classe du programme telle qu'utilisée 
                                         par le gestionnaire de fenêtres
    --name=NOM                         Nom du programme tel qu'utilisé 
                                         par le gestionnaire de fenêtres
    --gdk-debug=DRAPEAUX               Drapeaux de débogage GDK à définir
    --gdk-no-debug=DRAPEAUX            Drapeaux de débogage GDK à ne pas définir
    --gtk-module=MODULES               Charge des modules GTK+ additionnels
    --g-fatal-warnings                 Rend tous les avertissements fatals
    --gtk-debug=DRAPEAUX               Drapeaux de débogage GTK+ à définir
    --gtk-no-debug=DRAPEAUX            Drapeaux de débogage GTK+ à ne pas définir

  Options to open new windows or terminal tabs; more than one of these may be specified:
    --window                           Ouvre une nouvelle fenêtre contenant 
                                         un onglet avec le profil par défaut
    --tab                              Ouvre un nouvel onglet dans la dernière fenêtre 
                                         ouverte avec le profil par défaut

  Window options; if used before the first --window or --tab argument, sets the default for all windows:
    --show-menubar                     Afficher la barre de menu
    --hide-menubar                     Masquer la barre de menu
    --maximize                         Maximise la fenêtre
    --full-screen                      Fenêtre en plein écran
    --geometry=GÉOMÉTRIE               Définit la taille de la fenêtre ; 
                                         par exemple : 80x24 ou 80x24+200+200 
                                        (COLONNESxLIGNES+X+Y)
    --role=RÔLE                        Définit le rôle de la fenêtre
    --active                           Rend le dernier onglet indiqué actif

  Terminal options; if used before the first --window or --tab argument, sets the default for all terminals:
    -e, --command                      Exécute le paramètre de cette option dans le terminal
    --profile=NOMDUPROFIL              Utilise le profil indiqué au lieu du profil par défaut
    -t, --title=TITRE                  Définit le titre du terminal
    --working-directory=RÉPTRAVAIL     Définit le répertoire de travail
    --zoom=ZOOM                        Définit le facteur de zoom du terminal 
                                         (1.0 = taille normale)

  Options de l'application :
    --load-config=FICHIER              Charge un fichier de configuration pour le terminal
    --display=AFFICHAGE                Affichage X à utiliser

  Émulateur de terminal GNOME
"""

