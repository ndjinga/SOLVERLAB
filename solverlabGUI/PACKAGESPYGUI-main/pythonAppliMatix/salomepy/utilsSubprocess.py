#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


"""
RunCmd is to launch subprocess.Popen with timeout
"""

import os
import sys
import subprocess as SP
import threading
import time
import traceback
import pprint as PP

verbose = False

_expected_terminals = "mate-terminal gnome-terminal konsole".split()

# terminal colors
CSI = '\033['
BLACK           = CSI + str(30) + 'm'
RED             = CSI + str(31) + 'm'
GREEN           = CSI + str(32) + 'm'
YELLOW          = CSI + str(33) + 'm'
BLUE            = CSI + str(34) + 'm'
MAGENTA         = CSI + str(35) + 'm'
CYAN            = CSI + str(36) + 'm'
WHITE           = CSI + str(37) + 'm'
RESET           = CSI + str(39) + 'm'


def is_stdout_tty():
  # in salome sys.stdout.isatty() -> 'PyOut' object has no attribute 'isatty'
  try:
    res = sys.stdout.isatty()
  except:
    res = False
  if verbose: sys.stdout.write("\nsys.stdout.isatty() %s\n" % res)
  return res


class RunCmd(threading.Thread):
    def __init__(self, cmd, timeOut):
        super(RunCmd, self).__init__()
        self.cmd = cmd
        self.timeOut = timeOut
        self.timeOutReached = None

    def run(self):
        self.p = SP.Popen(self.cmd, shell=True)
        self.p.wait()

    def RunWithTimeOut(self):
        self.start()
        self.join(self.timeOut)
        if self.is_alive():
            #res = self.p.communicate()
            self.p.terminate()      #use self.p.kill() if process needs a kill -9
            self.p.kill()
            self.join()
            if verbose: print("timeOut reached")
            self.timeOutReached = True
        else:
            #res = self.p.communicate()
            if verbose: print("timeOut NOT reached")
            self.timeOutReached = False

def getRunCastemDefault():
  test = ["CASTEM4CASSIS", "CASTEMDEFAULT"]
  for k in test:
    cmd = os.getenv(k)
    if cmd != None: return cmd
  print("Castem command not found in environment variables: '%s'" % test)
  return None
  
def runCastem(name, srcDir, destDir, timeOut, runCastem=None, resultFile=None):
  """
  launch a case castem on name.dgibi
  copy all tree srcDir to destDir and run castem command 'runCastem name.dgibi'
  with timeOut in seconds
  outputs of run castem are in castem_stdout.log and castem_stderr.log
  """
  if runCastem == None:
    runCastemCurrent = getRunCastemDefault()
  else:
    runCastemCurrent = runCastem
  
  if runCastemCurrent == None:
    return "ko: run castem command not found"

  import salomepy.utilsWorkdir as UTW
  UTW.copyDir(srcDir, destDir)
  
  #have to create fin_tmp.dgibi input to solve if gibiane error and "lecture continue sur clavier"
  nameFile =  os.path.join(destDir, "fin_tmp.dgibi")
  with open(nameFile, 'w') as f:
    f.write("* FILE CREATED by utilsSubprocess.py.runCastem\nMESS 'EXIT ON GIBIANE ERROR';\nFIN;\nFIN;\n")
  
  #cmd = "cd %s ; $CASTEM/runCastem %s.dgibi < fin.dgibi 1> castem_stdout.log 2> castem_stderr.log" % (destDir, name)
  cmd = "cd %s ; %s %s.dgibi < fin_tmp.dgibi 1> castem_stdout.log 2> castem_stderr.log" % (destDir, runCastemCurrent, name)
  if verbose: print("\nlaunch castem: '%s'" % cmd)
  
  aRun = RunCmd(cmd, timeOut)
  aRun.RunWithTimeOut()
  if aRun.timeOutReached == True:
    return "ko: timeout reached: %s" % aRun.timeOut
  if verbose: print("runCastem: timeOutReached: ", aRun.timeOutReached)
  if resultFile == None:
    return "ok" #no test, supposed to be good
  else: #a result file defined
    #fileRes = os.path.join(destDir, name + ".sauv")
    fileRes = os.path.join(destDir, resultFile)
    res = os.path.isfile( fileRes )
    if verbose: print("runCastem: test result file: %s %s" % (res, fileRes))
    if not res:
      return "ko: Result file not found: %s" % fileRes
    return "ok"


def cmdline(cmd, verbose=False):
  process = SP.Popen( args=cmd, stdout=SP.PIPE, stderr=SP.PIPE, shell=True )
  out, err = process.communicate()
  if len(err) == 0: # python3 returns b''
    return out
  else:
    if verbose: print("ERROR: cmd line '%s' -> '%s'" % (cmd, err), len(err))
    return None


def getCurrentTerminal():
  """with linux could be mate-terminal gnome-terminal konsole"""
  terminal = os.getenv("MATIX_TERMINAL", None)
  for terminal in _expected_terminals:
    out = cmdline("which " + terminal, True)
    if verbose: print()
    if out is not None: return terminal
  print("ERROR: getCurrentTerminal not found in %s" % _expected_terminals)
  return _expected_terminals[0]  # as default to try on error


class RunCmdInTerminal(object):
  """
  launch a cmd/script in a linux terminal
  with log trace stdout/stderr as tail --follow file.log
  functionality, without freezing caller (as a PyQt mainwindow)
  """
  def  __init__(self, title="RunCmdInTerminal"):
    if verbose: print("RunCmdInTerminal __init__ TRACEBACK\n", "".join(traceback.format_stack()))
    self._logExternalTerminal = True
    self._fileExternalTerminal = None
    self._processExternalTerminal = None
    self._terminal = getCurrentTerminal() # mate-terminal gnome-terminal konsole
    self._title = title
    # 120x30+100+20 is 120 caracters x 30 lines + 100 pixelfromleft + 20 pixelfromup
    self._cmdExternalTerminal = "%s --geometry 120x30+100+20" % self._terminal
    if self._terminal in "mate-terminal".split(" "):
      self._cmdExternalTerminal += " --title=%s" % self._title

  def __del__(self):
    print((MAGENTA+"EXIT '%s'\n"+RESET) % self._title)

  def run(self, cmd, workdir="/tmp"):
    """
    launch cmd as bash file command(s) with stdout/stderr redirected in file.log
    this log file is displayed/followed with linux command 'tail -f'
    """
    if self._fileExternalTerminal is None:
      _workdir = os.path.expandvars(workdir)
      _TIMETAG = time.strftime('%y%m%d_%H%M%S')
      t = _TIMETAG
      f1 = _workdir + "/run_%s.log" % t
      self._fileExternalTerminal = f1

      # "{0} love {1}".format("GeeksforGeeks", "Geeks")

      fcmd = _workdir + "/terminal_cmd_%s.bash" % t
      flog = _workdir + "/terminal_log_%s.bash" % t

      cmdlog = """\
#!/usr/bin/bash

# could be launch by 'mate-terminal --geometry 120x50+200+100 --command {0}'

echo "current directory (may be there is symbolic link)"
pwd
echo "follow file {1}"
echo
tail --lines=100 --retry --follow=name {1}
# precaution for user eyes if unexpected exiting
sleep 10
""".format(flog, f1)

      with open(flog, "w") as f:
        f.write(cmdlog)
      st = os.stat(flog)
      os.chmod(flog, st.st_mode | 0o111) # executable

      cmduser = """\
#!/usr/bin/bash

echo "{1}#### Begin ExternalTerminal command(s){2}"

# user command(s)
{0}

echo "{1}#### End ExternalTerminal command(s){2}"
""".format(cmd, MAGENTA, RESET)

      with open(fcmd, "w") as f:
        f.write(cmduser)
      st = os.stat(fcmd)
      os.chmod(fcmd, st.st_mode | 0o111)  # executable

      if self._terminal in "konsole".split():
        cmd_terminal = self._cmdExternalTerminal + " -e " + flog
      else:  # mate-terminal gnome-terminal
        cmd_terminal = self._cmdExternalTerminal + " --command " + flog

      tit = "BEGIN '%s'" % self._title
      msg = "command log in terminal is\n  '%s'\nlog file is\n  '%s'" % \
            (cmd_terminal, f1)

      print(MAGENTA+"\n"+tit+RESET)  # tty colored for eyes
      # if is_stdout_tty(): print(MAGENTA+"\n"+tit+RESET)  # tty colored for eyes
      print(msg)

      self._processExternalTerminal = SP.Popen(cmd_terminal, shell=True, cwd=_workdir)
      if verbose: print("SP.Popen done for windows terminal")

      # redirect stdout and stderr on scruted file.log (tail -f)
      cmduser_piped = "{0} &> {1}".format(fcmd, f1) # redirect stdout and stderr on scruted file.log (tail -f)
      self._processExternalTerminal = SP.Popen(cmduser_piped, shell=True, cwd=_workdir)
      if verbose: print("SP.Popen done for command\n%s" % MAGENTA+cmd+RESET)

  def coutInExternalTerminal(self, mess):
    with open(self._fileExternalTerminal, "a") as f:
      f.write(mess + '\n')


def test():
  aCmd = RunCmdInTerminal("RunCmdInTerminal_test")
  cmd = """
set -x 
pwd
echo $USER
set +x 
"""
  aCmd.run(cmd, workdir="${HOME}")
  time.sleep(10)
  print("print red helllooo")
  aCmd.coutInExternalTerminal(RED+"HELLLOOOOOOO"+RESET)


"""
python
import utilsSubprocess as USP
USP.test()
"""
if __name__ == '__main__':
  # for test
  test()
