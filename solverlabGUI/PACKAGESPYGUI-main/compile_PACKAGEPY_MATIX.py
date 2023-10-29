#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END

"""
SALOME-SAT compile script for PACKAGESPY
"""

#import sys
import os
import subprocess
import common



def printLocalEnv(param="None"):
  with open('/tmp/environ.tmp', 'w') as f:
    f.write(param+'\n')
    for i in sorted(os.environ):
      f.write(i+" : "+os.environ[i]+"\n")

def printConfig(config):
  with open('/tmp/config_fftwcode.tmp', 'w') as f:
    for i in config:
      f.write("\n%s: %s\n" % (i, config.getByPath(i)))

def compil_example(config, helper, logger):
  printLocalEnv('aCompilation')

# This script is used to build the application module.
# First, it copies the content of the sources directory to the install directory.
# Then it runs 'lrelease' to build the resources.

def compil(config, helper, logger):
    helper.prepare()
    helper.results.buildconfigure = True
    helper.results.configure = True
    helper.results.make = True
    helper.results.cmake = True
    helper.results.ctest = True

    environ = helper.build_environ.environ.environ

    #it is done in SOURCE/PACKAGESPY to git push and have uptodate doc in sources
    command = "cd %s; makeDoc.sh" % config.TOOLS.common.module_info.PACKAGESPY.source_dir
    res = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=environ).communicate()
    #print "\nINFO PACKAGEPY:",config.TOOLS.common.module_info.PACKAGESPY.source_dir
    if res[1] != "": #an error occured
      if not "WARNING" in res[1]:
        helper.results.install = False
        helper.results.check = False
        print("ERROR:",res[1])
        helper.log_file.write(res[1]+"\n")
        return helper.results

    if not helper.source_dir.smartcopy(helper.install_dir):
        raise common.SatException(_("Error when copying %s sources to install dir") % helper.module)

    helper.results.install = True
    helper.results.check = True
    return helper.results

    """ obsolete
    # test lrelease #.pyconf needs in ..._APPLI pre_depend : ['qt']
    command = "which lrelease"
    res = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=environ).communicate()
    if res[1] != "": #an error occured
      helper.results.install = False
      helper.results.check = False
      print "ERROR:",res[1]
      helper.log_file.write(res[1]+"\n")
      return helper.results

    # run lrelease
    command = "lrelease *.ts"
    res = subprocess.call(command,
                          shell=True,
                          cwd=str(helper.install_dir + "resources"),
                          env=environ,
                          stdout=helper.log_file,
                          stderr=subprocess.STDOUT)

    helper.results.install = True
    helper.results.check = True

    return helper.results
    """
