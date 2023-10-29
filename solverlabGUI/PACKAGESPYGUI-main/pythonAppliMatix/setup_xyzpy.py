#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END

from distutils.core import setup
import os
from os.path import isfile, join

"""
https://docs.python.org/2/distutils/introduction.html#a-simple-example:

To create a source distribution for this module, 
you would create a setup script, setup.py, 
containing the above code, and run this command from a terminal:
python setup.py sdist

sdist will create an archive file (e.g., tarball on Unix, ZIP file on Windows) 
containing your setup script setup.py, 
and your module foo.py. 
The archive file will be named foo-1.0.tar.gz (or .zip), 
and will unpack into a directory foo-1.0.

If an end-user wishes to install your foo module, 
all she has to do is download foo-1.0.tar.gz (or .zip), 
unpack it, and—from the foo-1.0 directory—run

python setup.py install

which will ultimately copy foo.py to the appropriate directory for 
third-party modules in their Python installation.

This simple example demonstrates some fundamental concepts of the Distutils. 
First, both developers and installers have the same basic user interface, 
i.e. the setup script. 
The difference is which Distutils commands they use: 
the sdist command is almost exclusively for module developers, 
while install is more often for installers 
(although most developers will want to install their own code occasionally).

If you want to make things really easy for your users, 
you can create one or more built distributions for them. 
For instance, if you are running on a Windows machine, 
and want to make things easy for other Windows users, 
you can create an executable installer 
(the most appropriate type of built distribution for this platform) 
with the bdist_wininst command. For example:

python setup.py bdist_wininst

will create an executable installer, foo-1.0.win32.exe, in the current directory.

Other useful built distribution formats are RPM, 
implemented by the bdist_rpm command, Solaris pkgtool (bdist_pkgtool), 
and HP-UX swinstall (bdist_sdux). 
For example, the following command will create an RPM file called foo-1.0.noarch.rpm:

python setup.py bdist_rpm
"""

def mySetup(): #for sphynx autodoc -> set as def with if __name__ == '__main__'
  verbose = True
  directory, name = os.path.split( os.path.realpath( __file__ ) )
  if os.getcwd() != directory:
    raise Exception("ERROR: '%s' have to be launch in current directory '%s'" % (name, directory))
  if name != "setup.py":
    raise Exception("ERROR: '%s' have to be launch as renamed 'setup.py'" % (name))
  listDir = [x[0] for x in os.walk("xyzpy")] #have to be relative

  #dictFile = {k: ["*.txt"] for k in listDir} #to keep all files .txt etc... python  V2.7
  dictFile = dict( ( k, [ f for f in os.listdir(k) if isfile(join(k,f)) ] ) for k in listDir)

  if verbose:
    for i in listDir: print("package for dirs:", i)
    print("package for files:", dictFile)

  setup(
        name ='xyzpy', # Application name
        version = '1.0',
        description = 'Python package pythonAppliMatix xyzpy',
        license = os.path.join("..", "LICENSE.txt"),
        #long_description = open(os.path.join(directory, "README.txt")).read(),
        long_description = "subpackage of pythonAppliMatix",
        author = 'Christian Van Wambeke',
        author_email = 'christian.van-wambeke@cea.fr',
        packages = listDir,
        package_data = dictFile, #to keep all files .sh etc...
        # Dependent packages (distributions)
        install_requires = [ "numpy" ], #example
      )

if __name__ == '__main__':
  mySetup()