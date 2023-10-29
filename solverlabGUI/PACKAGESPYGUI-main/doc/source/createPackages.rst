ALL THIS IS OBSOLETE
===========================

Never implemented, may be in future


Création des packages
------------------------------------------

La création d'un package particulier consiste à créer un fichier xxx.zip
contenant un **sous-ensemble** de l'arborescence PACKAGEPY
telle que clonée de sa base git initiale.

Il suffit de détruire l'ensemble des répertoires et fichiers inutiles à l'application choisie.

Il y a un script python "createPackage.py" qui s'occupe de cela:

.. code-block:: bash

  > cd $MYDIR
  > packagespy/createPackage.py

  create a package python from a git repo
  delete .git and some not used directories and files

  example (in python):
    import createPackage
    params = {}
    params["packageName"] = "appliMatix"
    params["tagName"] = "HEAD"
    CR = createPackage.CreatePackage( params )
    CR.go()

  example (in shell):
    createPackage.py appliMatix HEAD

  existing known packages: ['appliMatix', 'oscarPandasTest']

  > packagespy/createPackage.py oscarPandasTest

  package with
  CreatePackage instance:
  params: {
   date: 2015-01-20 13:35:17
   deleteDirs: ['crescendopy', 'pythonHelpDistene', 'pythonTestImport', 'pythonCloneRename', '.eric4project']
   deleteDirsVersioning: ['.git', '.CVS']
   deleteFiles: ['setPathMaison.sh', 'packagespy.e4p', '.spyderproject']
   gitRepo: /home/matix/GitRepo/packagespy.git
   hostname: is206786
   makeArchive: oscarPandasTest
   packageName: oscarPandasTest
   skipDirs: ['.git', '.CVS']
   tagName: HEAD
   tempDir: /tmp
   user: matix }



  ERROR: cloneDir temporary target directory exists yet.
   rm -rf /tmp/createPackage (y/n)? y
  cmd: cd /tmp/createPackage; git clone /home/matix/GitRepo/packagespy.git
  Initialized empty Git repository in /tmp/createPackage/packagespy/.git/
  cmd: cd /tmp/createPackage/packagespy ; git checkout HEAD
  remove dir  /tmp/createPackage/packagespy/.git
  remove dir  /tmp/createPackage/packagespy/pythonCloneRename
  remove dir  /tmp/createPackage/packagespy/pythonHelpDistene
  remove dir  /tmp/createPackage/packagespy/pythonTestImport
  remove dir  /tmp/createPackage/packagespy/pythonAppliMatix/crescendopy
  remove file /tmp/createPackage/packagespy/.spyderproject
  remove file /tmp/createPackage/packagespy/packagespy.e4p
  remove file /tmp/createPackage/packagespy/setPathMaison.sh

  12 directories
  making archive of '/tmp/createPackage/packagespy' ...
  archive '/tmp/createPackage/oscarPandasTest.zip' done



Mise en oeuvre des packages
------------------------------------------

Voir l'exemple :ref:`OSCAR oscarPandasTest <oscarPPY>` Présentation de l’IHM OSCAR.
