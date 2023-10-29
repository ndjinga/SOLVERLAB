#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


"""
This class is used to launch salome in batch mode... or not
It launches the necessary server to use GEOM and co.
then close the servers in destructor

goal is easy one-line access functions
to salome, study, builder, geom, smesh, GUI etc...

example of use xsalomesession.py
import salomepy.xsalomesession as XSS
geompy = XSS.getGeompy()
smesh = XSS.getSmesh()
geomBuilder = XSS.getGeomBuilder()
smeshBuilder = XSS.getSmeshBuilder()
geomGUI = XSS.getGeomGUI() #None if not GUI mode
desktop = XSS.getDesktop() #None if not GUI mode

expert mode:
from salomepy.xsalomesession import XSalomeSession
x_salome_session = XSalomeSession(modules=["GEOM", "SMESH", "MED", etc...])
x_salome_session =  XSS.getXSalomeSession()
"""

import os
import glob
import pprint as PP
import xyzpy.loggingXyz as LOG
logger = LOG.getLogger()

verbose = False

_x_salome_session = None
_defaultModules = ["GEOM", "SMESH", "MEDCOUPLING"]
_nbSessions = 0
_nbStudies = 0
_messDone = []
_nbTriesLaunchSalome = 0


def set_OMNIORB_CONFIG():
  """
  set environment variable OMNIORB_CONFIG if not existing
  used for import pvsimple for paraview server SALOME 780
  something like export OMNIORB_CONFIG='/tmp/.omniORB_wambeke_is210329_2813.cfg'
  """
  OMNIORB_CONFIG = os.getenv("OMNIORB_CONFIG")
  if OMNIORB_CONFIG == None:
    user = os.getenv("USERNAME")
    hostname = os.getenv("HOSTNAME").split(".")[0] #is210329 from is210329.intra.cea.fr
    pathname = "/tmp/.omniORB_%s_%s_*.cfg" % (user, hostname)
    files = glob.glob(pathname)
    if len(files) == 0:
      logger.warning("OMNIORB_CONFIG not set because no files '%s'" % pathname)      
    os.environ["OMNIORB_CONFIG"] = files[0]
    logger.warning("set OMNIORB_CONFIG='%s'" % os.getenv("OMNIORB_CONFIG"))
  return
  
def get_default_OMNIORB_CONFIG(port=2811):
  """SALOME 780"""
  user = os.getenv("USERNAME")
  hostname = os.getenv("HOSTNAME").split(".")[0] #is210329 from is210329.intra.cea.fr
  name = "/tmp/.omniORB_%s_%s_%s.cfg" % (user, hostname, port)
  return name
  
  
class XSalomeSession(object):

    def __init__(self, modules=None, showDesktop=False):
        global _nbSessions
        if verbose: print("XSalomeSession.__init__")
        #precaution for the moment
        if _nbSessions > 0:
           raise Exception("only one XSalomeSession session %s" % _nbSessions)
        _nbSessions += 1
        if modules==None:
           self.modules = _defaultModules
        else:
           self.modules = modules
        self.showDesktop = showDesktop
        self.geompy = None
        self.smeshpy = None
        self.geomBuilder = None
        self.smeshBuilder = None
        self.geomGUI = None
        self.sgPyQt = None
        self.sg = None
        self.sgDesktop = None
        self.port = None
        if not isSalomeAlreadyLaunched() or self.needRelaunch():
           self.launchSalome()

    def __repr__(self):
        return "XSalomeSession(id: %s, port: %s)" % (id(self), self.port)

    def needRelaunch(self):
        need = False
        if "GEOM" in self.modules:
            try: #not try: important to know why...
                # imports pour Salome
                import salome
                salome.salome_init()
                if salome.myStudy == None:
                   self.createNewStudy()
                import GEOM
                from salome.geom import geomBuilder
                self.geomBuilder = geomBuilder
                self.geompy = geomBuilder.New()
                #self.geompy.MakeVertex
            except:
                print("Salome has been killed. It needs to be launched again")
                need = True
            #print "need: ", need
        return need

    def launchSalome(self):
        print("XSalomeSession.launchSalome")
        modules = self.modules
        # Salome params
        import sys
        argv_ini = sys.argv
        sys.argv  = ['bin/salome/runSalome.py']
        sys.argv += ["-t"]
        if self.showDesktop:
          sys.argv += ["--gui"]
        #do not work sys.argv += ["--embedded=study"]
        if len(modules):
            sys.argv += ['--modules=%s'%(",".join(modules))]
            pass
        if verbose: print("sys.argv: ", sys.argv)
        # Launch Salome
        import runSalome
        clt, d = runSalome.main()
        port = d['port']
        self.port = port
        print("Salome launched on port %s with sys.argv %s" % (port, str(sys.argv)))
        sys.argv = argv_ini
        return

    def getGeompy(self):
        if "GEOM" in self.modules:
            try:
                # imports pour Salome
                import salome
                salome.salome_init()
                if salome.myStudy == None:
                    self.createNewStudy()

                import GEOM
                from salome.geom import geomBuilder
                self.geomBuilder = geomBuilder
                self.geompy = geomBuilder.New()
                return self.geompy
            except:
                print("WARNING: Problem in Salome for GEOM")
                return None
        else:
            print("WARNING: XSalomeSession.getGeompy: Salome without GEOM")
            return None

    def getSmesh(self):
        if "SMESH" in self.modules:
            try:
                # imports pour Salome
                import salome
                salome.salome_init()
                if salome.myStudy == None:
                    self.createNewStudy()

                import SMESH
                from salome.smesh import smeshBuilder
                self.smeshBuilder = smeshBuilder
                self.smeshpy = smeshBuilder.New()
                return self.smeshpy
            except:
                print("WARNING: Problem in Salome for SMESH")
                return None
        else:
            print("WARNING: XSalomeSession.getSmesh: Salome without SMESH")
            return None

    def getGeomBuilder(self):
        import salome
        from salome.geom import geomBuilder
        self.geomBuilder = geomBuilder
        return self.geomBuilder

    def getSmeshBuilder(self):
        import salome
        from salome.smesh import smeshBuilder
        self.smeshBuilder = smeshBuilder
        return self.smeshBuilder

    def getGeomGUI(self):
        try:
          import salome
          self.geomGUI = salome.ImportComponentGUI("GEOM")
          return self.geomGUI
        except:
          return None

    def getSmeshGUI(self):
        try:
          import salome
          self.smeshGUI = salome.ImportComponentGUI("SMESH")
          return self.smeshGUI
        except:
          return None

    def getDesktop(self):
        if True: #try:
            import SalomePyQt
            self.sgPyQt = SalomePyQt.SalomePyQt()
            import libSALOME_Swig
            self.sg = libSALOME_Swig.SALOMEGUI_Swig()
            self.sgDesktop = self.sgPyQt.getDesktop()
            return self.sgDesktop
        else: #except:
            self.sgDesktop = None
            print("WARNING: XSalomeSession.getDesktop: not found")
            return None

    def createNewStudy(self, name=None):
      """
      called to clean active study: close old study before
      return the SALOMEDS._objref_Study instance
      risky because users of geompy etc...
      have to "renewed their geompy" if they are in global
      geompy = XSS.getGeompy() is not dynamic yet
      ?_x_salome_session.geompy = _x_salome_session.geomBuilder.New()
      may be better never change initial study? ... wait and see
      version SALOME 9: only one study
      """
      import salome

      tmpStudy = salome.myStudy
      if tmpStudy != None:
        tmpName = tmpStudy._get_Name()
      else:
        tmpName = "None"

      print("dir(salome)\n%s" % PP.pformat(dir(salome)))
      print("dir(tmpStudy)\n%s" % PP.pformat(dir(tmpStudy)))

      print("activeStudy to close:", tmpName)  # ,tmpStudy #mefiance: salome.myStudy sometimes is not renewed

      newStudy = stdmgr.NewStudy(aName)
      salome.setCurrentStudy(newStudy)
      # print "createNewStudy", newStudy._get_Name(), tmpStudy, newStudy, salome.myStudy
      if aName != salome.myStudyName:
        raise Exception("unexpected study name: %s expected:%s" % (salome.myStudyName, aName))
      print("createNewStudy", newStudy._get_Name())
      return newStudy

      # if tmpStudy != None:
      #    #print "activeStudy.IsEmpty", tmpStudy.IsEmpty()
      #    tmpStudy.Close()
      # aInt = salome.createNewStudy()
      # print "createNewStudy", aInt, salome.myStudy

      """stdmgr = salome.myStudyManager
      studies = stdmgr.GetOpenStudies()
      if len(studies) == 1:
          activeStudy = stdmgr.GetStudyByName(studies[0])
      stdmgr.Close(activeStudy)
      #bug in salome.py: salome.myStudy is not renewed

      _nbStudies += 1
      newStudy = stdmgr.NewStudy("newStudy%s" % _nbStudies)
      print "createNewStudy", newStudy._get_Name(), tmpStudy, newStudy
      salome.myStudy = newStudy #do it for bug?
      #reste au GUI a faire Menu-File-Connect: n'est pas present dans l'API python?
      return newStudy
      #salome.myStudy.attach()
      #return salome.myStudy
      #import salome_study
      #myStudyManager, myStudyId, myStudy, myStudyName = salome_study.salome_study_init()
      #print "myStudy",myStudy,myStudyId,myStudyName"""

    def createNewStudy_old(self, name=None):
      """
      called to clean active study: close old study before
      return the SALOMEDS._objref_Study instance
      risky because users of geompy etc...
      have to "renewed their geompy" if they are in global
      geompy = XSS.getGeompy() is not dynamic yet
      ?_x_salome_session.geompy = _x_salome_session.geomBuilder.New()
      may be better never change initial study? ... wait and see
      version SALOME < 9: may be some studies
      """
      global _nbStudies
      import salome
      tmpStudy = salome.myStudy
      if tmpStudy != None:
        tmpName = tmpStudy._get_Name()
      else:
        tmpName = "None"
      print("activeStudy to close:", tmpName)  # ,tmpStudy #mefiance: salome.myStudy sometimes is not renewed
      stdmgr = salome.myStudyManager
      studies = stdmgr.GetOpenStudies()
      # print "stdmgr.GetOpenStudies()",studies
      activeStudy = None
      if len(studies) == 1:
        activeStudy = stdmgr.GetStudyByName(studies[0])
      # for i in dir(activeStudy): print "activeStudy", i
      # print "before stdmgr.Close(activeStudy) IsEmpty", activeStudy.IsEmpty()
      if activeStudy != None:
        stdmgr.Close(activeStudy)
      salome.myStudy, salome.myStudyName = None, None
      _nbStudies += 1
      stdmgr = salome.myStudyManager
      if name == None:
        aName = "newStudy%s" % _nbStudies
      else:
        aName = name
      newStudy = stdmgr.NewStudy(aName)
      salome.setCurrentStudy(newStudy)
      # print "createNewStudy", newStudy._get_Name(), tmpStudy, newStudy, salome.myStudy
      if aName != salome.myStudyName:
        raise Exception("inexpected study name: %s expected:%s" % (salome.myStudyName, aName))
      print("createNewStudy", newStudy._get_Name())
      return newStudy

      # if tmpStudy != None:
      #    #print "activeStudy.IsEmpty", tmpStudy.IsEmpty()
      #    tmpStudy.Close()
      # aInt = salome.createNewStudy()
      # print "createNewStudy", aInt, salome.myStudy

      """stdmgr = salome.myStudyManager
      studies = stdmgr.GetOpenStudies()
      if len(studies) == 1:
          activeStudy = stdmgr.GetStudyByName(studies[0])
      stdmgr.Close(activeStudy)
      #bug in salome.py: salome.myStudy is not renewed

      _nbStudies += 1
      newStudy = stdmgr.NewStudy("newStudy%s" % _nbStudies)
      print "createNewStudy", newStudy._get_Name(), tmpStudy, newStudy
      salome.myStudy = newStudy #do it for bug?
      #reste au GUI a faire Menu-File-Connect: n'est pas present dans l'API python?
      return newStudy
      #salome.myStudy.attach()
      #return salome.myStudy
      #import salome_study
      #myStudyManager, myStudyId, myStudy, myStudyName = salome_study.salome_study_init()
      #print "myStudy",myStudy,myStudyId,myStudyName"""

    def cleanStudy_risky(self, verbose=None):
        """
        fait le menage des geom dans occ viewer
        
        http://www.salome-platform.org/forum/forum_10/366900504Pour un objet shape créé avec geompy, il faut utiliser:
        et Christophe..
        sb = salome.myStudy.NewBuilder()
        so = salome.ObjectToSObject(shape)
        sb.RemoveObjectWithChildren( so )
        shape.UnRegister()

        Pour un objet mesh créé avec smesh.Mesh(), il faut utiliser:
        meshSO = salome.ObjectToSObject(mesh.GetMesh())
        sb.RemoveObjectWithChildren( meshSO )

        L'appel à UnRegister() est fait dans RemoveObjectWithChildren pour les objets Mesh.
        Attention, la destruction des objets peut être longue. Il faut évaluer si le gain de mémoire obtenu vaut le coup de perdre du temps à détruire les objets.
        """
        """
        from geomtools.py deleteShape(self,shapeStudyEntry):
        This completly deletes a geom shape.
        WARNING: please be aware that to delete a geom object, you have
        three operations to perform:
        1. erase the shape from the viewers
        2. remove the entry from the study
        3. destroy the underlying geom object
        self.eraseShapeByEntry(shapeStudyEntry)
        shape = self.removeFromStudy(shapeStudyEntry)
        shape.Destroy()
        """
        def findGEOM_recursive(Study, SO, Builder, offset, studyEntries):
            import GEOM
            from salome.kernel.services import IDToObject, IDToSObject
            it = Study.NewChildIterator(SO)
            while it.More():
                CSO = it.Value()
                entry = CSO.GetID()
                a = offset * "--"  + " GEOM "+ CSO.GetID()
                find, AtName = Builder.FindAttribute(CSO, "AttributeName")
                if find:
                  a = a + ":" + AtName.Value()
                #find,AtIOR = Builder.FindAttribute(CSO, "AttributeIOR")
                #if find:
                #  a=a+":"+AtIOR.Value()
                #find,RefSO = CSO.ReferencedObject()
                #if find:
                #  a=a+":"+RefSO.GetID()
                ###print a
                findGEOM_recursive(Study, CSO, Builder, offset+2, studyEntries)
                geomObject=IDToObject(entry, Study)
                if geomObject != None:
                  shape = geomObject._narrow(GEOM.GEOM_Object)
                  studyEntries.append(("GEOM", entry, CSO, geomObject, shape))
                else:
                  #print 'study entry is empty:',a
                  studyEntries.append(("ALL", entry, CSO, None, None))
                it.Next()

        def findSMESH_recursive(Study, SO, Builder, offset, studyEntries):
            import SMESH
            from salome.kernel.services import IDToObject, IDToSObject
            it = Study.NewChildIterator(SO)
            while it.More():
                CSO = it.Value()
                entry = CSO.GetID()
                a = offset * "--" + " SMESH " + CSO.GetID()
                find, AtName = Builder.FindAttribute(CSO, "AttributeName")
                if find:
                  a = a + ":" + AtName.Value()
                #find,AtIOR = Builder.FindAttribute(CSO, "AttributeIOR")
                #if find:
                #  a=a+":"+AtIOR.Value()
                #find,RefSO = CSO.ReferencedObject()
                #if find:
                #  a=a+":"+RefSO.GetID()
                ###print a
                #findSMESH_recursive(Study, CSO, Builder, offset+2, studyEntries)
                geomObject=IDToObject(entry, Study)
                if geomObject != None:
                  shape = None #geomObject._narrow(GEOM.GEOM_Object)
                  studyEntries.append(("SMESH", entry, CSO, geomObject, shape))
                else:
                  #print 'study entry is empty:',a
                  studyEntries.append(("ALL", entry, CSO, None, None))
                it.Next()
        
        print("xsalomesession.py cleanStudy_risky")
        import salome
        study = salome.myStudy
        itcomp = study.NewComponentIterator()
        builder = study.NewBuilder()
        studyEntries = []
        while itcomp.More():
            SC = itcomp.Value()
            name = SC.ComponentDataType()
            if verbose: print("Clean ComponentDataType %s" % name)
            #builder.RemoveObjectWithChildren( SC ) #ne marche pas bien... plante apres...
            it = study.NewChildIterator(SC)
            offset = 1
            if name == "GEOM": findGEOM_recursive(study, SC, builder, offset, studyEntries)
            if name == "SMESH": findSMESH_recursive(study, SC, builder, offset, studyEntries)
            itcomp.Next()

        from salome.kernel.studyedit import getStudyEditor
        studyEditor = getStudyEditor()
        #from salome.geom import geomtools  ...not use...buggué?..
        #gst = geomtools.GeomStudyTools(studyEditor)

        geomGui = self.getGeomGUI()
        for name, entry, CSO, geomObject, shape in studyEntries:
          
          if name == "GEOM":
            if verbose: print("cleanStudy_risky GEOM Remove %s %s" % (name,entry))
            if geomGui != None:
              try:
                eraseFromAllWindows=True
                geomGui.eraseGO(entry, eraseFromAllWindows)  # if GUI only
              except:
                # print("cleanStudy_risky no geomGUI.eraseGO")
                # print("xsalomesession.py dir(geomGUI):\n%s" % PP.pformat(dir(geomGUI)))
                # print("geomGUI", geomGUI)
                pass

            #shape.UnRegister()
            builder.RemoveObjectWithChildren( CSO )
            ##withChildren = True
            ##studyEditor.removeItem(CSO, withChildren)
            
          if name == "SMESH":
            print("TODO cleanStudy_risky SMESH Remove", name, entry)

          if name == "ALL":
            if verbose: print("Remove ", name, entry, CSO)
            builder.RemoveObjectWithChildren( CSO )

        """from salome.smesh import smeshBuilder
        from salome.geom import geomBuilder
        #self.geomBuilder = geomBuilder
        self.geompy = geomBuilder.New()
        #self.smeshBuilder = smeshBuilder
        self.smeshpy = smeshBuilder.New()"""

        updateObjBrowser()
        return study.IsEmpty()

    def cleanStudy_risky2(self, verbose=None):
        """
        ne fait PAS le menage des geom dans occ viewer
        
        http://www.salome-platform.org/forum/forum_10/366900504Pour un objet shape créé avec geompy, il faut utiliser:
        et Christophe..
        sb = salome.myStudy.NewBuilder()
        so = salome.ObjectToSObject(shape)
        sb.RemoveObjectWithChildren( so )
        shape.UnRegister()

        Pour un objet mesh créé avec smesh.Mesh(), il faut utiliser:
        meshSO = salome.ObjectToSObject(mesh.GetMesh())
        sb.RemoveObjectWithChildren( meshSO )

        L'appel à UnRegister() est fait dans RemoveObjectWithChildren pour les objets Mesh.
        Attention, la destruction des objets peut être longue. Il faut évaluer si le gain de mémoire obtenu vaut le coup de perdre du temps à détruire les objets.
        """
        def find_recursive(Study, SO, Builder, offset, studyEntries, name):
            import GEOM
            from salome.kernel.services import IDToObject, IDToSObject
            it = Study.NewChildIterator(SO)
            while it.More():
                CSO = it.Value()
                entry = CSO.GetID()
                a = offset * "--" + name + " " + CSO.GetID()
                find, AtName = Builder.FindAttribute(CSO, "AttributeName")
                if find:
                  a = a + ":" + AtName.Value()
                #find,AtIOR = Builder.FindAttribute(CSO, "AttributeIOR")
                #if find:
                #  a=a+":"+AtIOR.Value()
                #find,RefSO = CSO.ReferencedObject()
                #if find:
                #  a=a+":"+RefSO.GetID()
                print(a)
                #only first level!!! find_recursive(Study, CSO, Builder, offset+2, studyEntries, name)
                studyEntries.append((name, entry, CSO))
                it.Next()

        import salome
        study = salome.myStudy
        itcomp = study.NewComponentIterator()
        builder = study.NewBuilder()
        studyEntries = []
        while itcomp.More():
            SC = itcomp.Value()
            name = SC.ComponentDataType()
            if verbose: print("Clean ComponentDataType %s" % name)
            #builder.RemoveObjectWithChildren( SC ) #ne marche pas bien... plante apres...
            it = study.NewChildIterator(SC)
            offset = 1
            find_recursive(study, SC, builder, offset, studyEntries, name)
            itcomp.Next()

        sb = salome.myStudy.NewBuilder()
        for name, entry, CSO in studyEntries:
            sb.RemoveObjectWithChildren( CSO )
        
        #study.UpdateIORLabelMap() #not swigged
        updateObjBrowser()
        return study.IsEmpty()

    def kill(self):
        print("\nXSalomeSession.kill port %s" % self.port)
        global _x_salome_session
        _x_salome_session = None
        if self.port != None:
            try:
                os.system('killSalomeWithPort.py %s' % self.port)
            except Exception as e:
                print("WARNING: XSalomeSession.kill %s" % e)
            return
        # print "XSalomeSession.kill without port, nothing to do."
        return

    def __del__(self):
        if verbose: print("XSalomeSession.__del__")
        # avoid th. SALOME_NamingService.cxx [1379] : Destroy_Directory(): CosNaming::NamingContext::NoEmpty /Study/ is not empty
        # closeActiveStudy()
        # Exception ImportError: 'No module named salome' in <bound method XSalomeSession.__del__ of <salomepy.xsalomesession.XSalomeSession object at 0x1083ed0>> ignored
        # self.kill() # seem tricky as deleted object, may be explicitely user call


"""
for EZ example of use:

import salomepy.xsalomesession as XSS
geompy = XSS.getGeompy()
smesh = XSS.getSmesh()
geomBuilder = XSS.getGeomBuilder()
smeshBuilder = XSS.getSmeshBuilder()
geomGUI = XSS.getGeomGUI() #None if not GUI mode
desktop = XSS.getDesktop()
newStudy = XSS.createNewStudy()

import salomepy.xsalomesession as XSS
newStudy = XSS.createNewStudy()
"""

def killXSalomeSession():
    global _x_salome_session
    print("killXSalomeSession %s" % _x_salome_session)
    if _x_salome_session != None:
      # erase all objects from the current view from previous test in same salome session
      import salome
      salome.sg.EraseAll()
      _x_salome_session.kill()


def isSalomeAlreadyLaunched():
    in_salome = False
    """
    ## Bug dans Salome 6.6.0 installé par les wizard, NSparam renvoie 2810
    ## même si Salome n'a pas été lancé
    import os
    ## Si $OMNIORB_CONFIG n'est pas dans $HOME, il n'a pas été fixé au lancement de Salome
    omniorb_cfg = os.environ["OMNIORB_CONFIG"]
    omniorb_cfg_dir = os.path.dirname(omniorb_cfg)
    home = os.environ["HOME"]
    if omniorb_cfg_dir != home:
        ## => On doit lancer Salome
        return False
    ## fin contournement bug
    """
    kernel = os.getenv("KERNEL_ROOT_DIR")
    if kernel == None: 
        raise Exception("No salome environment set.")
    try:
        import NSparam
    except:
        raise Exception("Impossible to import NSparam.")
    try:
        machine, port = NSparam.getNSparams()
        if port:
            in_salome = True
            if verbose: print("isSalomeAlreadyLaunched : %s port %i" % (in_salome, port))
            return True
    except Exception as e:
        print("problem isSalomeAlreadyLaunched NSparam.getNSparams %s" % e)
        pass
    # print("!!!! DEBUG TODO isSalomeAlreadyLaunched forced to False !!!")
    import salome as SS
    in_salome = SS.hasDesktop()
    if not(in_salome): print("isSalomeAlreadyLaunched forced to False as hasDesktop(), fix it to another trick")
    return in_salome

def closeActiveStudy():
    import salome
    tmpStudy = salome.myStudy

    print("activeStudy to close:", tmpStudy._get_Name()) #,tmpStudy #mefiance: salome.myStudy is not renewed
    print("dir(salome)\n%s" % PP.pformat(dir(salome)))
    print("dir(tmpStudy)\n%s" % PP.pformat(dir(tmpStudy)))

    stdmgr = salome.myStudyManager
    studies = stdmgr.GetOpenStudies()
    if len(studies) == 1: #only one could exist
        activeStudy = stdmgr.GetStudyByName(studies[0])
    else:
        raise Exception("inexpected: more than one study")
    stdmgr.Close(activeStudy)
    salome.myStudy, salome.myStudyName = None, None

def getXSalomeSession(*args, **kwargs):
    global _x_salome_session
    global _nbTriesLaunchSalome
    #print "getXSalomeSession",_x_salome_session
    if _x_salome_session != None:
        if _x_salome_session.needRelaunch():
            _nbTriesLaunchSalome +=1
            if _nbTriesLaunchSalome < 2:
              _x_salome_session.launchSalome()
            else:
              print("\n!!! too much tries launchSalome, abort !!!\n")
              self.kill()
              return None
        return _x_salome_session
    _x_salome_session = XSalomeSession(*args, **kwargs)
    return getXSalomeSession(*args, **kwargs)

def getSg():
    getXSalomeSession()
    return _x_salome_session.getGeompy()

def getGeompy():
    getXSalomeSession()
    #import salome
    #have to "renewed their geompy" if they are in global
    #_x_salome_session.geompy = _x_salome_session.geomBuilder.New()
    return _x_salome_session.getGeompy()

def getSmesh():
    getXSalomeSession()
    return _x_salome_session.getSmesh()

def getSmeshBuilder():
    getXSalomeSession()
    return _x_salome_session.getSmeshBuilder()

def getGeomBuilder():
    getXSalomeSession()
    return _x_salome_session.getGeomBuilder()

def getGeomGUI():
    getXSalomeSession()
    return _x_salome_session.getGeomGUI()

def getSmeshGUI():
  getXSalomeSession()
  return _x_salome_session.getSmeshGUI()

def getDesktop():
    if isSalomeAlreadyLaunched():
      getXSalomeSession()
      return _x_salome_session.getDesktop()
    return None

def createNewStudy():
    getXSalomeSession()
    return _x_salome_session.createNewStudy()

def getActiveStudy():
    try:
      import salome
      return salome.myStudy
    except:
      mess = "getActiveStudy can't import salome"
      if mess not in _messDone:
        logger.warning(mess)
        _messDone.append(mess)
      return None
      
def openHdf(name):
    try:
      import salome
    except:
      mess = "openHdf can't import salome"
      if mess not in _messDone:
        logger.warning(mess)
        _messDone.append(mess)
      return (False, "ERROR: openHdf: can't import salome")

    if not os.path.isfile(name):
      return (False, "ERROR: openHdf: inexisting file: %s" % name)
    try:
      closeActiveStudy()
      study = salome.myStudyManager.Open(name)
      salome.setCurrentStudy(study)
      return (True, "")
    except:
      return (False, "ERROR: openHdf: problem on file: %s" % name)

def cleanStudy_risky():
    getXSalomeSession()
    return _x_salome_session.cleanStudy_risky()

def cleanStudy_risky2():
    getXSalomeSession()
    return _x_salome_session.cleanStudy_risky2()

def updateObjBrowser():
    try:
      import salome
    except:
      mess = "updateObjBrowser can't import salome"
      if mess not in _messDone:
        logger.warning(mess)
        _messDone.append(mess)
      return
    if salome.sg.hasDesktop():
      salome.sg.updateObjBrowser()

def getSelectedCount():
    getXSalomeSession()
    desktop = getDesktop()
    if desktop != None:
      count = _x_salome_session.sg.SelectedCount()
      return count
    else:
      print("WARNING: getSelectedCount: no desktop: no selection")
      return None

def getAllSelected():
    getXSalomeSession()
    selNumber = getSelectedCount()
    if selNumber == None: return None #no desktop
    listSelected = []
    for i in range(selNumber):
        listSelected.append(_x_salome_session.sg.getSelected(i))
    return listSelected

def getSelected(i=0):
    """get first selected by default if multiple"""
    getXSalomeSession()
    selNumber = getSelectedCount()
    if selNumber == None: return None #no desktop
    return _x_salome_session.sg.getSelected(i)

def getSelectedGeom(Verbose=False, Load=False):
    from salome.kernel.services import IDToObject, IDToSObject
    import GEOM
    entries = getAllSelected()
    study = getActiveStudy()
    res = []
    for entry in entries:
      sobj = IDToSObject( entry, study)
      #obj = IDToObject( entry, study)
      #print "entry GEOM",entry,sobj,sobj.GetName(),obj
      if ( sobj ) != None:
        # check if it is a geom object
        ComponentFather = sobj.GetFatherComponent()
        typeComponent = ComponentFather.ComponentDataType()
        if str(typeComponent) == 'GEOM' :
           # retrieve geom object
           obj = sobj.GetObject()
           if Load and obj==None: obj = getOrLoadSobj(sob)
           if Verbose: print("geom object",sobj,obj) #,obj._narrow(GEOM.GEOM_Object)
           res.append(sobj)
    return res

def getSelectedSmesh(Verbose=False, Load=False):
    from salome.kernel.services import IDToObject, IDToSObject
    entries = getAllSelected()
    study = getActiveStudy()
    res = []
    for entry in entries:
      sobj = study.FindObjectID( entry )
      #obj = IDToObject( entry, study)
      #print "entry SMESH",entry,sobj,obj
      if ( sobj ) != None:
        # check if it is a smesh object
        ComponentFather = sobj.GetFatherComponent()
        typeComponent = ComponentFather.ComponentDataType()
        if str(typeComponent) == 'SMESH' :
           # retrieve smesh object
           obj = sobj.GetObject()
           if Load and obj==None: obj = getOrLoadSobj(sob)
           if Verbose: print("smesh object",sobj,obj)
           res.append(sobj)
    return res

def getGeomSobj(Verbose=False):
    """
    import salomepy.xsalomesession as XSS
    gg = XSS.getGeomSobj(True)
    ss = XSS.getSmeshSobj(True)
    """
    try:
      import salome
    except:
      mess = "getGeomSobj can't import salome"
      if mess not in _messDone:
        logger.warning(mess)
        _messDone.append(mess)
      return []
    sgeom = salome.myStudy.FindObject("Geometry")
    iter = salome.myStudy.NewChildIterator(sgeom) # initialize from the component
    iter.InitEx(True) # init recursive mode
    res = []
    while iter.More():
      c = iter.Value()
      if Verbose: print("Geometry",c, c.GetID(),c.GetObject()) #,getOrLoadObjectFromSObject(c)
      res.append(c)
      iter.Next()
    return res

def getSmeshSobj(Verbose=False):
    try:
      import salome
    except:
      mess = "getSmeshSobj can't import salome"
      if mess not in _messDone:
        logger.warning(mess)
        _messDone.append(mess)
      return []
    smesh = salome.myStudy.FindObject("Mesh")
    iter = salome.myStudy.NewChildIterator(smesh) # initialize from the component
    iter.InitEx(True) # init recursive mode
    res = []
    while iter.More():
      c = iter.Value()
      if Verbose: print("Mesh", c, c.GetID(),c.GetObject()) #,getOrLoadObjectFromSObject(c)
      res.append(c)
      iter.Next()
    return res

def getOrLoadSobj(sobjs):
    try:
      import salome
    except:
      mess = "getOrLoadSobj can't import salome"
      if mess not in _messDone:
        logger.warning(mess)
        _messDone.append(mess)
      return None
    from salome.kernel.studyedit import getStudyEditor
    studyEditor = getStudyEditor()
    if type(sobjs) == list:
      res = [studyEditor.getOrLoadObject(s) for s in sobjs]
    else:
      res = studyEditor.getOrLoadObject(sobjs)
    return res

def getSobjPath(sobjs):
    try:
      import salome
    except:
      mess = "getSobjPath can't import salome"
      if mess not in _messDone:
        logger.warning(mess)
        _messDone.append(mess)
      return None
    study =salome.myStudy
    if type(sobjs) == list:
      res = [study.GetObjectPath(s) for s in sobjs]
    else:
      res = study.GetObjectPath(sobjs)
    return res

def getOrLoadGeomPath(aPath):
    res = getGeomSobj()
    if res == None:
      mess = "getOrLoadGeomPath can't get Geom"
      if mess not in _messDone:
        logger.warning(mess)
        _messDone.append(mess)
      return None
    items = getSobjPath(res)
    if items == None:
      mess = "getOrLoadGeomPath no Geom items"
      if mess not in _messDone:
        logger.warning(mess)
        #_messDone.append(mess)
      return None
    items = [item.replace("/Geometry/","") for item in items]
    resOk = []
    ii=0
    if verbose: print("getOrLoadGeomPath",items)
    for item in items:
      if aPath in item:
        resOk.append(res[ii])
        getOrLoadSobj(res[ii])
      ii += 1
    return resOk

def saveStudy(aPath):
    try:
      import salome
    except:
      mess = "saveStudy can't import salome: no save"
      if mess not in _messDone:
        logger.error(mess)
        #_messDone.append(mess)
      return False  # ko
    tmpStudy = salome.myStudy
    if tmpStudy == None:
      mess = "saveStudy no Study to save"
      if mess not in _messDone:
        logger.warning(mess)
        #_messDone.append(mess)
      return False  # ko

    # print("dir(salome)\n%s" % PP.pformat(dir(salome)))
    # print("dir(tmpStudy)\n%s" % PP.pformat(dir(tmpStudy)))
    # help(tmpStudy)
    # help(tmpStudy.SaveAs)

    aName = aPath
    if aPath[-4:] != ".hdf": aName += ".hdf"
    if verbose: print("saveStudy: '%s'" % aName)

    tmpStudy.SaveAs(aName, False, False)  # name, bool theMultiFile, bool theASCII
    return True  # ok


def exportXAOGeomPath(aGeomPath, aDir, creator=None):
  """
  export to some files .xao some geometries in directory aDir
  to export all grainxxx set aGeomPath to grain
  """
  geompy = getGeompy()
  res = getOrLoadGeomPath(aGeomPath)
  savedFiles = ""
  if creator == None:
    aCreator = os.getenv("USERNAME")
  else:
    aCreator = str(creator)
  for s in res:
    aName = s.GetName()
    saveFile = os.path.join(aDir,aName + ".xao")
    #GEOM_Object shape, ListOfGO groups, ListOfFields fields, string author, string fileName
    obj = s.GetObject()
    groups = geompy.GetGroups(obj)
    fields = geompy.GetFields(obj)
    ok = geompy.ExportXAO(obj, groups, fields, aCreator, saveFile)
    savedFiles += saveFile + "\n"
  return savedFiles

def ObjectToSObject(obj):
  import salome
  return salome.ObjectToSObject(obj)


"""
sandbox...

def loadComponentGeom():
    import salome
    c = salome.lcc.FindOrLoadComponent("FactoryServer", "GEOM")
    print "ComponentGEOM",c
    #c.SetCurrentStudy( salome.myStudy )
    #c.ActivateModule()
    return c

def loadComponentSmesh():
    import salome
    c = salome.lcc.FindOrLoadComponent("FactoryServer", "SMESH")
    print "ComponentSMESH",c
    #c.SetCurrentStudy( salome.myStudy )
    #c.ActivateModule()
    return c


import salomepy.xsalomesession as XSS
import salome
c = salome.lcc.FindOrLoadComponent("FactoryServer", "GEOM")
print "ComponentGEOM",c
ComponentGEOM <salome.geomBuilder object at 0x1bbbaf0>
dir(c)

c = salome.lcc.FindOrLoadComponent("FactoryServer", "SMESH")
dir(c)

XSS.getGeomGUI()
<libGEOM_Swig.GEOM_Swig; proxy of <Swig Object of type 'GEOM_Swig *' at 0x25148c0> >

myStudy=salome.myStudy
myBuilder = myStudy.NewBuilder()
mySComp = myStudy.FindComponent("SMESH")
print mySComp
<SALOMEDS._objref_SComponent instance at 0x339fb170>
myBuilder.LoadWith(mySComp, c)

"""
