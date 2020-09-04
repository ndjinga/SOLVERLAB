from PyQt4 import QtGui, QtCore
from PyQt4.uic import loadUi
from utils import completeResPath

import CoreFlows as cf
import cdmath as cm

class MainCFWidget(QtGui.QTabWidget):
  def __init__(self):
    QtGui.QTabWidget.__init__(self)
    loadUi(completeResPath("MainCFWidget.ui"), self)
    self._python_dump = []
    
  def scanWidgets(self):
    print self.tabModel
    dictCF={}
    for k in self.__dict__:
      att = self.__dict__[k] 
      if isinstance(att, QtGui.QWidget):
        name = str(att.objectName())
        if name != "":
#          print(name)
          if name.endswith("RadioButton"):
            assert(isinstance(att, QtGui.QRadioButton))
            if att.isChecked() :
              # parse name
              name=name[:len(name)-len("_RadioButton")]#On retire le suffixe _Radiobutton
              if name.startswith("Dim") :
                dictCF["spaceDim"]=int(name[len("Dim")])
              elif name=="DriftModel" or "name==SinglePhase" or name.endswith("Equation") or name.endswith("TwoFluid") :
                dictCF["ModelName"]=name
              elif name=="InputFileName" :
                dictCF["InputFileName"]=True
              elif name=="MeshCreation" :
                dictCF["MeshCreation"]=True
          elif name.endswith("spinBox") :
            assert(isinstance(att, QtGui.QSpinBox))
            val = att.value()
            # parse name
            name=name[:len(name)-len("_spinBox")]#On retire le suffixe _SpinBox
            if name=="Nx" :
              dictCF["Nx"]=val
            elif name=="Ny" :
              dictCF["Ny"]=val
            elif name=="Nz" :
              dictCF["Nz"]=val
            elif name=="NO_MaxNbOfTimeStep" :
              dictCF["MaxNbOfTimeStep"]=val
            elif name=="NO_FreqSave" :
              dictCF["FreqSave"]=val
          elif name.endswith("doubleSpinBox"):
            assert(isinstance(att, QtGui.QDoubleSpinBox))
            val = att.value()
            # parse name
            name=name[:len(name)-len("_doubleSpinBox")]#On retire le suffixe _doubleSpinBox
            if name.endswith("Conductivity") :
              dictCF["Conductivity"]=val
            elif name.endswith("Conductivity_Phase1") :
              dictCF["Conductivity1"]=val
            elif name.endswith("Conductivity_Phase2") :
              dictCF["Conductivity2"]=val
            elif name.endswith("Viscosity") :
              dictCF["Viscosity"]=val
            elif name.endswith("Viscosity_Phase1") :
              dictCF["Viscosity1"]=val
            elif name.endswith("Viscosity_Phase2") :
              dictCF["Viscosity2"]=val
            elif name.endswith("HeatSource") :
              dictCF["HeatSource"]=val
            elif name.endswith("FrictionCoefficients") :
              dictCF["FrictionCoefficients"]=val
            elif name.endswith("Gravity_1d") :
              dictCF["Gravity_1d"]=val
            elif name.endswith("Gravity_2d") :
              dictCF["Gravity_2d"]=val
            elif name.endswith("Gravity_3d") :
              dictCF["Gravity_3d"]=val
            elif name=="Xinf" :
              dictCF["Xinf"]=val
            elif name=="Xsup" :
              dictCF["Xsup"]=val
            elif name=="Yinf" :
              dictCF["Yinf"]=val
            elif name=="Ysup" :
              dictCF["Ysup"]=val
            elif name=="Zinf" :
              dictCF["Zinf"]=val
            elif name=="Zsup" :
              dictCF["Zsup"]=val
            elif name=="NO_TimeMax" :
              dictCF["TimeMax"]=val
            elif name=="NO_Precision" :
              dictCF["Precision"]=val
            elif name=="NO_CFL" :
              dictCF["CFL"]=val
            elif name.endswith("_IC") :#Initial conditions
              name=name[:len(name)-len("_IC")]#On retire le suffixe _IC
              if name.endswith("Concentration") :
                dictCF["InitialConcentration"]=val
              elif name.endswith("Alpha") :
                dictCF["InitialAlpha"]=val
              elif name.endswith("Pressure") :
                dictCF["InitialPressure"]=val
              elif name.endswith("Temperature") :
                dictCF["InitialTemperature"]=val
              elif name.endswith("Velocity_1d") :
                dictCF["InitialVelocity_1d"]=val
              elif name.endswith("Velocity_2d") :
                dictCF["InitialVelocity_2d"]=val
              elif name.endswith("Velocity_3d") :
                dictCF["InitialVelocity_3d"]=val
            elif name.endswith("_BC") :#Boundary conditions
              name=name[3:len(name)-len("_BC")]#On retire le préfixe SP_ ou DM_ au début et le suffixe _BC à la fin
              dictCF[name]=val
          elif name.endswith('comboBox'):
            assert(isinstance(att, QtGui.QComboBox))
            val = att.currentText()
            # parse name
            name=name[:len(name)-len("_comboBox")]#On retire le suffixe _comboBox
            if name.endswith("1barOr155bar") :
              dictCF["pressureEstimate"]=val
            elif name.endswith("GasOrLiquid") :
              dictCF["fluidType"]=val
            elif name.endswith("_BC") :#Boundary conditions
              name=name[4:len(name)-len("_BC")]#On retire le préfixe SP_ ou DM_ au début et le suffixe _BC à la fin
              dictCF[name]=val
            elif name=="NO_Method" :
              dictCF["Method"]=val
            elif name=="NO_LS" :
              dictCF["LS"]=val #Linear solver
            elif name=="NO_Scheme" :
              dictCF["Scheme"]=val 
            elif name=="NO_Scheme_type" :
              dictCF["Scheme_type"]=val 
            elif name=="NO_Preconditioner" :
              dictCF["Preconditioner"]=val 
          elif name.endswith('lineEdit'):
            assert(isinstance(att, QtGui.QLineEdit))
            val = str(att.text())
            # parse name
            name=name[:len(name)-len("_lineEdit")]#On retire le suffixe _comboBox
            #print(name,val)
            if name=="NO_ResultFileName" :
              dictCF["ResultFileName"]=val
            elif name=="InputFileName" :
              dictCF["InputFileName"]=val

    return dictCF  
          
  def onLaunchSimu(self):
    print "Reading widgets"
    dictCF = self.scanWidgets()

    print "Setting Model and VV_Constant"
    ######## Setting Model and VV_Constant #########################
    if dictCF["ModelName"]=="SinglePhase" :
      exec "myProblem = cf.%s(cf.%s,cf.%s,%s)" % (dictCF["ModelName"],dictCF["fluidType"],dictCF["pressureEstimate"],dictCF["spaceDim"])
      nVar =  myProblem.getNumberOfVariables()
      VV_Constant =[0]*nVar
      VV_Constant[0] = dictCF["InitialPressure"]
      VV_Constant[1] = dictCF["InitialVelocity_1d"]
      if dictCF["spaceDim"] >1 :
        VV_Constant[2] = dictCF["InitialVelocity_2d"]
        if dictCF["spaceDim"] >2 :
          VV_Constant[3] = dictCF["InitialVelocity_3d"]
      VV_Constant[nVar-1] = dictCF["InitialTemperature"]
    elif dictCF["ModelName"]=="DriftModel" :
      exec "myProblem = cf.%s(cf.%s,%s)" % (dictCF["ModelName"],dictCF["pressureEstimate"],dictCF["spaceDim"])
      nVar =  myProblem.getNumberOfVariables()
      VV_Constant =[0]*nVar
      VV_Constant[0] = dictCF["InitialConcentration"]
      VV_Constant[1] = dictCF["InitialPressure"]
      VV_Constant[2] = dictCF["InitialVelocity_1d"]
      if dictCF["spaceDim"] >1 :
        VV_Constant[3] = dictCF["InitialVelocity_2d"]
        if dictCF["spaceDim"] >2 :
          VV_Constant[4] = dictCF["InitialVelocity_3d"]
      VV_Constant[nVar-1] = dictCF["InitialTemperature"]
    else :
      raise NameError('Model not yet handled', dictCF["ModelName"])  

    print "Setting initial data"
    ############ setting initial data ################################
    if dictCF["spaceDim"] ==1 :
      print "spaceDim= ", dictCF["spaceDim"]
      print "VV_Constant= ", VV_Constant
      print "Xinf= ", dictCF["Xinf"]
      print "Xsup= ", dictCF["Xsup"]
      print "Nx= ", dictCF["Nx"]
      myProblem.setInitialFieldConstant( dictCF["spaceDim"], VV_Constant, dictCF["Xinf"], dictCF["Xsup"], dictCF["Nx"],"Left","Right");
      print "Initial field set"
    elif dictCF["spaceDim"] ==2 :
      myProblem.setInitialFieldConstant( dictCF["spaceDim"], VV_Constant, dictCF["Xinf"], dictCF["Xsup"], dictCF["Nx"],"Left","Right", dictCF["Yinf"], dictCF["Ysup"], dictCF["Ny"],"Bottom","Top");
    elif dictCF["spaceDim"] ==3 :
      myProblem.setInitialFieldConstant( dictCF["spaceDim"], VV_Constant, dictCF["Xinf"], dictCF["Xsup"], dictCF["Nx"],"Back","Front", dictCF["Yinf"], dictCF["Ysup"], dictCF["Ny"],"Left","Right", dictCF["Zinf"], dictCF["Zsup"], dictCF["Nz"],"Bottom","Top");
    else :
      raise NameError('Dimension should be 1, 2 or 3', dictCF["spaceDim"])  

#    for k, v in dictCF.items():
#        line = "cf.set%s(%s)" % (k, v)
#        #exec line
#        self._python_dump.append(line)

    ######## 1D for the moment ######################
    if dictCF["ModelName"]=="SinglePhase" :
      myProblem.setInletBoundaryCondition("Left",dictCF["Temperature_Left"],dictCF["Concentration_Left"],dictCF["Velocity_1d_Left"])
      myProblem.setOutletBoundaryCondition("Right", dictCF["Pressure_Right"]);
    elif dictCF["ModelName"]=="DriftModel" :
      myProblem.setInletBoundaryCondition("Left",dictCF["DM_Temperature_Left"],dictCF["DM_Velocity_1d_Left"])
      myProblem.setOutletBoundaryCondition("Right", dictCF["DM_Pressure_Right"]);

    ########## Physical forces #################
    myProblem.setHeatSource(dictCF["HeatSource"]);
    gravite=[0]*dictCF["spaceDim"]
    gravite[0]=dictCF["Gravity_1d"]
    if dictCF["spaceDim"] >1 :
      gravite[1]=dictCF["Gravity_2d"]
      if dictCF["spaceDim"] >2 :
        gravite[2]=dictCF["Gravity_3d"]
    myProblem.setGravity(gravite)

    ########## Numerical options ###############
    exec "myProblem.setNumericalScheme(cf.%s, cf.%s)" % (dictCF["Scheme"],dictCF["Method"])
    myProblem.setWellBalancedCorrection(True);  

    myProblem.setCFL(dictCF["CFL"]);
    myProblem.setPrecision(dictCF["Precision"]);
    myProblem.setMaxNbOfTimeStep(dictCF["MaxNbOfTimeStep"]);
    myProblem.setTimeMax(dictCF["TimeMax"]);
    myProblem.setFreqSave(dictCF["FreqSave"]);
    myProblem.setFileName(dictCF["ResultFileName"]);
    myProblem.saveConservativeField(True);
    #myProblem.setVerbose(True);
    myProblem.usePrimitiveVarsInNewton(True);
    if(dictCF["spaceDim"]>1):
      myProblem.saveVelocity();
      pass

    myProblem.initialize()

    ok = myProblem.run()

    if (ok):
      print( "Simulation " + dictCF["ResultFileName"] + " is successful !" );
      pass
    else:
      print( "Simulation " + dictCF["ResultFileName"] + "  failed ! " );
      pass

    print( "------------ End of calculation !!! -----------" );

    myProblem.terminate()
    
    # TODO use a helper object here.
