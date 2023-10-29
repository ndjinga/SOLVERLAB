#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import os
import traceback
import fnmatch as FM
import pprint as PP

from PyQt5 import QtCore, QtGui, QtWidgets
import xml.etree.ElementTree as ET
import xyzpy.utilsXyz as UXYZ
import xyzpy.classFactoryXyz as CLFX
import xyzpy.helpsFactoryXyz as HLFX

from copy import deepcopy
import salomepy.iconsUser as IUSR
import xyzpy.loggingXyz as LOG

logger = LOG.getLogger()

verbose = False
debug = verbose
debugTooltip = False
verboseEvent = verbose

"""
note: future delegate stuff
http://stackoverflow.com/questions/7175333/howto-create-delegate-for-qtreewidget
http://www.riverbankcomputing.com/pipermail/pyqt/2011-August/030440.html
"""

"""the first model item of tree IS the model"""

####################################################
class TreeXmlXyzDelegate(QtWidgets.QItemDelegate):

  #notdone = True

  def paint(self, painter, option, index):
    """delegate paint for TreeXmlXyz"""

    """
    #for future
    if index.column() == -1: #Custom Draw Column 0
      icon = IUSR.getIconFromName("run")
      if icon != None:
        if self.notdone:
          for i in dir(option): print "option",option.rect.getRect(),i
          self.notdone=False
        ax, ay, aw, ah = option.rect.getRect()
        newRect = QtCore.QRect(ax, ay, ax, ah)
        #icon.paint(painter, option.rect) #You'll probably want to pass a different QRect
        icon.paint(painter, newRect) #You'll probably want to pass a different QRect
      return super(TreeXmlXyzDelegate, self).paint(painter, option, index)
    """

    """
    #debug
    if verbose: print 'TreeXmlXyzDelegate.paint',self,index.row(),index.column()
    if index.column() == 0 and index.row() == 0:
      print "option.rect.getRect()", option.rect.getRect()
    """

    modelItem = index.model() #PyQt5.QtCore.QAbstractItemModel
    treeWidget = modelItem.parent()
    treeXmlXyzItem = treeWidget.itemFromIndex(index)

    if index.column() == 0:
      #set _rect for contextmenu position ax, ay, aw, ah
      treeXmlXyzItem._rect = option.rect.getRect()
      return super(TreeXmlXyzDelegate, self).paint(painter, option, index)

    elif index.column() == 1: #Custom Draw Column 1
      if treeXmlXyzItem.isInBold:
        option.font.setWeight(QtGui.QFont.Bold)
        treeXmlXyzItem.setTextColor(1, treeWidget.colorInBold)
        treeXmlXyzItem.setTextColor(0, treeWidget.colorInBold)
      return super(TreeXmlXyzDelegate, self).paint(painter, option, index)

    else: #else Use the standard routine for other columns
      return super(TreeXmlXyzDelegate, self).paint(painter, option, index)

  def createEditor(self, parent, option, index):
    """delegate create editor for TreeXmlXyz, only on column 1"""

    #parent is PyQt5.QtWidgets.QWidget
    #option is PyQt5.QtWidgets.QStyleOptionViewItemV4
    #index is PyQt5.QtCore.QModelIndex

    if index.column() == 0: return None #no edition for column 0

    if verbose: logger.info('TreeXmlXyzDelegate.createEditor')
    if index.column() == 1: #free edition
      #return super(TreeXmlXyzDelegate, self).createEditor(parent, option, index) #default editor
      #search from factory for editor
      modelItem = index.model() #PyQt5.QtCore.QAbstractItemModel
      treeWidget = modelItem.parent()
      treeXmlXyzItem = treeWidget.itemFromIndex(index)

      try:
        xlmElement = treeWidget.mapXlmData[treeXmlXyzItem] #xlmElement of treeWidget.xmlData
        #if verbose: print "xlmElement.name", xlmElement.tag, xlmElement.text , str(xlmElement.attrib)
      except:
        logger.warning("no key in the local model (mapXlmData) for this item, no editor")
        return None
        #print "WARNING: no key in the local model (mapXlmData) for this item, set default editor"
        #anEditor = super(TreeXmlXyzDelegate, self).createEditor(parent, option, index) #default editor, no control
        #return anEditor

      aXyzClass = CLFX.getXyzClassFromName(xlmElement.attrib) #factory
      if verbose: logger.info("create editor for class '%s'" % xlmElement.attrib)
      if aXyzClass == None:
        logger.warning("no key in Xyz class factory for item '%s', no editor" % xlmElement.attrib)
        return None
        #anEditor = super(TreeXmlXyzDelegate, self).createEditor(parent, option, index) #default editor, no control
        #return anEditor

      #for i in dir(treeWidget): print "treeWidget",i
      #print "treeWidget.getController",treeWidget.getController()
      controller = treeWidget.getController()
      if controller != None:
        xlmNode = treeWidget.mapXlmData[treeXmlXyzItem]
        treePath = xlmNode.attrib["treePyName"]
        #print "treePath",treePath
        #WARNING anInstance have to be const in MVC model,
        #get anInstance in model only to get createEditor
        anInstance = controller._model.getValueByTreePyName(treePath)
        #print "anInstance",anInstance
      else:
        try:
          anInstance = aXyzClass(xlmElement.text)
        except: #if click on column 1 of type not immutable CrescendoCnf for example
          if verbose:
            logger.warning("can't set editor for '%s' with value '%s'" % (xlmElement.attrib, xlmElement.text))
          anInstance = aXyzClass() #default
          #return None #no!

      if verbose: logger.info("TreeXmlXyz create anInstance %s" % anInstance.__class__.__name__)

      if not hasattr(anInstance, "createEditor"):
        mess = "No direct edition of '%s'\ntry context menu EditDialog (mouse right click)"
        QtWidgets.QMessageBox.warning(parent, "warning", mess % anInstance.__class__.__name__)
        anEditor = None #default
        return None

      try:
        anEditor = anInstance.createEditor(parent)
      except:
        import traceback
        trace = traceback.format_exc()
        mess = "problem in method createEditor for '%s'" % anInstance.__class__.__name__
        logger.warning(mess + "\n" + trace)
        QtWidgets.QMessageBox.warning(parent, "warning", mess)
        anEditor = None #default

      if anEditor == None:
        #could be None, choosing or not...
        return None

      if type(anEditor) == str:
        QtWidgets.QMessageBox.warning(parent, "warning", anEditor)
        return None

      if verbose:
        logger.info("TreeXmlXyzDelegate.createEditor '%s' '%s' %s" % \
                    (anInstance, xlmElement.text.strip(), anEditor.__class__.__name__))
      #anInstance, if immutable is the value
      #treeXmlXyzItem.text is the xml text string value
      anEditor.setValue(anInstance)
      return anEditor

    if index.column() == 2: #ComboBox edition only column 1
      return None #no edition for column 2
      """
      #a QComboBox example... for debug
      combo = QtWidgets.QComboBox(parent)
      combo.addItems(["AnExemple 0", "AnExemple 1", "AnExemple etc..."])
      combo.setCurrentIndex(1)
      return combo
      """
    return None #no edition for other! columns.

  def setModelData(self, editor, model, index):
    #editor is PyQt5.QtWidgets.QLineEdit (for example)
    #model is PyQt5.QtCore.QAbstractItemModel
    #index is PyQt5.QtCore.QModelIndex
    if index.column() == 1:
      #TODO replace this stuff by XyzEditor.getValue()...
      if issubclass(editor.__class__, QtWidgets.QLineEdit):
        data = editor.text()
      elif issubclass(editor.__class__, QtWidgets.QComboBox):
        data = editor.currentText()
      elif issubclass(editor.__class__, QtWidgets.QSpinBox):
        data = editor.value()
      elif issubclass(editor.__class__, QtWidgets.QDoubleSpinBox):
        data = editor.value()
      else:
        logger.warning("setModelData: problem reaching result of editor %s" % editor)
    else: #no exec modifications
      return
    if verbose: logger.info('TreeXmlXyzDelegate.setModelData %s' % data)
    model.setData(index, QtCore.QVariant(data))


####################################################
class TreeXmlXyzItem(QtWidgets.QTreeWidgetItem):

  def __init__(self, parent=None, treeWidget=None):
    super(TreeXmlXyzItem, self).__init__(parent)
    self._className = self.__class__.__name__ #shortcut
    #self.setFlags( self.flags() | QtCore.Qt.ItemIsEditable | QtCore.Qt.ItemIsSelectable)
    self.setFlags( self.flags() | QtCore.Qt.ItemIsEditable | QtCore.Qt.ItemIsSelectable)
    self.setExpanded(True)
    self.isEmpty = True
    self.treeWidget = treeWidget
    self.isInBold = False
    #self.setDisabled(True) #andself.isDisabled()
    #self.setFirstColumnSpanned(True) #True to hide columns[1:]
    self._rect = None #view rect of item (set from paint) #ax, ay, aw, ah
    self._along = None
    self._filename = None # tooltip file name with expandvar ${xxx}
    self.tooltip_1 = None # tootip column 1 value throught xml

  def __hash__(self):
    """
    make it hashable for python 3 (when used as dict key)
    https://stackoverflow.com/questions/1608842/types-that-define-eq-are-unhashable
    """
    return id(self)

  def setFromXml(self, data):
    if data == None:
      logger.warning("problem for data as None")
      return
    if type(data) == str:
      logger.warning("TreeXmlXyzItem problem unexpected data type as str '%s'" % data)
      return
    try:
      aText = "_%s_" % data.attrib["index"]
    except:
      try:
        aText = self._toStr(data.tag)
      except:
        logger.warning("TreeXmlXyzItem unexpected data type %s" % data.__class__.__name__)
        aText = "TreeXmlXyzItem unexpected data type %s" % data.__class__.__name__
    self.setText(0, aText)

    boo = {'True': True, 'False': False}
    if 'hidden' in data.attrib:
      self.setHidden(boo[data.attrib["hidden"]])

    if 'icon' in data.attrib:
      icon = IUSR.getIconFromName(data.attrib["icon"])
      self.setIcon(0, icon)

    import debogpy.debug as DBG
    #DBG.write("data.attrib", data.attrib, True)

    # tooltips functionalities are tricky and have to be fixed... (2019)
    nameAttribute = data.tag
    parentItem = self.parent()
    if parentItem != None:
      parentTag = self.treeWidget.mapXlmData[self.parent()]
      try:
        nameClass = parentTag.attrib["typeClass"] #typeClass can be absent
        tooltip = HLFX.getCommonToolTip(nameClass, nameAttribute)
        if tooltip == None and "file__" in data.tag: #case set tooltip file__xx as basename
          tooltip = os.path.basename(data.text)
      except:
        tooltip = None

      self.tooltip_1 = None
      if "tooltip_1" in data.attrib:
        self.tooltip_1 = data.attrib["tooltip_1"] #tooltip_1 can be absent

      if tooltip != None:
        self.setToolTip(0, tooltip)
        if debugTooltip:
          logger.warning("Tooltip for %s.%s is '%s'" % (nameClass, nameAttribute, tooltip))
      else:
        # self.setToolTip(0, "No tooltip")
        if debugTooltip:
          logger.warning("No tooltip for %s.%s:\n%s" % (nameClass, nameAttribute, PP.pformat(data.attrib)))

    ashort, along = self.filterToStr(data.text)
    self.setData(1, QtCore.Qt.EditRole, along) #self.setText(1, along)

    try:
      display = UXYZ.toStrForTreeView( along, self.treeWidget.formats_treeview )
      #print "treeWidget.formats_treeview", self.treeWidget.formats_treeview, along, display
      self.setData(1, QtCore.Qt.DisplayRole, display)
    except:
      pass

    ashort3, along3 = self.filterToStr(data.attrib)
    super(TreeXmlXyzItem, self).setData(2, QtCore.Qt.EditRole, along3)

    self._along = along

    if "$" in along:  # may be environ variable... expand for tooltip (as filename)
      self._filename = os.path.expandvars(along)
    else:
      self._filename = along

    filename = self._filename

    ###############################################
    # files/directories for Iradina...
    if os.path.isfile(filename):
      # print("setData is file '%s'" % filename)
      # super(TreeXmlXyzItem, self).setData(0, QtCore.Qt.ToolTipRole, filename)
      # super(TreeXmlXyzItem, self).setData(1, QtCore.Qt.ToolTipRole, along)
      # super(TreeXmlXyzItem, self).setData(0, QtCore.Qt.DisplayRole, os.path.basename(filename))
      # super(TreeXmlXyzItem, self).setData(1, QtCore.Qt.EditRole, along)  # editrole override displayrole
      # super(TreeXmlXyzItem, self).setData(1, QtCore.Qt.DisplayRole, "")
      basename = os.path.basename(filename)
      super(TreeXmlXyzItem, self).setData(1, QtCore.Qt.ToolTipRole, "%s\n%s" % (basename, filename))

    if os.path.isdir(filename):
      # print("setData is dir '%s'" % filename)
      basename = os.path.basename(filename)
      super(TreeXmlXyzItem, self).setData(1, QtCore.Qt.ToolTipRole, "%s\n%s" % (basename, filename))

    return

  def setData(self, column, role, value):
    """
    http://qt-project.org/doc/qt-4.8/qt.html#ItemDataRole-enum
    column is integer
    role is integer
    value is PyQt5.QtCore.QVariant
    """
    #if column == 1 and role == QtCore.Qt.BackgroundColorRole:
    #  return qVariantFromValue(QtCore.QColor(QtCore.Qt.red))

    #print("setData type: %s value '%s'" % (type(value), str(value)))

    if column == 0: #file__xx to check as pseudo-browser current file in disk
      try:
        v = value.toString()
      except:
        try:
          v = u"%s" % value # str(value) plante sur 'température' python2
        except Exception as e:
          msg = "problem setData from value type '%s':\n'%s'" % (type(value), value.__repr__())
          logger.critical(msg)
          raise Exception(msg)

      if "file__" == v[0:6]:
        #print "column 0 file role %i value %s" % (role, v)
        vv = v.replace("file_", "")
        return super(TreeXmlXyzItem, self).setData(column, role, vv)
      else:
        return super(TreeXmlXyzItem, self).setData(column, role, value)


    if column > 2 : #unexpected ... yet (sait-on jamais)
      return super(TreeXmlXyzItem, self).setData(column, role, value)

    if role != QtCore.Qt.EditRole:
      return super(TreeXmlXyzItem, self).setData(column, role, value)

    #####################################
    #QtCore.Qt.EditRole
    #####################################
    # print("setDataNew EditRole %s %s" % (type(value), str(value)))
    ashort, along = self.filterToStr(value)

    if column == 2:
      super(TreeXmlXyzItem, self).setData(column, QtCore.Qt.EditRole, along)
      super(TreeXmlXyzItem, self).setData(column, QtCore.Qt.DisplayRole, ashort)
      super(TreeXmlXyzItem, self).setData(column, QtCore.Qt.ToolTipRole, along)
      return

    if self.isEmpty and column == 1: #first init: no check
      #if verbose: print "item is empty, first init no check '%s'" % value.value()
      self.isEmpty = False
      # print("setDataNew EditRole %s %s" % (type(value), str(value)))

      super(TreeXmlXyzItem, self).setData(column, QtCore.Qt.EditRole, along)
      #eviter: pour qLineEdit on se retrouve avec le ashort=DisplayRole dans le createEditor
      #super(TreeXmlXyzItem, self).setData(column, QtCore.Qt.DisplayRole, ashort)
      if self.tooltip_1 is None:
        super(TreeXmlXyzItem, self).setData(column, QtCore.Qt.ToolTipRole, along)
      else:
        super(TreeXmlXyzItem, self).setData(column, QtCore.Qt.ToolTipRole, self.tooltip_1)

      """
      #automatic basename in tooltip if file... could be time consuming
      if os.path.isfile(along):
        #print "isfile",along
        super(TreeXmlXyzItem, self).setData(0, QtCore.Qt.ToolTipRole, os.path.basename(along))
      return
      """

      return
      #return super(TreeXmlXyzItem, self).setData(column, role, value)

    if column == 1:
      if verbose:
        logger.info("TreeXmlXyzItem.setData editRole column %i %s" % (column, along))

      treeWidget = self.treeWidget
      controller = treeWidget.getController()

      res, anIntermediateLocalModel = self.testNewValue(value)
      if debug: logger.info("%s setData anIntermediateLocalModel:\n%s" % (res, anIntermediateLocalModel))
      if res == "ko":
        return
      if res == "ok":
        if controller != None:
          controller.replaceModel(anIntermediateLocalModel)
      if res == "force":
        self.isInBold = True
        super(TreeXmlXyzItem, self).setData(column, QtCore.Qt.EditRole, along)
        #eviter ce qui suit:
        #pour qLineEdit on se retrouve avec le ashort=DisplayRole dans le createEditor
        #super(TreeXmlXyzItem, self).setData(column, QtCore.Qt.DisplayRole, ashort)
        #super(TreeXmlXyzItem, self).setData(column, QtCore.Qt.ToolTipRole, along)

      return

  def testNewValue(self, value):
    #test new value (for column 1) through xlm data tree before accept
    treeWidget = self.treeWidget
    xlmNode = treeWidget.mapXlmData[self]
    xmlValue = self._toStr(xlmNode.text)
    itemValue = str(self.text(1)).strip()
    try:          #QString not in qt5
      newValue = str(value.toString()) #QVariant
    except:
      newValue = str(value)
    if debug: logger.info("TreeXmlXyzItem.testNewValue newValue '%s'" % newValue)

    #modification check
    #controller = treeWidget.getController()
    if True: ######## controller ==  None: #no MVC pattern
      ashort, along = self.filterToStr(newValue)
      xlmNode.text = UXYZ.toStrForXml(along) #set local model Xml, xml reading previous blanks format
      if verbose:
        logger.info("modifications is checked through model: itemValue: '%s', '%s'" % (itemValue, along))

      #try to create intermediate local model, catch errors
      try:
        if debug: logger.info("modelTemporaryBad ? %s" % treeWidget.modelTemporaryBad)
        if along == "!!":
          print(aBugForTest) # fait expres inexisting var ?
        anIntermediateLocalModel = UXYZ.fromXml(treeWidget.xmlData)
        #no errors... accept modifications in local model xml
        if debug: logger.info("modifications OK ! %s" % anIntermediateLocalModel.__class__.__name__)
        treeWidget.modelTemporaryBad = False
        return "ok", anIntermediateLocalModel #have to set initial model of controller

      except Exception as e:
        errorMessage = str(e) #.split("\n")
        treePath = xlmNode.attrib["treePyName"]
        mess = "\nSome errors occurs modifying data:\n   on %s\n" % treePath
        if treeWidget.modelTemporaryBad == False:
          #this mess is for first error encountered in UXYZ.fromXml(treeWidget.xmlData)
          #it could be not the good one if multiples bad
          mess += "   %s\n\n" % errorMessage
        else:
          mess += "   at least one error\n\n"
        mess += "Force change value to '%s' (at your own risks) ?... " % along
        if verbose:
          traceback.print_exc() #better explicit verbose problem
        logger.error(mess)
        rep = QtWidgets.QMessageBox.question(self.treeWidget, "WARNING", mess, QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
        if rep == QtWidgets.QMessageBox.No:
          logger.info("Answer is No")
          logger.info("restore old value '%s'" % itemValue)
          xlmNode.text = itemValue #restore xmlData
          #no need to set initial model of controller...
          return "ko", None
        else:
          #local modification, model is temporary? bad...
          logger.info("Answer is Yes")
          logger.info("""not restoring old value '%s'
model is temporary? incorrect value '%s'""" % (itemValue, along))
          #have to set initial model of controller later
          treeWidget.modelTemporaryBad = True
          return "force", None


  def removeFirstEmptyLines(self, aText):
    """replace first empty lines by nothing (for tooltip)"""
    lines = aText.split("\n")
    idep = 0
    for line in lines:
      if self.isEmptyString(line):
        idep += 1
        continue
      break
    if idep == 0: return aText
    res = ""
    for line in lines[idep:]:
      res += line + "\n"
    return res[:-1] #without last rc

  def removeFirstBlanks(self, aText):
    """replace first blancks by nothing (for tooltip)"""
    return aText.lstrip()

  def isEmptyString(self, aLine):
    """replace tab and rc and blancks by nothing to test rest"""
    ifEmpty = aLine.replace("\t", "").replace("\n", "").replace(" ", "")
    if ifEmpty == "":
      return True
    return False

  def filterToStr(self, theTextOrDict):
    """
    avoid to long strings, with to long lines to be in tree cells
    return tuple (shortString, longString) to use in
    (text cell, tooltip cell)
    """
    if theTextOrDict == None: return ("", "")

    if issubclass(theTextOrDict.__class__, int) or \
       issubclass(theTextOrDict.__class__, float):
      textOrDict = str(theTextOrDict)
      #print("filterToStr %s->'%s'" % (theTextOrDict.__class__, textOrDict))
    else:
      textOrDict = theTextOrDict

    if type(textOrDict) == str:
      atext = textOrDict
      if atext == "": return ("", "")
      if atext[-1] == "\n": atext = atext[:-1]
      if self.isEmptyString(atext): return ("", "")
      if atext==None: return ("", "")
      atext = self.removeFirstEmptyLines(atext)
      atext = self.removeFirstBlanks(atext)
      if len(atext) > 1500: #trunck it for tooltip also not too long...
        atext = atext[0:1500] + "\n...etc..."
      if atext.count("\n") < 1 : #accept 1 or 2 lines
        if len(atext) <= 30:
          return (atext, atext)
        else:
          return (atext[0:30] + "...",  atext)
      else:
        return (atext.split("\n")[0][0:30] + "...", atext)

    try: # possibly pyqt5 returns unicode, QString obsolete
      if type(textOrDict) == str:
        # print("text unicode: %s -> %s" % (type(textOrDict), str(textOrDict)))
        return self.filterToStr(str(textOrDict))
    except:
      pass

    if type(textOrDict) == dict: #attributes
      if textOrDict == {}: return ("", "")
      res = ""
      for k in sorted(textOrDict.keys()):
        res += str(k) + "='" + str(textOrDict[k]) + "'\n"
      return (self.filterToStr(res)[0], res[:-1])

    if issubclass(textOrDict.__class__, QtCore.QVariant):
      #return self.filterToStr(textOrDict.toString().toUtf8())
      return self.filterToStr(str(textOrDict.toString()))

    mess = "problem only text or unicode or dict or QVariant: %s -> %s" % (type(textOrDict), str(textOrDict))
    logger.warning( mess)

    res = ("?type text?: '%s'" % str(type(textOrDict)), str(textOrDict))
    return res

  def newfilterToStr(self, textOrDict):
    """
    avoid to long strings, whith to more lines to be in tree cells
    return tuple (shortString, longString) to use in (cell, tooltip)
    """
    if textOrDict == None: return ("", "")
    if type(textOrDict) == str:
      atext = textOrDict
      ifEmpty = atext.replace("\t", "").replace("\n", "").replace(" ", "")
      if ifEmpty == "": return ("", "")
      if atext[-1] == "\n": atext = atext[:-1]
      if atext == "None":
        return ("", "")
      if atext.count("\n") < 2 : #accept 1 or 2 lines
        if len(atext) <= 30:
          return (atext, atext)
        else:
          return (atext[0:30] + "...",  atext)
      else:
        return (atext.split("\n")[0] + "...", atext)
    if type(textOrDict) == dict: #attributes
      if textOrDict == {}:
        return ("", "")
      res = ""
      for k in sorted(textOrDict.keys()):
        res += str(k) + "='" + str(textOrDict[k]) + "'\n"
      return (self.filterToStr(res)[0], res[:-1])
    if issubclass(textOrDict.__class__, QtCore.QVariant):
      #return self.filterToStr(textOrDict.toString().toUtf8())
      return self.filterToStr(str(textOrDict.toString()))
    res = ("?type text?: " + str(type(textOrDict)), str(textOrDict))
    return res

  def _toStr(self, text):
    """avoid to long strings, whith to more lines to be in tree cells"""
    ashort, along = self.filterToStr(text)
    return ashort

_expanded = [True]

####################################################
class TreeXmlXyz(QtWidgets.QTreeWidget):
  """
  Is Widget Item of a View of pattern MVC

  It receive events/request via pyqt signal from his
  Controller with his method receiveRequestToView(self, aRequestFromController)
  The aRequestFromController is a string of XML syntax
  (have to parse and scrute to know that View have do do).

  It send events/request to his Model (via his Controller)
  with self.getController().sendRequestToController(aRequestToController)
  The aRequestToController is a string of XML syntax
  an example is in TreeXmlXyzItem.setData
  """

  refresh = QtCore.pyqtSignal()
  quickEditFiles = QtCore.pyqtSignal(list)
  index = [0] #unambigous objectName

  class COLS:
    labels = ['Tag', 'Text', 'Attributes']
    nbLabels = len(labels)+1
    Tag = 0
    Text = 1
    Attributes = 2

  def __init__(self, parent=None):
    super(TreeXmlXyz, self).__init__(parent)
    self._className = self.__class__.__name__ #shortcut
    self._verboseEvent = verboseEvent #quick existing
    #default QTreeWidgetItem class, could be modified in inherided
    self._TreeXmlXyzItemClass = TreeXmlXyzItem #default QTreeWidgetItem
    #objectName = "TreeXmlXyz"+str(self.index)
    objectName = str(self._className)+str(self.index)
    self.index[0] += 1 #unambigous objectName
    self.setObjectName(objectName)
    self.setWindowTitle(self.objectName())
    self.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
    self._lastReceiveRequest = None #firstly only for unittest
    #really could be None if no use in view without MVC pattern
    self._controller = None

    #self.header().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)
    #treeView->header()->setStretchLastSection(false);
    #treeView->header()->setSectionResizeMode(1, QHeaderView::Stretch);
    self.header().setStretchLastSection(True)
    #self.header().setSectionResizeMode(1, QtWidgets.QHeaderView.Stretch)
    self.setColumnWidth(0, 300)

    self.setHeaderLabels(self.COLS.labels)
    """
    #set in app.setFont
    self.setFont(QtGui.QFont("Monospace", 9))
    font = QtGui.QFont(self.font()) #copy: self.font() is const
    font.setFamily("Monospace")
    font.setPointSize(9)
    self.setFont(font)
    """
    self.setAlternatingRowColors(True)
    pal=self.palette()
    pal.setColor(pal.Base, QtGui.QColor(200,230,230))
    # alternatingRowColors do not work on S990, but ok on conda activate py3qt5 oops
    # pal.setColor(pal.AlternateBase, QtGui.QColor(100,130,130))
    # print("TreeXmlXyz AlternatingRowColors %s" % self.alternatingRowColors() )
    pal.setColor(pal.Text, QtGui.QColor(0,0,0))
    self.setPalette(pal)

    self.itemClicked.connect(self._handleItemClicked)
    self.itemChanged.connect(self._handleItemChanged)
    self.setItemDelegate(TreeXmlXyzDelegate())
    self.xmlData = None
    self.mapXlmData = None #mapXlmData[modelItem]=xlmElement of self.xmlData
    self.prefixShortcut = "Ctrl+"
    whoami = os.getenv("USERNAME")
    self._thirdColumnShown = True # not False as next flip/flop
    # logger.info("tree view third column shown %s as %s" % (not self._thirdColumnShown, whoami))
    self.showHideThirdColumnAction() # flip/flop
    self._createActions()
    self.contextMenus = self._createContextMenus()
    self.rootModelItem = None
    self.theHeader=self.header()
    self.theHeader.setSectionsClickable(True) #qt5

    self.theHeader.sectionClicked.connect(self.headerContextMenuEvent)
    self.currentDir = os.getenv("HOME") #to open xml file
    self.modelTemporaryBad = False
    self.colorInBold = QtGui.QColor("red")
    self.formats_treeview = UXYZ.FORMATS_TREEVIEW
    """
    self.theHeaderItem=self.headerItem()
    print "headerItem flags", self.theHeaderItem.flags()
    self.itemClicked.connect(self.headerContextMenuEvent)
    print "Item flags", self.itemUnderMouseAtContextMenuEvent.flags()
    self.theHeaderItem.setFlags(self.itemUnderMouseAtContextMenuEvent.flags())
    print "self.headerItem", self.theHeaderItem, self.itemUnderMouseAtContextMenuEvent
    """

  def setController(self, controller):
    """really could be None if no use in view without MVC pattern"""
    if self._controller == None:
      if verbose: logger.info("setController for %s" % self.objectName())
      self._controller = controller
      return
    raise Exception("TreeXmlXyz.setController done yet for %s as %s" % (self.objectName(), self.getController().objectName()))

  def getController(self):
    """to get (for example) sendRequest method of controller api"""
    return self._controller

  def lastReceiveRequest(self):
    """firstly in base class only for unittest"""
    return self._lastReceiveRequest

  def receiveRequestToView(self, strXmlRequest):
    if verboseEvent:
      logger.info("%s %s receiveRequestToView:\n%s" % (self.__class__.__name__, self.objectName(), strXmlRequest))
    self._lastReceiveRequest = strXmlRequest # firstly only for unittest
    aRequest = UXYZ.fromXml(strXmlRequest)
    if True: #TODOtry:
      typeRequest = aRequest.typeRequest
      if typeRequest == "newModel":
        xmlData = self.getController().getModel()
        if debug:
          logger.info("new model:\n %s" % UXYZ.prettyPrintET(xmlData))
        self.setFromXml(xmlData)
      if typeRequest == "clearModel":
        self.clearModel()
    else: #except:
      logger.warning("something wrong in request \n'%s'" % strXmlRequest)
    return True

  def getSelectedItems(self):
    """as selectedItems() for convenience"""
    return self.selectedItems()

  def getExpandedAndSelectedItems(self):
    itemCurrent = self.rootModelItem
    itemsExpanded=[]
    itemsSelected=[]
    self._getExpandedAndSelectedItems(itemCurrent, itemsExpanded, itemsSelected)
    return (itemsExpanded, itemsSelected)

  def _getExpandedAndSelectedItems(self, item, itemsExpanded,  itemsSelected):
    """Return the expanded items."""
    if item==None:
      return
    else:
      if item.isExpanded():
        itemsExpanded.append(item)
      if item.isSelected():
        itemsSelected.append(item)
      for i in range(item.childCount()):
        self._getExpandedAndSelectedItems(item.child(i), itemsExpanded,  itemsSelected)
      return

  def setExpandedRecursive(self, item, value):
    item.setExpanded(value)
    for i in range(item.childCount()):
      self.setExpandedRecursive(item.child(i), value)

  def expandSelected(self):
    items = self.selectedItems()
    for item in items:
      self.setExpandedRecursive(item, True)

  def collapseSelected(self):
    items = self.selectedItems()
    for item in items:
      self.setExpandedRecursive(item, False)

  def _setDictOfTreePath(self, items):
    """dict is faster than list to find keys"""
    res = {}
    for item in items:
      xlmNode = self.mapXlmData[item]
      treePath = xlmNode.attrib["treePyName"]
      if treePath in res:
        if treePath == '': continue #TODO duplicate same treePath key '' exists! why?
        logger.warning("duplicate same treePath key '%s':\npossible with pure Xml, impossible with Xml from structure Xyz" % treePath)
        #print res
        #traceback.print_stack()
        #print plante
      else:
        #print "treePath expanded",treePath
        res[treePath] = True
    return res

  def setExpandedAsPrevious(self, dictOfTreePathExpanded, item=None):
    """set expanded items if in TreePath keys corresponds"""
    if item==None:
      itemCurrent = self.rootModelItem
    else:
      itemCurrent = item
    xlmNode = self.mapXlmData[itemCurrent] #treeWidget.mapXlmData[itemCurrent]
    try: #pb if xlmNode == None or no attrib "treePyName"
      treePath = xlmNode.attrib["treePyName"]
      expanded = dictOfTreePathExpanded.get(treePath, False) #default=False
    except:
      expanded = False  #case xml without treePyName
    itemCurrent.setExpanded(expanded) #expand or collapse current item
    #recursive walk in tree
    for i in range(itemCurrent.childCount()):
      self.setExpandedAsPrevious(dictOfTreePathExpanded, itemCurrent.child(i))
    return

  def userExpandModel(self, aList, item=None, Verbose=False):
    """
    set expanded items if in TreePath keys corresponds
    filename patterns '*,??',
    warning not for '[' ']' as 'alist[*]'
    no found way to quote meta-character '[' and ']'
    https://docs.python.org/2/library/fnmatch.html
    """
    if item==None:
      itemCurrent = self.rootModelItem
    else:
      itemCurrent = item
    xlmNode = self.mapXlmData[itemCurrent] #treeWidget.mapXlmData[itemCurrent]
    treePath = xlmNode.attrib["treePyName"].split(".")[-1] #basename treePyName
    expanded = False
    for k in [''] + aList: #root have "" as treePyName
      if FM.fnmatch(treePath, k):
        #print "treePath expanded ok ", treePath, k
        expanded = True
        break
    if Verbose: logger.info("treePath expanded %s for '%s'" % (expanded, treePath))
    itemCurrent.setExpanded(expanded) #expand or collapse current item
    #recursive walk in tree
    for i in range(itemCurrent.childCount()):
      self.userExpandModel(aList, itemCurrent.child(i), Verbose)
    return

  def ExpandItems(self, aList, item=None, Verbose=False):
    """aList of exact TreePyName(s), not patterns"""
    if item==None:
      itemCurrent = self.rootModelItem
    else:
      itemCurrent = item
    xlmNode = self.mapXlmData[itemCurrent] #treeWidget.mapXlmData[itemCurrent]
    # treePath = xlmNode.attrib["treePyName"].split(".")[-1] #basename treePyName
    treePath = xlmNode.attrib["treePyName"]  # treePyName
    if Verbose: logger.info("ExpandItems '%s' %s" % (treePath, aList))

    expanded = None  # as do not change
    if treePath in aList:
      # logger.info("treePath expanded ok %s" % treePath)
      expanded = True  # to expand

    if expanded is not None:
      if Verbose: logger.info("treePath expanded %s for '%s'" % (expanded, treePath))
      itemCurrent.setExpanded(expanded) # expand current item

    #recursive walk in tree
    for i in range(itemCurrent.childCount()):
      self.ExpandItems(aList, itemCurrent.child(i), Verbose)
    return

  def CollapseItems(self, aList, item=None, Verbose=False):
    """aList of exact TreePyName(s), not patterns"""
    if item==None:
      itemCurrent = self.rootModelItem
    else:
      itemCurrent = item
    xlmNode = self.mapXlmData[itemCurrent] #treeWidget.mapXlmData[itemCurrent]
    # treePath = xlmNode.attrib["treePyName"].split(".")[-1] #basename treePyName
    treePath = xlmNode.attrib["treePyName"]  # treePyName
    if Verbose: logger.info("CollapseItems %s %s" % (treePath, aList))

    expanded = None  # as do not change
    if treePath in aList:
      # logger.info("treePath collapsed ok %s" % treePath")
      expanded = False  # to collapse

    if expanded is not None:
      if Verbose: logger.info("treePath expanded %s for '%s'" % (expanded, treePath))
      itemCurrent.setExpanded(expanded) #expand or collapse current item

    #recursive walk in tree
    for i in range(itemCurrent.childCount()):
      self.CollapseItems(aList, itemCurrent.child(i), Verbose)
    return

  def setSelectedAsPrevious(self, dictOfTreePathSelected, item=None):
    """set expanded items if in TreePath keys corresponds"""
    if item==None:
      itemCurrent = self.rootModelItem
    else:
      itemCurrent = item
    xlmNode = self.mapXlmData[itemCurrent] #treeWidget.mapXlmData[itemCurrent]
    try: #pb if xlmNode == None or no attrib "treePyName"
      treePath = xlmNode.attrib["treePyName"]
      selected = dictOfTreePathSelected.get(treePath, False) #default=False
    except:
      selected = False  #case xml without treePyName
    if selected:
      selModel=self.selectionModel()
      index = self.indexFromItem(itemCurrent)
      selModel.select(index, selModel.Select | selModel.Rows)
    #recursive walk in tree
    for i in range(itemCurrent.childCount()):
      self.setSelectedAsPrevious(dictOfTreePathSelected, itemCurrent.child(i))
    return

  def refreshModel(self):
    if self._controller != None:
      controller = self._controller
      aDataXml =  controller.getModel() #get copy as xml
      if aDataXml != None: #no model
        self.setFromXml(aDataXml)
        if verboseEvent: logger.info("refreshModel %s" % self.objectName())
      else:
        if verboseEvent: logger.info("refreshModel empty model %s" % self.objectName())
    else:
      logger.info("refreshModel with no controller")

  def editDialog_withoutController(self):
    """obsolete this treeXml.editDialog is only for case without controller"""
    import xyzpy.guiXyz.dialogXmlXyz as DXYZ

    items = self.selectedItems()
    if len(items) == 0:
      QtWidgets.QMessageBox.warning(self.treeWidget, "warning", "select one item please")
      return
    if len(items) > 1:
      QtWidgets.QMessageBox.warning(self.treeWidget, "warning", "select only one item please")
      return
    item = items[0]
    xlmNode = self.mapXlmData[item]

    if self._controller != None:
      mess = "treeXmlXyz.editDialog is only for case without controller"
      QtWidgets.QMessageBox.warning(self, "warning", mess)
      return
      #controller = self._controller
      #aDataXml =  controller.getModel() #get copy as xml
      #self.setFromXml(aDataXml)
    else:
      logger.info("treeXmlXyz.editDialog without controller")
    self.widDialog = DXYZ.DialogXmlXyz()
    dataXyz = UXYZ.fromXml(xlmNode)
    try:
      dataStrXml = dataXyz.toStrXml()
    except:
      mess = "No context menu EditDialog for '%s'\ntry direct edition (mouse double left click)"
      QtWidgets.QMessageBox.warning(self.treeWidget, "warning", mess % dataXyz.__class__.__name__)
      return

    if debug: print("dialogXmlXyz.setFromXml dataStrXml:\n",dataStrXml)

    self.widDialog.setFromXml(dataStrXml) #(self.xmlData)
    self.widDialog.resize(500, 500)
    self.widDialog.show()

  def expandToDepth0(self):
    self.expandToDepth(0)

  def clearModel(self):
    """only raz local model, not controller model"""
    if verboseEvent: logger.info("TreeXmlXyz clearModel")
    self.xmlData = None
    self.mapXlmData = {}
    self.clear()
    self.rootModelItem = None
    #print "clearModel done"
    #TODO for i in dir(self): print "treeXml",  i

  def setFromBase(self, aBase):
    """
    set data aBase by value with conversion to Elementtree
    """
    data = ET.fromstring(aBase.toStrXml())
    self.setFromXml(data)
    return

  def setFromXml(self, data):
    """initial populate the tree with QTreeWidgetItem items, raz model"""
    if verboseEvent: logger.info("TreeXmlXyz setFromXml")
    verb = self._verboseEvent
    self._verboseEvent = False #avoid many print
    currentRoot = self
    currentColumnWidths = [self.columnWidth(i) for i in [0,1]]
    itemsExpanded, itemsSelected = self.getExpandedAndSelectedItems()
    #if verboseEvent:
    #  print "setFromXml itemsSelected", len(itemsSelected)
    #  print "setFromXml itemsExpanded", len(itemsExpanded)
    dictOfTreePathSelected = self._setDictOfTreePath(itemsSelected)
    dictOfTreePathExpanded = self._setDictOfTreePath(itemsExpanded)
    self.clearModel() #TODO verify it is good...
    self.xmlData = data
    self._setFromXml(self.xmlData, currentRoot)
    self.setExpandedAsPrevious(dictOfTreePathExpanded)
    self.setSelectedAsPrevious(dictOfTreePathSelected)
    self._verboseEvent = verb
    for i in [0,1] : self.setColumnWidth(i, currentColumnWidths[i])

  def _setFromXml(self, data, item):
    """
    recursive populate the tree with QTreeWidgetItem items
    with TreeXmlXyzItem by default, or else, set in self._TreeXmlXyzItemClass
    """
    rowItem = self._TreeXmlXyzItemClass(parent=item, treeWidget=self)
    if self.mapXlmData == {}:  #first insert in new tree
      self.rootModelItem = rowItem

    self.mapXlmData[rowItem] = data
    if data == None:
      logger.warning("problem for data as None")
      return
    rowItem.setFromXml(data)
    for row in data:
      self._setFromXml(row, rowItem)
    return

  def _handleItemClicked(self, item, column):
    if self._verboseEvent:
      logger.info("TreeXmlXyz._handleItemClicked(Item, column) %s %s" % (item.text(column), column))
    return

  def _handleItemChanged(self, item, column):
    if self._verboseEvent:
      logger.info("TreeXmlXyz._handleItemChanged(Item, column) %s %s" % (item.text(column), column))
    return

  def xxx_event(self, event):
    if self._verboseEvent:
      logger.info("TreeXmlXyz.event %s %s" % (self.strEvent(event), event.type()))
    try:
      self.itemUnderMouseAtContextMenuEvent = self.itemAt(event.pos())
      #print "self.event with event pos '%s' '%s'"% (self.itemUnderMouseAtContextMenuEvent.text(0), self.theHeaderItem.text(0))
      return super(TreeXmlXyz, self).event(event)
    except:
      #print "****self.event no event pos"
      return super(TreeXmlXyz, self).event(event)

  def xxx_mousePressEvent(self, event):
    logger.info("mousePressEvent %s %s" % (event.x(),event.y()))
    return super(TreeXmlXyz, self).mousePressEvent(event)

  def strEvent(self, event):
    """catch readable explicit type of event"""
    mapEvent={}
    anInt = event.type()
    #for i in dir(event): print "event", i, getattr(event, i).__class__.__name__, getattr(event, i)
    for i in dir(event):
      anAttr = getattr(event, i)
      if anAttr.__class__.__name__ == "Type":
        #print "event", i, anAttr.__class__.__name__, anAttr
        mapEvent[int(anAttr)] = str(i)
    try:
      res = mapEvent[int(event.type())]
    except:
      res = "unknown event '%i'" % anInt
    #print "strEvent ****** ", event.type(), res
    return res

  def headerContextMenuEvent(self, arg):
    #print "TreeXmlXyz.headerContextMenuEvent at column %i" % arg
    pos = QtGui.QCursor.pos()
    theMenu = self.contextMenus[-1]
    #theMenu.exec_(self.mapToGlobal(QtCore.QPoint(10,10)))
    theMenu.exec_(pos)
    return

  def contextMenuEvent(self, event):
    self.columnAtContextMenuEvent= self.columnAt(event.pos().x())
    item = self.itemAt(event.pos())

    #print "x,x", event.pos().x(), item._rect[0] #ax, ay, aw, ah
    #painter.boundingRect(QRectF(0,0,0,0),Qt.AlignLeft,self.displayText())
    #ax, ay, aw, ah = option.rect.getRect()

    self.itemUnderMouseAtContextMenuEvent = self.itemAt(event.pos())
    self.itemsSelectedAtContextMenuEvent  = self.getSelectedItems()
    self.posAtContextMenuEvent = event.pos()
    self.posXAtContextMenuEvent = event.pos().x()
    self.theHeaderItem=self.headerItem()
    if verboseEvent:
      logger.info("TreeXmlXyz.contextMenuEvent at column %i position %i %i" % \
                   (self.columnAtContextMenuEvent, event.pos().x(), event.pos().y()))

    theMenu = self._getContextMenuModelItem()
    if theMenu != None:
      theMenu.exec_(self.mapToGlobal(event.pos()))
      return
    else:
      #logger.warning("can't launch context menu at column %i" % self.columnAtContextMenuEvent)
      QtWidgets.QTreeView.contextMenuEvent(self, event) #standard super menu event
      return

  def _getTreePathsOfListItems(self, aList):
    #treePathsSelected = [self.mapXlmData[item].attrib["treePyName"] for item in self.itemsSelectedAtContextMenuEvent]
    res = []
    for item in aList:
      try:
        res.append(self.mapXlmData[item].attrib["treePyName"])
      except:
        logger.warning("problem get treePyName on '%s'" % (str(item)))
        pass
    return res

  def _getContextMenuModelItem(self, ):
    theMenu = None
    controller = self.getController()
    treePathsSelected = self._getTreePathsOfListItems(self.itemsSelectedAtContextMenuEvent)
    if verbose: logger.info("treePathsSelected %s" % treePathsSelected)
    if controller == None:
      logger.warning("no controller: no menu from controller")
      return self.contextMenus[self.columnAtContextMenuEvent]

    #print "**xevent, xitem", self.posAtContextMenuEvent.x(), self.itemUnderMouseAtContextMenuEvent._rect[0]

    #menu expand etc... only on beginning of rows

    #if self.columnAtContextMenuEvent == 1: #obsolete
    if self.itemUnderMouseAtContextMenuEvent == None: #if happends sometimes
      return None
    if self.itemUnderMouseAtContextMenuEvent._rect == None: #if happends sometimes
      logger.warning("itemUnderMouseAtContextMenuEvent._rect is None: no menu")
      return None
    if self.posAtContextMenuEvent.x() > self.itemUnderMouseAtContextMenuEvent._rect[0]+10:
      #+10 to avoid collapse/expand in clicking "+" item menu expand
      #here search only if controller and model with actions
      nbSel = len(treePathsSelected)
      if nbSel == 1:
        treePath = treePathsSelected[0]
        if verboseEvent: logger.info("_getMenuModelItem single selection for menu of '%s'" % treePath)
        return controller.getDefaultContextMenuForItem(treePath)
      elif nbSel == 0:
        if verboseEvent: logger.info("_getMenuModelItem zero selection")
        return controller.getDefaultContextMenuZeroSelection()
      else:
        if verboseEvent: logger.info("_getMenuModelItem multiple selection")
        return controller.getDefaultContextMenuMultipleSelection(treePathsSelected)

    if theMenu == None:
      try:
        return self.contextMenus[self.columnAtContextMenuEvent]
      except:
        logger.warning("can't get context menu at column %i" % self.columnAtContextMenuEvent)
        return None
    else:
      return theMenu

  def _createActions(self):
    self.actionEdit_withoutController = self._createAction('EditDialog', None, 'Edit dialog widget', self.editDialog_withoutController, 'dialogwidget')
    self.expandActions = [
      #self.obsolete_actionEdit,
      self._createAction('Expand all', None, 'expand all items', self.expandAll, 'expand'),
      self._createAction('Expand first level', None, 'expand all items', self.expandToDepth0, None),
      self._createAction('Expand selected', None, 'expand items selected', self.expandSelected, None),
      self._createAction('Collapse all', None, 'collapse all items', self.collapseAll, 'collapse'),
      self._createAction('Collapse selected', None, 'collapse items selected', self.collapseSelected, None),
      ]
    self.actionShowHideThirdColumn = self._createAction( 'ShowHideThirdColumn', '3', 'Show/Hide third column', self.showHideThirdColumnAction, 'addColumn')
    self.actionRefresh = self._createAction('Refresh', '0', 'Refresh View', self.refreshModel, 'refresh')
    self.otherActions = [
      #self._createAction( 'LoadFileXmlAction', '1', 'Load File Xml', self.loadFileXmlAction, 'openxml'),  #shortcut exists but no explicit, use TreeXmlXyzMainWidget to get menu
      self.actionShowHideThirdColumn,
      #self.actionRefresh, #shortcut exists but no explicit, use TreeXmlXyzMainWidget to get menu
      ]
    for a in self.otherActions:
      self.addAction(a)

  def _createAction(self, Name, Shortcut, ToolTip, Call, Icon=None, Enable=True):
    # print("%s._createAction '%s'" % (self.__class__.__name__, Name))
    action = QtWidgets.QAction(Name, self)
    action.setIconVisibleInMenu(True)
    if Shortcut!=None: action.setShortcut(self.prefixShortcut+Shortcut)
    action.setToolTip(ToolTip)
    if Icon!=None:
      action.setIcon(IUSR.getIconFromName(Icon))
    action.setEnabled(Enable)
    action.triggered.connect(Call)
    return action

  def _createContextMenus(self):
    """each column have a different menu"""
    menus = {}

    #outside column event
    menuGeneral = QtWidgets.QMenu("General", self)
    for a in self.otherActions:
      menuGeneral.addAction(a)
    menus[-1] = menuGeneral

    #column 0 event
    menuExpand = QtWidgets.QMenu("Expand", self)
    for action in self.expandActions:
      menuExpand.addAction(action)
    menus[0] = menuExpand

    #column 1 event: the default menu action column 1: see _getMenuModelItem
    menuItemAction = QtWidgets.QMenu("ItemActions", self)
    #menuItemAction.addAction(self.obsolete_actionEdit)
    menus[1] = menuItemAction

    #column 2 event
    #do not work: runtimeError: super-class __init__() of type QMenu was never called
    #menus[2] = deepcopy(menuGeneral)
    menuAttributes = QtWidgets.QMenu("Attributes", self)
    menuAttributes.addAction(self.actionShowHideThirdColumn)
    menus[2] = menuAttributes
    return menus

  def showHideThirdColumnAction(self):
    """flip flop hide/show"""
    if self._thirdColumnShown == True: #hide third
      self.header().setStretchLastSection(False)
      #self.header().setSectionResizeMode(1, QtWidgets.QHeaderView.Stretch)
      self.hideColumn(2)
      #self.resizeColumnToContents(1)
      #self.setColumnWidth(1, self.columnWidth(1)+10) #a marge...
      self.setColumnWidth(1, 110)
      self._thirdColumnShown = False

    else: #show third
      self.header().setStretchLastSection(True)
      self.header().setSectionResizeMode(2, QtWidgets.QHeaderView.Stretch)
      self.showColumn(2)
      #self.resizeColumnToContents(1)
      #self.setColumnWidth(1, self.columnWidth(1)+10) #a marge...
      self.setColumnWidth(1, 110)
      self._thirdColumnShown = True

  def loadFileXmlAction(self):
    #self.setFromXml(tag)
    if verboseEvent: logger.info("LoadFileXml")
    #if self._model != None: #TODO question save if existing
    nameFile = QtWidgets.QFileDialog.getOpenFileName(self, 'Load file xml', self.currentDir, "(*.xml *.XML)" )[0]
    nameFile = str(nameFile)
    if nameFile == "": return True #cancel
    realPath = os.path.realpath(nameFile)
    if not os.path.isfile(realPath):
      QtWidgets.QMessageBox.warning(self._desktop, "warning",
          "Load xml file: not a file \n'%s'" % realPath)
      return False
    try:
      #aXml = UXYZ.getEtFromFileXml(realPath)
      aXml = ET.parse(realPath).getroot()
      logger.warning("reading file xml done..., set widget now...")
      self.setFromXml(aXml)
      #self.setModel(aData)
      #self.RefreshCrescendoModelSignal.emit(None)
      return True
    except:
      trace = traceback.format_exc()
      QtWidgets.QMessageBox.warning(self, "error",
          "Load xml file: problem loading from \n'%s'\n\n%s" % (realPath, trace))
      return False

####################################################
class TreeXmlXyzMainWidget(QtWidgets.QWidget):

  refresh = QtCore.pyqtSignal()
  quickEditFiles = QtCore.pyqtSignal(list)

  def __init__(self, parent=None):
    super(TreeXmlXyzMainWidget, self).__init__(parent)

    self.widget = TreeXmlXyz()
    self.__createActions()
    self.__addToolBars()
    hbox = self.toolBars[0]
    # Set the layout
    vbox = QtWidgets.QVBoxLayout()
    vbox.addWidget(hbox)
    vbox.addWidget(self.widget)

    self.setLayout(vbox)
    self.resize(600,700)
    self.setWindowTitle("TreeXmlXyzMainWidget")

  def __addToolBars(self):
    self.toolBars = []
    tb = QtWidgets.QToolBar("EditXml")
    for action in self.actions:
      tb.addAction(action)
    self.toolBars.append(tb)

  def __createActions(self):
    """create local widget actions"""
    self.actions = []
    self.actions.append(self.__createAction(
      "ReadXmlFileTreeXmlXyz", None, "Read file xml", self.widget.loadFileXmlAction, "openxml" ) )
    self.actions.append(self.__createAction(
      "ShowHideThirdColumn", None, "Show/Hide Attributes Xml", self.widget.showHideThirdColumnAction, "addColumn" ) )
    self.actions.append(self.__createAction(
      "RefreshTreeXmlXyz", None, "Refresh xml", self.widget.refreshModel, "refresh" ) )
    self.actions.append(self.__createAction(
      "ClearTreeXmlXyz", None, "Clear xml", self.widget.clearModel, "clearModel" ) )

  def __createAction(self, Name, Shortcut, ToolTip, Call, Icon=None):
    action = QtWidgets.QAction(Name, self)
    if Shortcut!=None: action.setShortcut(self.prefixShortcut+Shortcut)
    action.setToolTip(ToolTip)
    if Icon!=None:
      action.setIcon(IUSR.getIconFromName(Icon))
    action.triggered.connect(Call)
    return action


####################################################
def example_createTreeXmlXyz():
  treeWid = TreeXmlXyz()
  treeWid.setWindowTitle("Elementary Xml Viewer")
  treeWid.resize(500, 500)
  treeWid.show()
  treeWid.loadFileXmlAction()
  return treeWid

def example_createTreeXmlXyzMainWidget():
  treeWid = TreeXmlXyzMainWidget()
  treeWid.setWindowTitle("Xml Viewer")
  treeWid.resize(500, 500)
  treeWid.show()
  #treeWid.widget.loadFileXmlAction()
  return treeWid

if __name__ == '__main__':
  # python ...xyzpy/quiXyz/treeXmlXyz.py
  from salomepy.onceQApplication import OnceQApplication
  app = OnceQApplication()
  #treeWid = example_createTreeXmlXyz()
  treeWid = example_createTreeXmlXyzMainWidget()
  app.exec_()
