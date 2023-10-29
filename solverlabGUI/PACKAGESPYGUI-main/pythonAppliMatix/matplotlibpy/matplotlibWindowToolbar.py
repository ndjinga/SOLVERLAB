#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


"""
from embedding_in_qt4_wtoolbar.py in matplotlib http://matplotlib.org
"""

import os
import sys

from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *

from widgetpy.comboBoxDialog import ComboBoxDialog
from salomepy.strEvent import *
from time import sleep
import xyzpy.loggingXyz as LOG
import pprint as PP

import numpy as np
from matplotlib.figure import Figure
#from matplotlib.backend_bases import key_press_handler
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT # as NavigationToolbar

import matplotlib
#print("matplotlib", matplotlib.__file__)

#from matplotlib.backends import qt_compat
use_pyside = False

import matplotlib.pyplot as plt

import xyzpy.loggingXyz as LOG

########################interactive on
plt.ion()

logger = LOG.getLogger()
progname = os.path.basename(__file__)
verbose = False


##############################################################################
class NavigationToolbarCustom(NavigationToolbar2QT):

  def sizeHint(self):
    # For some reason, self.setMinimumHeight doesn't seem to carry over to
    # the actual sizeHint, so override it instead in order to make the
    # aesthetic adjustments noted above.
    size = super(NavigationToolbar2QT, self).sizeHint()
    if verbose:
      print("NavigationToolbarCustom sizeHint as ", size)
    return size


##############################################################################
class MatplotlibWindowToolbar(QMainWindow):
  def __init__(self, parent=None):
    super(MatplotlibWindowToolbar, self).__init__(parent)
    if verbose: print("!!! MatplotlibWindowToolbar NEW !!!")
    self._create_main_frame()
    self.set_curve()
    self.on_draw()
    self.docks = [] #no docks yet
    self.toolBars = [] #no toolBars yet
    #self._timer = Nkey_press_handlerone
    self._deltaTime = 1000
    self._controller = None
    self.setWindowTitle("MatplotlibWindowToolbar")

  def _create_main_frame(self, figure=None):
    main_frame = QWidget()
    
    if figure == None:
      fig = Figure((8.0, 5.0), dpi=100)
    else:
      fig = figure
    main_frame.canvas = FigureCanvas(fig)
    main_frame.canvas.setParent(main_frame)
    main_frame.canvas.setFocusPolicy(Qt.StrongFocus)
    main_frame.canvas.setFocus()
    
    main_frame.mpl_toolbar = NavigationToolbarCustom(main_frame.canvas, main_frame)
    main_frame.mpl_toolbar.setMinimumHeight(20)
    # main_frame.mpl_toolbar.setMaximumHeight(15)
    main_frame.mpl_toolbar.setIconSize(QSize(20, 20))

    main_frame.canvas.mpl_connect('key_press_event', self.on_key_press)
    #main_frame.canvas.mpl_connect('pick_event', self.on_pick)

    vbox = QVBoxLayout()
    vbox.addWidget(main_frame.mpl_toolbar)
    vbox.addWidget(main_frame.canvas)  # the matplotlib canvas
    main_frame.setLayout(vbox)
    self.setCentralWidget(main_frame)
    main_frame.figure = fig
    self.main_frame = main_frame
  
  def setFromPlotPandas(self, aPlot):
    """
    aPlot is from convenience pandas DataFrame.plot(),
    aPlot is in fact a matplotlib Axes
    """
    figure = aPlot.get_figure()
    widget = figure.canvas.parentWidget()
    self.setFigure(figure)
    #because widget automatic show with pandas
    #then hide because figure is in self now
    widget.hide()
  
  def setFigure(self, aFigure):
    self._create_main_frame(aFigure)
    
  def getFigure(self):
    return self.main_frame.figure
    
  def setController(self, controller):
    """goal is get trace of plotted figures in a browser through controller"""
    self._controller = controller
    
  def setNewFigureToController(self):
    if self._controller != None:
      if verbose: print("setNewFigure")
      try:
        self._controller.addNewFigure(self.getFigure())
      except:
        logger.error("problem for controller.addNewFigure")
    
  def _getDefaultGridSize(self, nbc):
    """get a default grid from a total number of plot"""
    if nbc == 1: return (1, 1)
    if nbc == 2: return (1, 2)
    if nbc == 3: return (2, 2)
    if nbc == 4: return (2, 2)
    if nbc == 5: return (2, 3)
    if nbc == 6: return (2, 3)
    if nbc >= 7 and nbc <= 9: return (3, 3)
    if nbc > 9: return (3, nbc/3 + 1)
    
  def _getDefaultColor(self):
    """get a default sequential colors of plot in series"""
    return ["blue", "red", "green", "cyan", "magenta", "yellow", "black"]
    
  def set_curve(self, x=None, y=None, 
                suptitle=None, title=None, 
                title_x=None, title_y=None, 
                color=None, linewidth=.5, comment=None):
    """simple for only one curve"""
    self._create_main_frame() #create new figure
    if x is None and y is None: return
    if color == None: 
      col = self._getDefaultColor()[0]
    else:
      col = color
    #not so evident that is better, invalidate configure subplots
    #self.getFigure().set_tight_layout(True)
    ax = self.getFigure().add_subplot(1, 1, 1)
    if x is not None: ax.plot(x, y, color=col, linewidth=linewidth, label='curve_0')
    if suptitle is not None: self.getFigure().suptitle(suptitle, fontsize=14)
    if title is not None: ax.set_title(title)
    if title_x is not None: ax.set_xlabel(title_x)
    if title_y is not None: ax.set_ylabel(title_y)
    if comment is not None:
      # https://matplotlib.org/3.1.1/gallery/text_labels_and_annotations/annotation_demo.html?highlight=text
      ax.annotate(comment, xy=(5, 5), xycoords='figure pixels', fontsize=8)
    self.setNewFigureToController()

  def set_2_curves(self, xy=None, 
                   suptitle="", title=None, titles_xy=None, 
                   color=None, linewidth=[1., .5], comment=None):
    """for 2 curves on same graph with 2 independent y_scales"""
    self._create_main_frame() #create new figure
    if xy == None: return
    if color == None: 
      col = self._getDefaultColor()
    else:
      col = color
    
    ax1 = self.getFigure().add_subplot(1, 1, 1)
    ax1.plot(xy[0], xy[1], color=col[0], linewidth=linewidth[0], linestyle="-", label='curve_0')
    ax1.set_xlabel(titles_xy[0])
    ax1.set_ylabel(titles_xy[1], color=col[0])
    #ax1.set_label('left')

    ax2 = ax1.twinx()
    ax2.plot(xy[0], xy[2], color=col[1], linewidth=linewidth[1], linestyle="-", label='curve_1')
    ax2.set_ylabel(titles_xy[2], color=col[1])
    #ax1.set_label('right')

    self.getFigure().suptitle(suptitle) #, fontsize=14)
    if title is not  None: ax1.set_title(title)
    if comment is not None:
      # https://matplotlib.org/3.1.1/gallery/text_labels_and_annotations/annotation_demo.html?highlight=text
      ax1.annotate(comment, xy=(5, 5), xycoords='figure pixels', fontsize=8)
    self.setNewFigureToController()

  def nbins_scott(self, x):
    """
    scott, or else for histogram
    https://docs.scipy.org/doc/numpy/reference/generated/numpy.histogram.html#numpy.histogram
    only version numpy v1.11, here 2017 numpy v1.9
    """
    xmax = max(x)
    xmin = min(x)
    nb = len(x)
    sigma = np.std(x)
    if sigma == 0:
      logger.warning("histogram nbins_scott with sigma 0, default nbins=10")
      return 10
    nbins = (xmax-xmin)/(3.5*sigma)*(nb**(1./3.))
    if np.isnan(nbins):
      print("!!!! nbins_scott calculus is NaN %s" % sigma)
      nbins = 10 #default
    return int(nbins)


  def set_n_histograms(self, xy=None, 
                   suptitle=None, title=None, titles_xy=None, 
                   color=None, linewidth=None, title_y=None):
    """
    for n histograms on same graph with 1 y_scale
    automatic binning with Scott
    """
    self._create_main_frame() #create new figure
    if xy == None: return
    if color == None: 
      col = self._getDefaultColor()
    else:
      col = color

    ax = self.getFigure().add_subplot(1, 1, 1)
    nbc = len(xy)

    for i in range(nbc):
      y = xy[i]
      """
      https://docs.scipy.org/doc/numpy/reference/generated/numpy.histogram.html#numpy.histogram
      only version numpy v1.11
      hist, bins = np.histogram(y, bins='auto')
      """
      nbins = max(10, self.nbins_scott(y)) #Scott      
      ax.hist(y, bins=nbins, color=col[i], alpha=0.5, label=titles_xy[i])
      #ax.hist(y, color=col[i], alpha=0.5, label=titles_xy[i])
      
    ax.grid(True)
    ax.legend(loc='upper right') #, fontsize=8)
    ax.set_ylabel("nb")
    #ax.set_xlabel(titles_xy[0])
    
    if suptitle != None: self.getFigure().suptitle(suptitle) #, fontsize=14)
    if title != None: ax.set_title(title)
    if title_y != None: ax.set_ylabel(title_y)
    self.setNewFigureToController()

  
  def set_n_curves(self, xy=None, 
                   suptitle=None, title=None, titles_xy=None, 
                   color=None, linewidth=None, title_y=None):
    """for n curves on same graph with 1 y_scale"""
    self._create_main_frame() #create new figure
    markers = ['x', 'o']
    if xy == None: return
    if color == None: 
      col = self._getDefaultColor()
    else:
      col = color

    ax = self.getFigure().add_subplot(1, 1, 1)
    nbc = len(xy)-1
    if linewidth == None: #regressive width for hidden superposed identicals curves
      lw = [.5+float(i)/2. for i in reversed(list(range(nbc)))]
    else:
      if issubclass(linewidth.__class__, list):
        lw = linewidth
      else:
        lw = [linewidth]*nbc #as constant

    x = xy[0]
    for i in range(nbc):
      y = xy[i+1]
      ax.plot(x, y, color=col[i], linewidth=lw[i], label=titles_xy[i+1])

    ax.legend(loc='upper left') #, fontsize=8)
    ax.set_xlabel(titles_xy[0])
    if suptitle != None: self.getFigure().suptitle(suptitle) #, fontsize=14)
    if title != None: ax.set_title(title)
    if title_y != None: ax.set_ylabel(title_y)
    self.setNewFigureToController()
  
  def set_n_plots_from_files( self, thefiles ):
    #a menu for what file in thefiles
    if type(thefiles) == list:
      if len(thefiles) > 1:
        self.combo = ComboBoxDialog("Select your file", thefiles, self)
        self.combo.exec_()
        afile = self.combo.choice
      else:
        if len(thefiles) == 1:
          afile = thefiles[0]
        else:
          afile = None
    else:
      afile = str(thefiles) #could be str for one file, not list
    
    if afile == None:
      print("set_n_plots_from_files: No file choice, cancel plot")
      return
    xy, titles_xy, atitle = self.readSimpleTableFile( afile )
    if xy != None:
      self.set_n_plots(xy=xy, suptitle=atitle, title=None, titles_xy=titles_xy, 
                       linewidth=None, title_y=None, grid_size=None)
    return

  def readSimpleTableFile(self, afile):
    """
    read a simple file with float values in columns, lines
    with ou without columns name header beginning with '#'
    separator blanks
    """
    import salomepy.readSimpleTableFile as RSTF
    return RSTF.readSimpleTableFile(afile)
  
  def set_n_plots(self, xy=None, 
                  suptitle=None, title=None, titles_xy=None, 
                  linewidth=None, title_y=None, grid_size=None):
    """for n distinct plots on same figure with 1 first x_scale"""
    
    """
    for Matplotlib 1.1.x only
    plot.subplot2grid((5, 2), (2, 0), rowspan=1, colspan=1)
    The first parameter of subplot2grid is the size of the grid. 
    The second parameter is the coordinate of one element of the grid.
    That element will be filled with a plot.
    Finally, we can specify the span of the plot,
    so that its spans over several elements of the grid.
    This gives a great deal of control over the layout, 
    without too much headaches. 
    There is a small gotcha however: the coordinates are in row, col order.
    Here's an example of subplot2grid in action.
    http://blog.marmakoide.org/?p=94
    """
    self._create_main_frame() #create new figure

    if xy == None: return
    nbc = len(xy)-1
    if grid_size == None:
      grid_sizec = self._getDefaultGridSize(nbc)
    else:
      grid_sizec = grid_size
    #print "grid_size", grid_sizec
    
    if linewidth == None:
      lw = [1 for i in range(nbc)]
    else:
      if issubclass(linewidth.__class__, list):
        lw = linewidth
      else:
        lw = [linewidth]*nbc #as constant

    x = xy[0]
    for i in range(nbc):
      y = xy[i+1]
      ax = self.getFigure().add_subplot(grid_sizec[0], grid_sizec[1] , i+1)
      ax.plot(x, y, linewidth=lw[i], label=titles_xy[i+1])
      ax.set_title(titles_xy[i+1])
      
    if suptitle != None: self.getFigure().suptitle(suptitle) #, fontsize=14)
    self.setNewFigureToController()
    return
    
  def set_appended_plots(self, xy=None, 
                         suptitle=None, title=None, titles_xy=None, color=None, 
                         linewidth=None, title_y=None, 
                         grid_size=None, position=None):
    """for n contiguous plots on same figure with 1 y_scale and n x_scale, like parex"""
    
    #https://scipy-lectures.github.io/intro/matplotlib/matplotlib.html#axes

    self._create_main_frame() #create new figure

    if xy == None: return
    nbc = len(xy)
    if grid_size == None:
      grid_sizec = (1, nbc)
    else:
      grid_sizec = grid_size
    #print "grid_size", grid_sizec
    
    if color == None: 
      col = self._getDefaultColor()
    else:
      col = color
    
    if linewidth == None: #regressive width for hidden superposed identicals curves
      lw = [.5+float(i)/2. for i in reversed(list(range(nbc)))]
    else:
      if issubclass(linewidth.__class__, list):
        lw = linewidth
      else:
        lw = [linewidth]*nbc #as constant
        
    if position == None:
      pos = [.1, .1, .85, .8] #[x0, y0, width ,height]
    else:
      pos = [float(i) for i in position]
    
    x0 = pos[0]
    y0 = pos[1]
    width = pos[2]/nbc
    height = pos[3]
    positions = []
    for i in range(nbc):
      xi = x0 + i * width
      positions.append([xi, y0, width, height])
    
    x = xy[0]
    axs = []
    for i in range(nbc):
      x = xy[i][0]
      y = xy[i][1]
      ax = self.getFigure().add_subplot(grid_sizec[0], grid_sizec[1] , i+1)
      axs.append(ax)
      ax.set_position(positions[i])
      ax.plot(x, y, color=col[i], linewidth=lw[i], label=titles_xy[i+1])
      ax.set_title(titles_xy[i+1], color=col[i])
      if i > 0: 
        ax.get_yaxis().set_visible(False)
      else:
        try:
          ax.yaxis.set_tick_params(length=0, width=1) #, which='major')
        except:
          #do not work in some version matplotlib
          logger.warning("problem for matplotlib version set_tick_params")
          pass
    
    if suptitle != None: self.getFigure().suptitle(suptitle) #, fontsize=14)
    if title_y != None: axs[0].set_ylabel(title_y)
    
    #or ...may be... shared_axis_demo.py
    yMinMax = axs[0].get_ylim()
    for ax in axs[1:]:
      yr = ax.get_ylim()
      if yMinMax[0] > yr[0]: yMinMax[0] = yr[0]
      if yMinMax[1] < yr[1]: yMinMax[1] = yr[1]
    for ax in axs: ax.set_ylim(yMinMax)
    
    for ax in axs[1:]:
      for xlabel_i in ax.get_xticklabels():
        #xlabel_i.set_fontsize(0.0)
        xlabel_i.set_visible(False) #it is reference pointer, not ax.get_xticklabels()[0]
        break

    """
    #do not work on pan/zoom
    self.canvas.draw() #have to be draw to set labels
    for ax in axs[1:]:
      labels = [item.get_text() for item in ax.get_xticklabels()]
      labels[0] = ''
      ax.set_xticklabels(labels)
    """
    self.setNewFigureToController()
        
  def on_draw(self):
    try:
      self.main_frame.canvas.draw()
    except:
      raise Exception("Problem MatplotlibWindowToolbar.on_draw (may be C/C++ object deleted)")

  def event(self, event):
    if verbose: logger.info("MatplotlibWindowToolbar.event: there is event: %s" % strEvent(event))
    return super(MatplotlibWindowToolbar, self).event(event)
  
  def contextMenuEvent(self, event):
    #theMenu = self.contextMenus[self.columnAtContextMenuEvent]
    #theMenu.exec_(self.mapToGlobal(event.pos()))
    logger.info("MatplotlibWindowToolbar.contextMenuEvent %s" % strEvent(event))
    return
      
  def on_key_press(self, event):
    logger.info("you pressed '%s' but key_press_handler is inactive" % event.key)
    #self.set_curves(self.shift) #for example...
    #self.on_draw()
    # implement the default mpl key press events described at
    # http://matplotlib.org/users/navigation_toolbar.html#navigation-keyboard-shortcuts
    #key_press_handler(event, self.canvas, self.mpl_toolbar)"""
  
  def on_pick(self, event):
    """
    The event received here is of the type
    matplotlib.backend_bases.PickEvent
    It carries lots of information, of which we're using
    only a small amount here.
    """
    logger.info("you picked %s but key_press_handler is inactive" % event.key)
    #box_points = event.artist.get_bbox().get_points()
    #msg = "You've clicked on a bar with coords:\n %s" % box_points
    #QMessageBox.information(self, "Click!", msg)
  
  def allTests(self):
    self._timer = QTimer()
    self._timer.timeout.connect(self._test)
    self._noTest = 0
    self._timer.start(self._deltaTime)
    return
  
  def _test(self, noTest=None):
    if verbose: logger.info("MatplotlibWindowToolbar._test %s %s" % (self._noTest, noTest))
    x = np.arange(0.0, 2.0, 0.01)
    x2 = np.arange(10.0, 30.0, 0.1)
    if len(x) != len(x2): 
      logger.warning("MatplotlibWindowToolbar x an x2 needs same size...") #but it is not obligatoire
    y1= np.sin(2*np.pi*x)
    y2= .5*np.sin(3*np.pi*x)
    y3= .2*np.sin(6*np.pi*x)
    
    suptitle = "a_suptitle"
    title = "a_title"
    xlabel = "x_label"
    ylabel = "y_label"
    y1label = "y1_label"
    y2label = "y2_label"
    y3label = "y3_label"
    
    xy12 = [x, y1, y2]
    xy12labels = [xlabel, y1label, y2label ]
    xy13 = [x, y1, y2, y3]
    xy13labels = [xlabel, y1label, y2label , y3label ]
    xy13_appended = [[x, y1], [x2, y2], [x, y3]]
    
    if noTest == None:
      noTestc = self._noTest
      self._noTest += 1
    else:
      noTestc = noTest
    
    if noTestc == 0:
      self.set_curve(x, y1, suptitle, title, xlabel, ylabel)
      self.on_draw()
      return
    if noTestc == 1:
      self.set_2_curves(xy12, suptitle, title, xy12labels)
      self.on_draw()
      return
    if noTestc == 2:
      self.set_n_curves(xy13, suptitle, title, xy13labels, title_y="y_label_unity")
      self.on_draw()
      return
    if noTestc == 3:
      self.set_n_plots(xy13, suptitle, title, xy13labels, title_y="y_label_unity")
      self.on_draw()
      return
    if noTestc == 4:
      self.set_appended_plots(xy13_appended, suptitle, title, xy13labels, title_y="y_label_unity")
      self.on_draw()
      return
      
    self._timer = None #stop timer garbage collecting
    if verbose: logger.info("MatplotlibWindowToolbar: no test: %s" % str(noTest))
    return
