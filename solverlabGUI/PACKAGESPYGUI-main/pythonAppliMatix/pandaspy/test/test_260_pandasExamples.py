#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


"""
have to verify matplotlib backend...
python #Python 2.7.10

import matplotlib
import matplotlib.pyplot as plt
plt.get_backend()
plt.switch_backend('QT5Agg') #default on my system
plt.get_backend()


import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
plt.ion()
import numpy as np
x = np.arange(0, 5, 0.1);
y = np.sin(x)
aPlot = plt.plot(x, y) #[<matplotlib.lines.Line2D object at 0x41c00b0>]
plt.show()
aPlot.get_figure().canvas.parentWidget()



import pprint as PP
import matplotlib
matplotlib.use('Qt5Agg')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.get_backend()
plt.ion()

periods = 100
ts = pd.Series(np.random.randn(periods), index=pd.date_range('1/1/2000', periods=periods))
ts = ts.cumsum()
df = pd.DataFrame(np.random.randn(periods, 4), index=ts.index, columns=list('ABCD'))
df = df.cumsum()
aPlot = df.plot() #<matplotlib.axes._subplots.AxesSubplot object at 0x45dd240>
aPlot.legend(loc='upper left')
#print "aPlot", PP.pformat(dir(aPlot.get_figure().canvas))
aPlot.get_figure().canvas.parentWidget()

aMgr = plt.get_current_fig_manager()
#print "aMgr", PP.pformat(dir(aMgr))
aMgr.window.close()
plt.close("all")
"""

import fnmatch
import os
import unittest
from collections import OrderedDict
import pprint as PP

import numpy as np


from PyQt5 import QtCore, QtGui
from salomepy.onceQApplication import OnceQApplication
import salomepy.utilsWorkdir as UTW
runningDir = UTW.getWorkdirDefault("TESTS")

try:
  import matplotlib
  # matplotlib.use('Qt5Agg')
  # This call to matplotlib.use() has no effect
  # because the backend has already been chosen;
  # matplotlib.use() must be called *before* pylab, matplotlib.pyplot,
  # or matplotlib.backends is imported for the first time.
  import matplotlib.pyplot as plt
  plt.switch_backend('QT5Agg') #default on my system
  from matplotlibpy.matplotlibWindowToolbar import MatplotlibWindowToolbar
  matplotlibOk = True
except:
  print("WARNING: no matplotlib library")
  matplotlibOk = False


try:
  # The pandas.np module is deprecated and will be removed from pandas in a future version. Import numpy directly instead
  # from pandas import *   # avoid that where np is pandas.np
  import pandas as pd
  import pandaspy.utilsPandasOscar as UPO
  from pandaspy.pandasOscarMainWidget import PandasTabWidget, PandasMainWidget
  pandasOk = True
except:
  print("ERROR: no pandas library, no tests")
  pandasOk = False

#ImportError: 'numexpr' not found. Cannot use engine='numexpr' for query/eval if 'numexpr' is not installed
#http://stackoverflow.com/questions/33909782/numexpr-not-found-in-pandas
#http://matplotlib.org/users/installing.html
#python -m pip install -U pip setuptools
#python -m pip install matplotlib
#pd.show_versions()

"""
# skip for python 3
try:
  import numexpr #it fails using pandas matix
  engine = 'numexpr'
except:
  engine = 'python'
"""
engine = 'python'


if pandasOk:
  verbose = True
  setTimer = True
  aDebug = False
  timers = []
  deltaTime = 1000
  withShow = True

user = os.getenv("USERNAME")

#aFiles = "./tmp/ResulatsOSCAR_IHMSortie/*/*.csv"
#aFiles = "./oscar_original/Livraison_DM2S/ResulatsOSCAR_IHMSortie/Terme_Source/Terme_Source_PF"


class TestCase(unittest.TestCase):

  if pandasOk:
    d1 = {'one' : [1., 2., 3., 4.],
          'two' : [4., 3., 2., 1.],
          'three' : [3.0, 3.1, 3.2, None],}
    df1 = pd.DataFrame(d1)
    df2 = pd.DataFrame(d1, index=['a', 'b', 'c', 'd'])
  
  def recursive_glob(self, rootdir='.', pattern='*'):
      return [os.path.join(looproot, filename)
              for looproot, _, filenames in os.walk(rootdir)
              for filename in filenames
              if fnmatch.fnmatch(filename, pattern)]
              
  def getHeadFile(self, aNameFile, nbLines=15):
    """
    get head of file, firsts nbLines lines, trunk long lines, only for debug
    """
    lines = ""
    with open(aNameFile) as myFile:
      for x in range(nbLines):
        try:
          line = next(myFile)
        except: #eof
          return lines
        if len(line) <= 70: 
          lines += line
        else:
          lines += line[0:70]+" ...\n"
    return lines

  def launchTimerPlt(self):
    #only one timer for all plt show
    #timer survives in [timers] and continue close periodically problem -> see test_999
    if setTimer:
      if len(timers) < 1:
        app = OnceQApplication()
        timer = QtCore.QTimer();
        timer.timeout.connect(plt.close)
        timer.start(deltaTime)
        timers.append(timer)
      
  def launchTimerFen(self, wid):
    #only one timer for wid plt show
    #timer survives in [timers] and continue close periodically problem -> see test_999
    if setTimer:
        app = OnceQApplication()
        timer = QtCore.QTimer();
        timer.timeout.connect(wid.close)
        timer.start(deltaTime)
        timers.append(timer)
      
  def print_exec(self, anExpr, anEval):
    if aDebug:
      print('\n%s:' % anExpr)
      print(anEval)
    return
  
  def test_005(self):
    self.assertEqual(plt.get_backend(), "Qt5Agg") #if not there is problem

  def test_010(self):
    self.assertTrue(pandasOk)
    
    class TestEssaiDF(pd.DataFrame):
      def __init__(self, *args, **kwargs):
        super(TestEssaiDF, self).__init__(*args, **kwargs)
        self._info = None

    df1 = self.df1
    df2 = self.df2
    dfe1 = TestEssaiDF(self.d1)
    dfe2 = TestEssaiDF(self.d1, index=['a', 'b', 'c', 'd'])
    if aDebug:
      print("\n***test_010 verbose:")
      expr = 'df1' ; self.print_exec(expr, eval(expr)) #, engine=engine))
      expr = 'df2' ; self.print_exec(expr, eval(expr)) #, engine=engine))
      expr = 'dfe1' ; self.print_exec(expr, eval(expr)) #, engine=engine))
      expr = 'dfe2' ; self.print_exec(expr, eval(expr)) #, engine=engine))
      
  def test_015(self):
    if not pandasOk: return
    info, df = UPO.getExampleDataFrameElementary()
    if not aDebug: return
    
    print('df:\n', df)
    res = df[df["X"] > 2]
    print('df[df["X"] > 2]:\n', res)
    print('df:\n', df)
    dfx = df
    nbrows = len(dfx)
    nbcols = len(dfx.columns)
    print("***columns %i row %i" % (nbcols, nbrows))
    for r in range(nbrows):
      print("row %i:\n" % r, dfx.iloc[r, :])
    for c in range(nbcols):
      print("col %i:\n"% c, dfx.iloc[:, c])
    dfx = res
    nbrows = len(dfx)
    nbcols = len(dfx.columns)
    print("***columns %i row %i" % (nbcols, nbrows))
    for r in range(nbrows):
      print("row %i:\n" % r, dfx.iloc[r, :])
    for c in range(nbcols):
      print("col %i:\n"% c, dfx.iloc[:, c])
      
  def test_020(self):
    #http://manishamde.github.io/blog/2013/03/07/pandas-and-python-top-10
    
    if not pandasOk: return
    d3 = {'int_col' : [1, 2, 6, 8, -1], 
          'float_col' : [0.1, 0.2, 0.3, 10.1, None],
          'other_col' : [3, 3.1, 3.2, 3.3, 3.4],}
    #avoid cotes in eval("d3"), cause problem 'str_col' : ['a','b',None,'c','a'],}
    df = pd.DataFrame(d3)
    if True:
      #expr = '' ; self.print_exec(expr, eval(expr, engine=engine))
      expr = "df" ; self.print_exec(expr, eval(expr)) #, engine=engine))

      expr = "df['int_col'][0]" ; self.print_exec(expr, eval(expr)) #, engine=engine))
      expr = "df['float_col'][0]" ; self.print_exec(expr, eval(expr)) #, engine=engine))
      
      res = df.loc[:,['float_col','int_col']]
      self.print_exec("df.loc[:,['float_col','int_col']]", res)

      res = df[['float_col','int_col']]
      self.print_exec("df[['float_col','int_col']]", res)

      #res = df[[1, 2]] #not pandas Salome8
      #self.print_exec("df[[1, 2]]", res)

      res = df.loc[[2,3],['float_col','int_col']]
      self.print_exec("df.loc[[2,3],['float_col','int_col']]", res)

      res = df.iloc[[2,3],[1, 2]]
      self.print_exec("df.iloc[[2,3],[1, 2]]", res)
      
      #Conditional indexing
      res = df[df['float_col'] > 0.15]
      self.print_exec("df[df['float_col'] > 0.15]", res)
      
      res = df[df['float_col'] == 0.1]
      self.print_exec("f[df['float_col'] == 0.1]", res)
      
      res = df[(df['float_col'] > 0.1) & (df['int_col']>2)]
      self.print_exec("df[(df['float_col'] > 0.1) & (df['int_col']>2)]", res)
      
      res = df[~(df['float_col'] > 0.1)]
      self.print_exec("df[~(df['float_col'] > 0.1)]", res)
      
      #Renaming columns. It copies the data to another DataFrame.
      res = df.rename(columns={'int_col' : 'some_other_name'})
      self.print_exec("df.rename(columns={'int_col' : 'some_other_name'})", res)
      self.print_exec("df", df)
      
      res.rename(columns={'some_other_name' : 'int_col'}, inplace = True)
      self.print_exec("df2.rename(columns={'some_other_name' : 'int_col'}, inplace = True)", res)
      
      #Handling missing values
      #The dropna can used to drop rows or columns with missing data (NaN). By default, it drops all rows with any missing entry.
      res = df.dropna()
      self.print_exec("df.dropna()", res)
      
      res = df.copy()
      mean = res['float_col'].mean()
      res = res['float_col'].fillna(mean)
      self.print_exec("res['float_col'].fillna( res['float_col'].mean() )", res)
      
      if aDebug: print("""
  Map, Apply
  Forget writing for loops while using pandas. 
  One can do beautiful vectorized computation by applying function over rows and columns
  using the map, apply and applymap methods.
    -The map operation operates over each element of a Series.
    -The apply is a pretty flexible function which, 
      as the name suggests, applies a function along any axis of the DataFrame. 
      The examples show the application of the sum function over columns. 
      (Thanks to Mindey in the comments below to use np.sum instead of np.sqrt in the example)
    -The applymap operation can be used to apply the function to each element of the DataFrame.
      """)
      res = df['float_col'].dropna().map(lambda x : 100 + x)
      self.print_exec("df['float_col'].dropna().map(lambda x : 100 + x)", res)

      self.print_exec("df", df)
      res = df.loc[:,['int_col','float_col']].apply(np.sqrt)
      self.print_exec("df.loc[:,['int_col','float_col']].apply(np.sqrt)", res)
      res = df.loc[:,['int_col','float_col']].apply(np.sum)
      self.print_exec("df.loc[:,['int_col','float_col']].apply(np.sum)", res)
      res = df[['int_col','float_col']].apply(np.sum)
      self.print_exec("df[['int_col','float_col']].apply(np.sum)", res)
      
      def some_fn(x):
        #print type(x)
        if type(x) is str:
          return 'applymap_' + x
        elif type(x) is np.int64:
          return -x
        else:
          return 100+x
      
      res = df.applymap(some_fn)
      self.print_exec("df.applymap(some_fn)", res)
      
      #Vectorized mathematical and string operations
      df = pd.DataFrame(data={"A":[1,2,3], "B":[1.1,1.2,1.3], "Z":["x","y","z"] })
      self.print_exec("df", df)
      
      df["C"] = df["A"]+df["B"]
      self.print_exec('df["C"] = df["A"]+df["B"]', df)
      
      df["Cbis"] = df.A-df.B*2
      self.print_exec('df["Cbis"] = df.A-df.B*2', df)
      
      df["D"] = df["A"]*3
      self.print_exec('df["D"] = df["A"]*3', df)
      
      df["E"] = np.sqrt(df["A"])
      self.print_exec('df["E"] = np.sqrt(df["A"])', df)
      
      df["F"] = df.Z.str.upper()
      self.print_exec('df["F"] = df.Z.str.upper()', df)
      
      df["G"] = df.A**2
      self.print_exec('df["G"] = df.A**2', df)
      
      if aDebug: print("""
  The groupby method let’s you perform SQL-like grouping operations.
  The example below shows a grouping operation performed with str_col columns entries as keys. 
  It is used to calculate the mean of the float_col for each key. 
  For more details, please refer to the split-apply-combine description on the pandas website.
  http://pandas.pydata.org/pandas-docs/dev/groupby.html
  By “group by” we are referring to a process involving one or more of the following steps
    -Splitting the data into groups based on some criteria
    -Applying a function to each group independently
    -Combining the results into a data structure
      """)
      grouped = df[['A']].groupby(df['Z']).mean()
      self.print_exec("df['A'].groupby(df['Z']).mean()", grouped)
      grouped = df[['A','B']].groupby(df['Z']).mean()
      self.print_exec("df['A','B'].groupby(df['Z']).mean()", grouped)
      
      if aDebug: print("""
      multiple columns as a function of a single column
      n.b. 'zip(*df...' as *args **kwargs
      """)
      def two_operations(x):
        return x*2, x**2
      
      df = pd.DataFrame({"A":[1.,2.,3.]})
      self.print_exec('df', df)
      #print df['A'].map(two_operations)
      
      #n.b. "zip(*df" as *args **kwargs
      df['*2'], df['**2'] = list(zip(*df['A'].map(two_operations)))
      self.print_exec("df['*2'], df['**2'] = zip(*df['A'].map(two_operations))", df)
      
      if aDebug: print("""
      single column as a function of multiple columns
      """)
      def sum_two_cols(series):
        return series['A'] + series['*2']
        
      df['*3'] = df.apply(sum_two_cols, axis=1)
      self.print_exec("df['*3'] = df.apply(sum_two_cols, axis=1)", df)
      
      if aDebug: print("""
      multiple columns as a function of multiple columns
      """)
      def squares(series):
        #return series['A'] + series['*2']
        return pd.Series({'pow1' : series['A'], 'pow2' : series['A']**2, 'pow4' : series['**2']**2, 'pow6' : series['**2']**3})
        
      res1 = df.apply(squares, axis = 1)
      self.print_exec("df.apply(squares, axis = 1)", res1)
      
      res = res1.describe()
      self.print_exec("quick stats: res.describe()", res)

      res = res1.cov()
      self.print_exec("covariance between suitable columns: res.cov()", res)

      res = res1.corr()
      self.print_exec("correlation between suitable columns: res.corr()", res)
      
      if aDebug: print("""
  Pandas supports database-like joins which makes it easy to link data frame
  has full-featured, high performance in-memory join operations idiomatically 
  very similar to relational databases like SQL. 
  These methods perform significantly better 
  (in some cases well over an order of magnitude better) 
  than other open source implementations (like base::merge.data.frame in R). 
  The reason for this is careful algorithmic design and internal layout of the data in DataFrame
    -left LEFT OUTER JOIN Use keys from left frame only
    -right RIGHT OUTER JOIN Use keys from right frame only
    -outer FULL OUTER JOIN Use union of keys from both frames
    -inner INNER JOIN Use intersection of keys from both frames
      """)
      
      other = pd.DataFrame({'A' : [1, 3, 4], 'some_val' : [91, 93, 94]})
      self.print_exec('df', df)
      self.print_exec('other', other)
      
      res = pd.merge(df,other,on='A',how='inner')
      self.print_exec("merge(df,other,on='A',how='inner')", res)
      res = pd.merge(df,other,on='A',how='outer')
      self.print_exec("merge(df,other,on='A',how='outer')", res)
      
      res = pd.merge(df,other,on='A',how='left')
      self.print_exec("merge(df,other,on='A',how='left')", res)
      res = pd.merge(df,other,on='A',how='right')
      self.print_exec("merge(df,other,on='A',how='right')", res)
      
  def test_025(self):
    if not pandasOk: return
    app = OnceQApplication()
    #http://pandas.pydata.org/pandas-docs/stable/visualization.html
    periods = 100
    ts = pd.Series(np.random.randn(periods), index=pd.date_range('1/1/2000', periods=periods))
    ts = ts.cumsum()
    df = pd.DataFrame(np.random.randn(periods, 4), index=ts.index, columns=list('ABCD'))
    df = df.cumsum()
    aPlot = df.plot()
    aPlot.legend(loc='upper left')
    #print "aPlot", PP.pformat(dir(aPlot.get_figure().canvas))
    self.launchTimerPlt() #aPlot.get_figure().canvas.parentWidget()
    app.exec_()

  def test_026(self):
    if not pandasOk: return
    app = OnceQApplication()
    #http://pandas.pydata.org/pandas-docs/version/0.15.0/visualization.html
    periods = 1000
    df = pd.DataFrame({"teta": np.linspace(0, np.pi*2*5, periods)} )
    #vectorized
    r = np.linspace(1, 0.1, periods)
    df["X"] = np.cos(df.teta)*r
    df["Y"] = np.sin(df.teta)*r
    aPlot = df.plot(x="X", y="Y", label="ellipse", legend=False)
    lim = 1.5
    aPlot.set_xlim(-lim,lim)
    aPlot.set_ylim(-lim,lim)
    self.launchTimerPlt()
    app.exec_()

    '''
    # patch AttributeError: 'FigureCanvasQTAgg' object has no attribute '_is_drawing'
    # in /volatile/common/miniconda3/envs/pyIradinaGUI3.6/lib/python3.6/site-packages/matplotlib/backends/backend_qt5.py
        def draw_idle(self):
        """Queue redraw of the Agg buffer and request Qt paintEvent.
        """
        # The Agg draw needs to be handled by the same thread matplotlib
        # modifies the scene graph from. Post Agg draw request to the
        # current event loop in order to ensure thread affinity and to
        # accumulate multiple draw requests from event handling.
        # TODO: queued signal connection might be safer than singleShot
        try:
          self._draw_pending = self._draw_pending
        except:
          print('inexisting _draw_pending')
          self._draw_pending = False
          self._is_drawing = False
        if not (self._draw_pending or self._is_drawing):
            self._draw_pending = True
            QtCore.QTimer.singleShot(0, self._draw_idle)
    '''

  def test_027(self):
    if not pandasOk: return
    app = OnceQApplication()
    #http://pandas.pydata.org/pandas-docs/version/0.15.0/visualization.html
    periods = 100
    df = pd.DataFrame({"teta": np.linspace(0, np.pi*2)} )
    df["X"] = np.sin(df.teta)
    df["Y"] = np.cos(df.teta)
    lim = 1.5
    aPlot = df.plot(x="X", y="Y", title="circle", legend=False, xlim=[-lim,lim], ylim=[-lim,lim])
    fen = MatplotlibWindowToolbar()
    """
    #something like that, see setFromPlotPandas()
    figure = aPlot.get_figure()
    widget = figure.canvas.parentWidget()
    fen = MatplotlibWindowToolbar()
    fen.setFigure(figure)
    widget.hide()
    """
    fen.setFromPlotPandas(aPlot)
    self.launchTimerPlt()
    app.exec_()
    
  def test_105(self):
    if not pandasOk: return
    app = OnceQApplication()
    fen = PandasMainWidget()
    fen.show()
    if not aDebug:
      self.launchTimerFen(fen)
    else:
      print("WARNING: test_105 no launchTimer as debug")
    app.exec_()

  def xxxxtest_120(self):
    if not pandasOk: return
    if aDirOscar==None:
      print("\nWARNING: test_120 oscar skipped")
      return
    #aFiles = aDirOscar+"/*/*/*.csv"
    #fileNames = glob.glob(aFiles)
    fileNames = self.recursive_glob(aDirOscar, "*.csv")
    #fileNames = ["/export/home/wambeke/oscar_original/Livraison_DM2S/ResulatsOSCAR_IHMSortie/IHMS/BilansCircuitPrimaire/Actinides_Primaire_Isotopes_Synthese.csv"]
    if verbose: print("nb files .csv found:", len(fileNames))
    for fileName in fileNames:
      if verbose: print("test read",fileName)
      try:
        infos, contents = UPO.read_oscar_csv(fileName)
      except:
        print("problem reading",fileName)
        print("firsts lines:\n%s" % self.getHeadFile(fileName))
        if verbose: #explicit problem
          infos, contents = UPO.read_oscar_csv(fileName)
      if False: #verbose: #wait for user eye
        print("infos:\n%s" % infos)
        x = input("\n***************************************\nsuite?...")
        
  def xxxxtest_130(self):
    if not pandasOk: return
    if aDirOscar==None:
      print("\nWARNING: test_130 oscar skipped")
      return
    app = OnceQApplication()
    fen = PandasTabWidget()

    #aFiles = aDirOscar+"/*/*/*.csv"
    fileNames = self.recursive_glob(aDirOscar, "*.csv")
    #print fileNames
    #verbose = True
    for fileName in fileNames:
      #if verbose: print "test read",fileName
      fen = PandasTabWidget()
      try:
        fen.read_oscar_csv(fileName)
      except:
        print("problem reading",fileName)
      fen.setCurrentIndex(1)
      fen.show()
      self.launchTimer(fen)
      app.exec_()

  def test_300(self):
    import math
    d1 = OrderedDict( ( # not a simple dict, have to be ordered for python 2
          ('one', [1, 2, 3, 4]),
          ('two', [2., 2.1, 2.2, 2.3]),
          ('three', [3.0, 3.1, 3.2, None]),
          ('four', [4.0, 4.1, 4.2, 4.3]),) )
    df1 = pd.DataFrame(d1)
    
    icols = [1,2]
    names = df1.columns.values.tolist()

    titles_xy = [names[icol] for icol in icols]
    #print df1['one'].values.tolist()
    #[1, 2, 3, 4]
    xy = [df1[name].values.tolist() for name in titles_xy]
    if aDebug: 
      print("titles_xy %s" % titles_xy)
      print("test_300 d1:\n%s" % d1)
      print("test_300 df1:\n%s" % df1)
      print(xy)

    self.assertEqual(xy[0][0], 2.0)
    self.assertEqual(xy[1][1], 3.1)
    self.assertTrue( math.isnan(xy[1][3]) )

  def test_999(self):
    #kill timers
    for i in range(len(timers)):
      timers.pop() 

if __name__ == '__main__':
  aDebug = True
  user = os.getenv("USERNAME")
  unittest.main()
  
