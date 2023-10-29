#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
examples_pandas.py
"""

import os
import sys
from PyQt5 import QtGui, QtCore, QtWidgets

import pandas as pd
import pandas.tools.plotting as pdp
import numpy as np

import matplotlib as mpl
#print mpl.matplotlib_fname()
import matplotlib.pyplot as plt
import pandaspy.uranieTools as URAT


def catFile(name):
  with open(name, 'r') as f: contents = f.read()
  print("\n%s:\n%s\n\n" % (name, contents))
  return

def test_010():
  ts = pd.Series(np.random.randn(1000), index=pd.date_range('1/1/2000', periods=1000))
  ts = ts.cumsum()
  #with pd.plot_params.use('x_compat', True): #ImportError: No module named openpyxl
  #plt.xticks(rotation=60)
  ts.plot()

  df = pd.DataFrame(np.random.randn(1000, 4), index=ts.index, columns=['A', 'B', 'C', 'D'])
  df = df.cumsum()
  df.plot()

  name = "test_pandas_1_tmp.hdf5"
  #df.to_hdf(name, 'df')
  #pd.read_hdf(name, 'df')

  name = name.replace("hdf5", "xlsx")
  df.to_excel(name, sheet_name='Sheet1')
  #os.system("ooffice %s" % name)

def test_020():
  df = pd.DataFrame(np.random.rand(10, 4), columns=['a', 'b', 'c', 'd'])
  #df.plot(kind='area')
  #df.plot(kind='area', stacked=False)
  df = pd.DataFrame(np.random.randn(1000, 4), columns=['a', 'b', 'c', 'd'])
  pdp.scatter_matrix(df, alpha=0.2, figsize=(6, 6), diagonal='kde')

def test_022():
  df = pd.DataFrame(np.random.randn(1000, 4), columns=['a', 'b', 'c', 'd'])
  pdp.scatter_matrix(df, alpha=0.2, figsize=(6, 6), diagonal='kde')

def test_030():
  ser = pd.Series(np.random.randn(1000))
  print(ser.plot(kind='kde'))
  #print ser.plot()

def test_040():
  data = pd.read_csv('iris.data')
  plt.figure()
  print(pdp.andrews_curves(data, 'Name'))
  
def test_050():
  data = pd.read_csv('iris.data')
  plt.figure()
  print(pdp.parallel_coordinates(data, 'Name'))

def test_051():
  #http://pandas.pydata.org/pandas-docs/stable/timeseries.html
  dates = [np.datetime64("2016-02-01 01:02:03"), np.datetime64("2016-03-01 01:02:03")]
  dates
  ddates = pd.DatetimeIndex(dates)
  ddates.year
  ddates.month
  ddates.day
  ddates.hour
  ddates.minute
  ddates.second
  ddates.dayofyear   #The ordinal day of year
  ddates.weekofyear  #The week ordinal of the year
  ddates.week	       #The week ordinal of the year
  ddates.dayofweek   #The day of the week with Monday=0, Sunday=6
  ddates.weekday     #The day of the week with Monday=0, Sunday=6
    
  
def test_052():
  #http://pandas.pydata.org/pandas-docs/stable/generated/pandas.read_csv.html
  dtype = {"Date": np.datetime64} #, ‘Heures’: np.} 
  parse_dates = [["Date", "Heures"]]
  #/data/tmpletr/gaudier/iteSEE/
  df = pd.read_csv('eCO2mix_RTE_Annuel-Definitif_2013.csv', sep=";", engine="c", parse_dates=parse_dates )
  print("\n**types\n%s" % df.dtypes)
  print("\n**head\n%s" % df.head())
  print("\n**tail\n%s" % df.tail(3))
  print("\n**describe\n%s" % df.describe())
  print("\n**slicing\n%s" % df[:8])
  plt.figure()
  print("pdp.__file__ %s" % pdp.__file__)
  pdp.parallel_coordinates(df[:8], 'Date_Heures')

def test_054():
  #http://pandas.pydata.org/pandas-docs/stable/generated/pandas.read_csv.html
  dtype = {"Date": np.datetime64} #, ‘Heures’: np.} 
  parse_dates = [["Date", "Heures"]]
  #/data/tmpletr/gaudier/iteSEE/
  print("reading csv")
  df = pd.read_csv('eCO2mix_RTE_Annuel-Definitif_2013.csv', sep=";", engine="c", parse_dates=parse_dates )
  print("reading csv done")
  plt.figure()
  ddates = pd.DatetimeIndex(df['Date_Heures'])
  df["month"] = ddates.month
  df["dayofmonth"] = ddates.day
  df["hour"] = ddates.hour
  df["dayofweek"] = ddates.dayofweek
  URAT.parallel_coordinates(df[:2000], 'Date_Heures')

def test_060():
  plt.figure()
  data = pd.Series(0.1 * np.random.rand(1000) + 0.9 * np.sin(np.linspace(-99 * np.pi, 99 * np.pi, num=1000)))
  pdp.lag_plot(data)

def test_070():
  fig, ax = plt.subplots(1, 1)
  df = pd.DataFrame(np.random.rand(5, 3), columns=['a', 'b', 'c'])
  ax.get_xaxis().set_visible(False)   # Hide Ticks
  df.plot(table=True, ax=ax)

def test_080():
  fig, ax = plt.subplots(1, 1)
  df = pd.DataFrame(np.random.rand(5, 3), columns=['a', 'b', 'c'])
  pdp.table(ax, np.round(df.describe(), 2), loc='upper right', colWidths=[0.2, 0.2, 0.2])
  df.plot(ax=ax, ylim=(0, 2), legend=None)

if __name__ == '__main__':
  try:
    mpl.style.use('ggplot')
  except:
    print("problem mpl.style.use('ggplot')")
  plt.ion()
  
  print("matplotlib rc.name: %s" % mpl.matplotlib_fname())
  print("matplotlib backend: %s" % mpl.rcParams['backend'])
  app = QtWidgets.QApplication(sys.argv)
  #test_010()
  test_020()
  #test_022()
  #test_030()
  #test_040()
  #test_050()
  test_052()
  test_054()
  #test_060()
  #test_070()
  #test_080()
  #mw = PandasMainWidget()
  #mw.show()
  app.exec_()
  print("END of %s" % __file__)

