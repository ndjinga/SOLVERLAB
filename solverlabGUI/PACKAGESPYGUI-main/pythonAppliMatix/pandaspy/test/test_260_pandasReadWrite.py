#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END

"""\
test module pandas read/write on miscellanous types of files
rm -rf ~/pandasReadWriteTest_tmp*
head -40 ~/pandasReadWriteTest_tmp*
"""

import os
import unittest
import time
import pprint as PP
import xyzpy.stringIO as IOX

try:
  import pandas as pd
  pandasOk = True
except:
  print("ERROR: no pandas library, no test")
  pandasOk = False


verbose = False

#thanks http://chrisalbon.com/python/pandas_dataframe_importing_csv.html

raw_data = {
  'first_name': ['Jason', 'Molly', 'Tina', 'Jake', 'Amy'],
  'last_name': ['Miller', 'Jacobson', ".", 'Milner', 'Cooze'],
  'age': [42, 52, 36, 24, 73],
  'preTestScore': [4, 24, 31, ".", "."],
  'postTestScore': ["25,000", "94,000", 57, 62, 71.5]
  }

df = pd.DataFrame(raw_data, columns=sorted(raw_data.keys()))

pd_readers = [i for i in dir(pd) if "read_" in i[0:6] ]
pd_writers = [i for i in dir(df) if "to_" in i[0:4] ]
excludes = ["to_dense",   "to_clipboard", 
            "to_wide",    "to_dict", 
            "to_panel",   "to_period", 
            "to_sparse",  "to_excel", 
            "to_gbq",     "to_hdf", 
            "to_sql",     "to_stata", 
            "to_string",  "to_timestamp",
            "to_msgpack", "to_pickle",
            "to_records", "to_feather", 
            "to_parquet", "to_xarray",
           ]
pd_writers = [i for i in pd_writers if i not in excludes]


class TestCase(unittest.TestCase):

  testDir=os.path.join(os.getenv("HOME"))
  aTestFile = os.path.join(testDir, "pandasReadWriteTest_tmp_")

  def test_001(self):
    if verbose:
      print("pandas version: %s" % pd.__version__)
      print("testDir: %s" % self.testDir)
      #print("dir(pd):\n%s" % PP.pformat(dir(pd)))
      print("pd_readers:\n%s" % PP.pformat(pd_readers))
      print("pd_writers:\n%s" % PP.pformat(pd_writers))
      
  def test_005(self):
    df = pd.DataFrame()
    #print dir(self)
    self.assertTrue("Empty DataFrame" in df.__repr__())
    self.assertTrue("Empty DataFrame" in df.to_string())
    
    
  def test_010(self):
    for i in pd_writers:
      if verbose: print("!!!!!!!! pd_writers ?", i)
    for i in "to_csv to_html to_json to_latex".split():
      testFile = self.aTestFile + i
      cmd = "df.%s('%s')" % (i, testFile)
      if verbose: print(cmd)
      exec( cmd , {"df": df} ) # python 2-3

  def test_012(self):
    df.to_csv(self.aTestFile + "to_csv2", sep=";")
    
  def test_014(self):
    if verbose: print("to string:",df.to_string())
    
  def xtest_016(self):
    #io.excel.xlsx.writer
    #ImportError: No module named openpyxl
    #ImportError: No module named xlwt
    df.to_excel(self.aTestFile + ".xlsx")
    
  def test_100(self):
    # Writing to a buffer
    output = IOX.StringIO()
    output.write('This goes into the buffer. ')
    output.write('And so does this.')
    self.assertEqual(output.getvalue(), 'This goes into the buffer. And so does this.')
    # Retrieve the value written
    #print output.getvalue()
    output.close() # discard buffer memory

    # Initialize a read buffer
    input = IOX.StringIO('Inital value for read buffer')
    # Read from the buffer
    #print input.read()
    self.assertEqual(input.read(), 'Inital value for read buffer')
    
  def test_110(self):
    output = IOX.StringIO()
    df.to_csv(output, sep=";")
    expected=""";age;first_name;last_name;postTestScore;preTestScore
0;42;Jason;Miller;25,000;4
1;52;Molly;Jacobson;94,000;24
2;36;Tina;.;57;31
3;24;Jake;Milner;62;.
4;73;Amy;Cooze;71.5;.
"""
    self.assertEqual(output.getvalue(), expected)
    
  def test_120(self):
    output = IOX.StringIO()
    df.to_csv(output, sep=";")
    output.seek(0)
    res = pd.read_csv(output, sep=";")
    if verbose:
      print("init:\n", df.to_string())
      print("read:\n", res.to_string())
  
  def test_130(self):
    output = IOX.StringIO()
    df.to_csv(output)
    output.seek(0)
    res = pd.read_csv(output) #create Unnamed: 0
    if verbose:
      print("init:\n", df.to_string())
      print("read:\n", res.to_string())
    
    res.drop(res.columns[[0]], axis=1, inplace=True) #remove Unnamed: 0
    self.assertEqual(df.to_string(), res.to_string())
    
    output = IOX.StringIO()
    df.to_csv(output)
    output.seek(0)
    res2 = pd.read_csv(output, index_col=0) #no create Unnamed: 0
    if verbose:
      print("init:\n", res.to_string())
      print("read2:\n", res2.to_string())
    
    self.assertEqual(df.to_string(), res2.to_string())
    
  

if __name__ == '__main__':
  #verbose = True #False
  unittest.main()
  pass

 
