#!/usr/bin/env python
#-*- coding:utf-8 -*-

"""
see https://pypi.org/project/python-fit/
python 2-3 compliant
"""

import os
import sys
import unittest
import pprint as PP

import fitpy.fit as FIT
import fitpy.test.essai_fit_1 as FIT_1

verbose = False # True
isplot = False # True

class TestCase(unittest.TestCase):
  "Test the fit.py"""

  def xtest_000(self):
    # one shot setUp() for this TestCase
    return

  def test_011(self):
    res = FIT_1.essai_1(isplot=isplot)

  def test_012(self):
    res = FIT_1.essai_2(isplot=isplot)

  def xtest_999(self):
    # one shot tearDown() for this TestCase
    return

if __name__ == '__main__':
    unittest.main(exit=False)
    pass

