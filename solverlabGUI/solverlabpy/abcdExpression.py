#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# Copyright (C) 2008-2023  CEA
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
#
# See http://www.salome-platform.org or email : webmaster.salome@opencascade.com
# %% LICENSE_END


"""
from '3(a2bc)' to 'abbcabbcabbc'
without regexp, not recursive
for not smart poor people
"""

import os
import sys

import unittest
import pprint as PP


_authorized = list("()0123456789abcdefghijklmnopqrstuvwxyz")
_alphabet = list("abcdefghijklmnopqrstuvwxyz")
_numeric = list("0123456789")


def toAbcd(aStr, verbose=False, details=False):
  """
  '10(abc)' to return  '10*(a+b+c)'
  raise exception if problem
  """
  e = ""  # result as python expression
  st = "b"  # begin status
  if verbose: print("****** toAbc('%s')" % aStr)
  for c in aStr.lower():
    if c == " ":
      continue
    if details: print("add '%s' to '%s' status '%s'" % (c, e, st))
    if c not in _authorized:
      raise Exception("forbidden character '%s' in abcdExpression '%s'#" % (c, aStr))

    if c in _alphabet:
      if st == "a":
        # st = "a"
        e += "+"
        e += c
        continue
      if st == "n":
        st = "a"
        e += "*"
        e += c
        continue
      if st == "b":
        st = "a"
        # e += "+"
        e += c
        continue

    if c in _numeric:
      if st == "a":
        st = "n"
        e += "+"
        e += c
        continue
      if st == "n":
        st = "n"
        # e += "+"
        e += c
        continue
      if st == "b":
        st = "n"
        # e += "+"
        e += c
        continue

    if c == "(":
      if st == "a":
        st = "b"
        e += "+"
        e += c
        continue
      if st == "n":
        st = "b"
        e += "*"
        e += c
        continue
      if st == "b":
        # st = "b"
        # e += "+"
        e += c
        continue

    if c == ")":
      if st == "a":
        # st = "a"
        # e += "+"
        e += c
        continue
      if st == "n":  # error syntax
        st = "b"
        e += "*"
        e += c
        continue
      if st == "b":  # error syntax
        st = "b"
        # e += "+"
        e += c
        continue

    raise Exception("unexpected character '%s' in abcdExpression '%s'" % (c, aStr))

  if verbose: print("toAbcd '%s' -> '%s' (python expression)" % (aStr, e))
  return e


def toEvalAbcd(aStr, verbose=False):
  """
  '3(a2bc)' to return  'abbcabbcabbc'
  raise exception if problem
  """
  if verbose: print("\n****** toEvalAbcd('%s')" % aStr)
  a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z = list("abcdefghijklmnopqrstuvwxyz")
  pyExpr = toAbcd(aStr)
  try:
    res = eval(pyExpr)
  except Exception as e:
    raise Exception("syntax error in abcdExpression '%s': %s" % (aStr, e))

  if verbose: print("toEvalAbc '%s' -> '%s'" % (pyExpr, res))
  return res

def getIndiceFromChar(aChar):
  """returns 0 for 'a', 1 for 'b' etc., max is 'z'"""
  lower = aChar.lower() # accept 'A' as 'a' etc.
  if lower < 'a' or lower > 'z':
    raise Exception("only 'a' to 'z' accepted in abcdExpression, got '%s'" % aChar)
  res = ord(lower) - ord('a') # ord(aChar) - 97
  return res

def toEval0123(aStr, verbose=False):
  """
  '3(a2bc)' to return [0,1,1,2,0,1,1,2,0,1,1,2] for 'abbcabbcabbc'
  raise exception if problem
  """
  if verbose: print("\n****** toEval0123('%s')" % aStr)
  res_c = toEvalAbcd(aStr, verbose=verbose)
  res = [getIndiceFromChar(c) for c in res_c]
  return res

def toEvalAbcdForTooltip(aStr, length=20):
  """
  set results in lenght characters lines
  """
  tmp = toEvalAbcd(aStr)
  nb = len(tmp)
  posdep = 0
  res = []
  for i in range(0, nb, length)[1:]:
    # print("toEvalAbcdForTooltip %i %s" % (i, tmp[posdep:i))
    res.append(tmp[posdep:i])
    posdep = i
  res.append(tmp[posdep:])
  if len(res) > 20: # lines
    res = res[0:10] + ["..."] + res[-10:]
  res = ["cell_count_x = %i" % nb] + res
  return '\n'.join(res)


