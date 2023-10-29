#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2008-20xx  CEA/DEN
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

import unittest

import solverlabpy.abcdExpression as ABCD
from solverlabpy.abcdExpression import toAbcd, toEvalAbcd

class TestCase(unittest.TestCase):
  "Test the toAbcd"""

  def test_020(self):
    a, b, c, d = list("abcd")
    self.assertEqual(a + b + c + d, "abcd")
    self.assertEqual(3 * a + b, "aaab")

  def test_030(self):
    self.assertEqual(toAbcd("ABcD"), "a+b+c+d")
    self.assertEqual(toAbcd("abcd"), "a+b+c+d")
    self.assertEqual(toAbcd("(abcd)"), "(a+b+c+d)")
    self.assertEqual(toAbcd("ab(cd)"), "a+b+(c+d)")
    self.assertEqual(toAbcd("ab2(cd)"), "a+b+2*(c+d)")
    self.assertEqual(toAbcd("3a"), "3*a")
    self.assertEqual(toAbcd("30a"), "30*a")
    self.assertEqual(toAbcd("3abc"), "3*a+b+c")
    self.assertEqual(toAbcd("3a2bc"), "3*a+2*b+c")
    self.assertEqual(toAbcd("3a2b4c"), "3*a+2*b+4*c")
    self.assertEqual(toAbcd("3(abc)"), "3*(a+b+c)")
    self.assertEqual(toAbcd("10(abc)"), "10*(a+b+c)")
    self.assertEqual(toAbcd("10(a3(bc))"), "10*(a+3*(b+c))")

  def test_040(self):
    self.assertEqual(toEvalAbcd("ABcD"), "abcd")
    self.assertEqual(toEvalAbcd("  AB  c    D"), "abcd")
    self.assertEqual(toEvalAbcd("abcd"), "abcd")
    self.assertEqual(toEvalAbcd("(abcd)"), "abcd")
    self.assertEqual(toEvalAbcd("ab(cd)"), "abcd")
    self.assertEqual(toEvalAbcd("ab2(cd)"), "abcdcd")
    self.assertEqual(toEvalAbcd("3a"), "aaa")
    self.assertEqual(toEvalAbcd("30a"), 30 * "a")
    self.assertEqual(toEvalAbcd("3abc"), "aaabc")
    self.assertEqual(toEvalAbcd("3a2bc"), "aaabbc")
    self.assertEqual(toEvalAbcd("3a2b4c"), "aaabbcccc")
    self.assertEqual(toEvalAbcd("3(abc)"), "abcabcabc")
    self.assertEqual(toEvalAbcd("10(abc)"), 10 * "abc")
    self.assertEqual(toEvalAbcd("10(a3(bc))"), 10 * "abcbcbc")
    self.assertEqual(toEvalAbcd("ab(3((cd)))"), "abcdcdcd")
    self.assertEqual(toEvalAbcd("ab( ( 2(dcd) ) )"), "abdcddcd")

  def test_060(self):
    # to see exception message(s) or else
    with self.assertRaises(Exception): toAbcd("-abcd")
    with self.assertRaises(Exception): toAbcd("ab[cd")
    with self.assertRaises(Exception): toAbcd("10(+abc)")
    with self.assertRaises(Exception): toAbcd("a+bc")
    with self.assertRaises(Exception): toAbcd("a*bc")
    with self.assertRaises(Exception): toAbcd("a*bc")

  def test_080(self):
    with self.assertRaises(Exception): toEvalAbcd("abcd(")
    with self.assertRaises(Exception): toEvalAbcd("ab+(cd)")
    with self.assertRaises(Exception): toEvalAbcd("10(+abc)")
    with self.assertRaises(Exception): toEvalAbcd("a+bc")
    with self.assertRaises(Exception): toEvalAbcd("a*bc")
    with self.assertRaises(Exception): toEvalAbcd("a*bc")



if __name__ == '__main__':
  unittest.main(exit=False)
  pass

