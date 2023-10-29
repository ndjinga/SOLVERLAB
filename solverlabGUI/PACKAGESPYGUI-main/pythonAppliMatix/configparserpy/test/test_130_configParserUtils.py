#!/usr/bin/env python
# -*- coding:utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


"""
utilities for best use ConfigParser

| see:
| https://wiki.python.org/moin/ConfigParserExamples
| https://docs.python.org/2/library/configparser.html
"""

import os
import sys
import unittest
import pprint as PP
import configparserpy.configParserUtils as CPAU

# import initializeTest  # set PATH etc for test

import debogpy.debug as DBG  # Easy print stderr (for DEBUG only)
import xyzpy.stringIO as IOX


verbose = False

_examples= {
  "1_OK": u"""\
# a big comment
[SectionOne]
Status: Single # is NOT another comment
Name: Derek
Value: Yes
Age: 30
Single: True

[SectionTwo]
FavoriteColor = Green
[SectionThree]
FamilyName: Johnson

[Others]
Route: 66
""",

  "2_OK": u"""\
[SectionOne]
Param1: Hello
Param2: World

[SectionTwo]
Param1: ${SectionOne:Param1} ${SectionOne:Param2}

[SectionThree]
Alpha: One
Bravo: Two
Charlie: ${Alpha} Mississippi
""",

  "3_OK": u"""\
[SectionOne]
Param1: AnotherHello
Param2: World
Param3: Bonjour
""",

  "4_KO": u"""\
[SectionOne]
Param1: AnotherHello
Param2: World
Param1: Oops Bonjour # !!!DUPLICATE Param1
""",

  "5_default_OK": u"""\
[SectionOne]
Param1: HelloDefault
Param2: WorldDefault
Param3: BonjourDefault

[SectionTwo]
Param1: HelloDefault
Param2: WorldDefault
Param3: BonjourDefault
  """,

  "6_user_OK": u"""\
[SectionOne]
Param1: HelloUserSectionOne

[SectionTwo]
Param1: HelloUserSectionTwo

[SectionThree]
Param1: ${SectionOne:Param1}
""",

}

_res_default_user_OK = """\
UtSafeConfigParser(
 'SectionOne': {'param1': 'HelloUserSectionOne',
                'param2': 'WorldDefault',
                'param3': 'BonjourDefault'},
 'SectionThree': {'param1': 'HelloUserSectionOne'},
 'SectionTwo': {'param1': 'HelloUserSectionTwo',
                'param2': 'WorldDefault',
                'param3': 'BonjourDefault'}
)
"""

class TestCase(unittest.TestCase):
  """Test the configParserUtils.py"""

  def test_000(self):
    # one shot setUp() for this TestCase
    DBG.push_debug(verbose)
    DBG.write("CPAU __file__", CPAU.__file__)
    DBG.write("dir CPAU", dir(CPAU))
    return

  def test_999(self):
    # one shot tearDown() for this TestCase
    # SAT.setLocale() # end test english
    DBG.pop_debug()
    return

  def test_010(self):
    for iexa in sorted([k for k in list(_examples.keys()) if "OK" in k]):
      cfg = CPAU.UtSafeConfigParser()
      name = "example_%s" % iexa
      stream = IOX.StringIO(_examples[iexa])
      cfg.readfp(stream, name)
      DBG.write("cfg %s" % name, cfg)
      self.assertIn("SectionOne", cfg.sections())
    # DBG.write("dir ConfigParser", dir(cfg))

  def test_020(self):
    if CPAU._ExtendedInterpolation is None: return # old version configparser
    for iexa in sorted([k for k in list(_examples.keys()) if "KO" in k]):
      cfg = CPAU.UtSafeConfigParser()
      name = "example_%s" % iexa
      stream = IOX.StringIO(_examples[iexa])
      with self.assertRaises(Exception):
        cfg.readfp(stream, name)
        DBG.write("cfg %s" % name, cfg)

  def test_030(self):
    # test merge user in defaults as user preferences
    if CPAU._ExtendedInterpolation is None: return # old version configparser
    cfg = CPAU.UtSafeConfigParser()
    iexa = "5_default_OK"
    name = "example_%s" % iexa
    stream = IOX.StringIO(_examples[iexa])
    cfg.readfp(stream, name)

    iexa = "6_user_OK"
    name = "example_%s" % iexa
    stream = IOX.StringIO(_examples[iexa])
    cfg.readfp(stream, name)

    DBG.write("cfg merged %s" % name, cfg)

    res = _res_default_user_OK

    #print("\nxxx %s" % cfg.__repr__())
    #print("\nyyy %s" % res)
    # shit unicode 2.7 -> 3.5
    self.assertEqual(cfg.__repr__().replace("u'","'").split(), res.replace("u'","'").split())

  def test_032(self):
    if CPAU._ExtendedInterpolation is None: return # old version configparser
    # test merge user in defaults as user preferences
    cfg = CPAU.UtSafeConfigParser()
    iexadefa = "5_default_OK"
    iexauser = "6_user_OK"
    cfg = CPAU.UtSafeConfigParser()
    cfg.readFromStr(_examples[iexadefa])
    cfg.readFromStr(_examples[iexauser], merge=True)

    DBG.write("cfg merged 32", cfg)

    res = _res_default_user_OK
    self.assertEqual(cfg.__repr__().replace("u'","'").split(), res.replace("u'","'").split())

  def test_034(self):
    if CPAU._ExtendedInterpolation is None: return # old version configparser
    # test merge user in defaults as user preferences
    cfg = CPAU.UtSafeConfigParser()
    iexadefa = "5_default_OK"
    iexauser = "6_user_OK"
    cfg = CPAU.UtSafeConfigParser()
    cfg.readDefaultAndUser(_examples[iexadefa], _examples[iexauser])

    DBG.write("cfg merged 34", cfg)

    res = _res_default_user_OK
    self.assertEqual(cfg.__repr__().replace("u'","'").split(), res.replace("u'","'").split())

  def test_100(self):
    if CPAU._ExtendedInterpolation is None: return # old version configparser
    # permits class attribute writings,
    # (but raise on accentuation and avoid spaces in section names)
    cfg = CPAU.UtSafeConfigParser()
    cfg.readFromStr(u'''\
[General]
reporter = tintin

[Se rie]
titre = tintin et milou etc...
    ''')
    config = cfg.toCatchAll()
    DBG.write("cfg 100", cfg)
    DBG.write("cfg 100", str(cfg.writeToStr()))
    DBG.write("config 100", config)
    self.assertEqual(config.General.reporter, "tintin")
    self.assertEqual(config.Se_rie.titre, "tintin et milou etc...")



if __name__ == '__main__':
  unittest.main(exit=False)
  pass

