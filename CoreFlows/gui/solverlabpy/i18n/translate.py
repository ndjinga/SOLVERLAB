#!/bin/env python

"""
create fr/LC_MESSAGES/solverlabGUI.mo from r/LC_MESSAGES/solverlabGUI.po
"""

import polib
txt = open('fr/LC_MESSAGES/solverlabGUI.po', 'r').read()
po = polib.pofile(pofile=txt, encoding='utf-8')
po.save_as_mofile('fr/LC_MESSAGES/solverlabGUI.mo')
print("OK translate.py seems to be done")
