# -*- coding: latin-1 -*-
#  Copyright (C) 2007-2010  CEA/DEN, EDF R&D, OPEN CASCADE
#
#  Copyright (C) 2003-2007  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
#  CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
#  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
#
# Author : A. Bruneton
#
import sys, os
from PyQt4.QtGui import QApplication 
from PyQt4.QtCore import SIGNAL, SLOT

from CFDesktop import CFDesktop
import SalomePyQt_MockUp
from PyQt4.QtCore import QTimer, QTranslator, Qt

desktop = None

def activate():
    """This method mimicks SALOME's module activation """
    global desktop
    fv = desktop.showCentralWidget()
    return True

def main(args) :
    global desktop
      
    app = QApplication(args)
    
    desktop = CFDesktop(None)
    sgPyQt = SalomePyQt_MockUp.SalomePyQt(desktop)
    desktop._sgPyQt = sgPyQt 
    desktop.initialize()
    activate()
    desktop.show()
    
    
    #
    app.connect(app,SIGNAL("lastWindowClosed()"),app,SLOT("quit()"))
    app.exec_()

if __name__ == "__main__" :
    main(sys.argv)
