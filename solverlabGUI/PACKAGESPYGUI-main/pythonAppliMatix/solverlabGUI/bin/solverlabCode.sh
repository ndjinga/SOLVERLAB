#!/bin/bash

#user may customize... other releases...

echo
echo "This script is used to launch SOLVERLABGUI"
echo "current directory: "`pwd`
echo "current solverlab code executable: "`which solverlab.exe`
echo

sleep 1 #to synchronize stdout/stderr echo in qt run log window 

${SOLVERLABGUI_ROOT_DIR}/solverlabGUI --gui"

echo
echo "END of SOLVERLABGUI"
echo
