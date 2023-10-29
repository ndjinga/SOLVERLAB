#!/usr/bin/env bash

echo "//
// to launch 2 independent processes with pipe communication
//
// first process write_pipefifo_001.py will send '! HIIYOO !'
//"
./write_pipefifo_001.py '! HIIYOO !' &

sleep 1

ls -alt pipe* 
#created  prw-r--r--  1 wambeke lgls    0 30 oct.  14:30 pipefifo_001
echo "//
// have you seen named pipe: 'prw-r--r-- ... pipefifo_001'
//
// second process read, GUI, and click 'launchParentFifo' to see received '! HIIYOO !'
//"
./example_pyqt_textedit_thread_named_pipe.py

echo "// 'ps -edf | grep pipe' #have to be empty..."
ps -edf | grep pipe
echo "//
// end of launchExampleNamedPipe.sh
//"

