
.. include:: ../rst_prolog.rst


.. _iraBackground:

====================
Launch in background
====================

When you choose to launch a simulation in background you can have more control on your case.

Launch in Background was thought to have the possiblity to replay case whenever you want in addition to not blocking the GUI_ when a simulation is launch.
Everything needed to relaunch the case is written in the output directory. And you can choose to copy your mesh here too so you can just zip the directory and send it for exemple.

Every case launched generate a save of the Solverlab tree in a .xml file. It means you will always have a readable save of all your simulation parameters.
You can also choose to save the mesh locally to always have it near the parameters.

One Bash script is generated. The scrpit generated will save the latest path used for PACKAGESPY_.
Without it you can't relaunch easily your case.
One python file is copied from PACKAGESPY_. It is needed to launch the simulation again because it need the Python class of your SOLVERLAB_ Tree.
This mean you can modify it locally and launch your simulation with almost similar parameters.

If the environment variable of SOLVERLAB_ is not set the script will not work. You need to first locate and source the environment file in a terminal

.. code-block:: bash

    source SOLVERLAB_install/env_solverlab.sh 

then launch the bash script in background

.. code-block:: bash

    ./LaunchSolverlab.sh.

