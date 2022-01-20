
.. include:: ./rst_prolog.rst


.. _iraBackground:

Background
-----------

When you choose to launch a simulation in background you can have more control on your case.

Launch in Background was thought to have the possiblity to replay case whenever you want in addition to not blocking the GUI when a simulation is launch.
Everything needed to relaunch the case is written in the output directory. And you can choose to copy your mesh here too so you can just zip the directory and send it for exemple.

Every case launched generate a save of the Solverlab tree in ..xml. It means your always have a readable save of all your parameters.
You can also choose to save the mesh locally to always have it near the parameters.

One Bash script is generated. The scrpit generated will save the latest path use for PACKAGESPY.
Without it you can't relaunch easily your case.
One python file is copied from PACKAGESPY. It's needed to launch the simulation again because it need the Python class of your Solverlab Tree.
This mean you can modify it locally and launch your simulation with almost similar parameters.

If environment variable of Solverlab and Packagespy are not set the script will not work. You need to source SOLVERLAB_install/env_solverlab.sh then ./LaunchSolverlab.sh.

