Readme File for RESINVM3D.v1		
Jan 16/06.

Copyright (c) 2007 by the Society of Exploration Geophysicists.
For more information, go to http://software.seg.org/2007/0001 .
You must read and accept usage terms at:
http://software.seg.org/disclaimer.txt before use.


RESINVM3D is a MATLAB package for inverting 3D Dc Resistivity and Electrical 
Resistivity Tomography data. The Package contains 35m files along with two 
demonstration data sets.

The demos that are included in the package serve as templates for the user to
adapt the code for their given problem. The demos, DCdriverS.m and DCdriverBH.m
run demonstrations for a surface based survey and a borehole survey,
respectively. 

These demo codes first create a parameter structure ('para') that controls the
amount of regularization, the number of iterations, convergence criteria, and
the tolerances of the internal solvers. Second, the demo codes create a data
structure ('MTX') from a series of seven user defined input files. These files
contain source and receiver locations, the measured data and errors, the
discretization of the model space, model weighting, the reference model.
In addition, one of the input files toggles certain model parameters 'on' 
or 'off'.   

Details regarding the format and naming of the input files can be found in
either of the documented demonstration files. 

Following the creation of the two structures, 'para' and 'MTX', the
demonstrations will call 'InvMainN.m' which is the inversion code. The code
will run until the program has converged or the maximum number of iterations
is reached. Once the inversion has finished, results for the demonstrations
are automatically plotted.  

When adapting this code to work with a given data set, it is suggested that
the user copy and modify one of the driver files to work with their specific
problem. Both DCdriverS.m and DCdriverBH.m are commented. The format on the
input files, and the contents of 'para' are detailed within these files. We
refer the user to these files for further information.


RUNNING THE DEMOS

To run either of the demos, first, unzip RESINVM3D. Set the MATLAB working
directory to RESINVM3D. From this directory you can call either DCdriverS.m
or DCdriverBH.m. This call can be done by typing either 'DCdriverS' or
'DCdriverBH' in the MATLAB command line. Either call uses a series of input
data files that are located in the respective directories:

         RESINVM3D/modeldata/surface/
        
         RESINVM3D/modeldata/borehole/ 
        
Note: Do not add the RESINVM3D folder and subfolders to the MATLAB path. The 
input files for DCdriverS.m and DCdriverBH.m have the same names, and thus
the code can call the incorrect files depending on the order of folders in
the path. Both of the driver files will temporarily add the appropriate
folders to the path.  
