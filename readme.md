# Original README

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

----------



This document may be found at:

http://software.seg.org/disclaimer.txt

MAKE SURE YOU READ AND UNDERSTAND THIS DOCUMENT FULLY.
By using material from the software.seg.org GEOPHYSICS
source-code archive, you agree to abide by these terms.

TERMS OF USAGE.
---------------
The Society of Exploration Geophysicists (SEG) owns worldwide copyright
to all material posted on the software.seg.org GEOPHYSICS source-code
archive unless the documentation associated with portions of it
specifically states otherwise. All material to which SEG owns copyright
may be used freely under the following terms of use:

1) If the software is used in the creation of a publication or software
   distribution, or in any substantive way in a formal presentation, the
   GEOPHYSICS paper associated with the software must be referenced, and
   the URL for the code at software.seg.org must be referenced.

2) Notices in the source code indicating the origin of the code may not be
   removed. Disclaimers, copyright, and patent notices in the code may not
   be removed.
   Source-code files should contain:
	A one-line SEG copyright notice.
	A link to this disclaimer online at software.seg.org.
	A link to the URL for the source code at software.seg.org.

3) If any portion of a code from software.seg.org is imported into another
   routine in either its original form or in modified form, the block of
   code containing that material must be clearly noted and attributed.
   The attribution must include, at a minimum, a reference to the
   GEOPHYSICS paper associated with the software, the URL for the
   code at software.seg.org, and a pointer to this disclaimer
   (http://software.seg.org/disclaimer.txt). If any part of the algorithm
   is restricted in its use by patents, the patents must be referenced
   and the terms of use documented.

4) If any material from the software.seg.org GEOPHYSICS source-code
   archive is modified and then redistributed, clear notice must be
   included indicating that the code has been modified. The notice
   should include the nature of the modification, the extent of the
   modification, the purpose of the modification, when the modification
   was done, who made the modification, and the affiliation of the party
   that made the modification.
   Such notices should not be removed during successive redistributions.

NO WARRANTIES.
--------------
You accept all the materials provided "as is."
You assume all responsibility for the results or use of the materials.
SEG makes no representations or warranties of any kind including, but
not limited to, the accuracy, reliability, or usability of these materials.
Any use that you make of the materials is at your own risk.

DISCLAIMER.
-----------
SEG provides no warranties to you, expressed, implied, or statutory,
including any implied warranties of fitness for a particular purpose.

DAMAGE WAIVER.
--------------
In no event will SEG be liable for any damages, including direct,
indirect, special, incidental, or consequential damages, arising out of
anyone's use of or inability to use these materials, or any copies prepared
from these materials, even if SEG has been advised as to the possibility
of such damages.
