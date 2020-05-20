# TODOs
* ./pyEPR/ansys.py
  * LINE 46:	 : Replace `win32com` with Linux compatible package.
  * LINE 795:	 : check if variable does not exist and quit if it doesn't?
  * LINE 1857: : make mesh tis own  class with preperties
  * LINE 1980: : create Wirebond class
  * LINE 2012: : Add option to modify these
  * LINE 2376: : Add a rotated rectangle object.
  * LINE 2425: : points = collection of points
  * LINE 2428: : find the plane of the polyline for now, assume Z
  * LINE 2453: 
  * LINE 2503: : find the plane of the polyline for now, assume Z

* ./pyEPR/ansys.py
  * LINE 46:	 : Replace `win32com` with Linux compatible package.
  * LINE 795:	 : check if variable does not exist and quit if it doesn't?
  * LINE 1857: : make mesh tis own  class with preperties
  * LINE 1980: : create Wirebond class
  * LINE 2012: : Add option to modify these
  * LINE 2376: : Add a rotated rectangle object.
  * LINE 2425: : points = collection of points
  * LINE 2428: : find the plane of the polyline for now, assume Z
  * LINE 2453: 
  * LINE 2503: : find the plane of the polyline for now, assume Z

* ./pyEPR/core_distributed_analysis.py
  * LINE 149: 	: turn into base class shared with analysis!
  * LINE 253: 	: replace this method with the one below, here because osme funcs use it still
  * LINE 339: 	: maybe sort column and index? # todo: maybe generalize
  * LINE 488: 	: change to integer?
  * LINE 548: 	: These should be common function to the analysis and here!
  * LINE 849: 	: Update make p saved sep. and get Q for diff materials, indep. specify in pinfo
  * LINE 1046:	: maybe load from data_file
  * LINE 1064:	
  * LINE 1139:	: Move inside of loop to funciton calle self.analyze_variation
  * LINE 1247:	: this should really be passed as argument  to the functions rather than a
  * LINE 1340:	: THis need to be changed, wont work in the future with updating result etc.
  * LINE 1513:	: Move to class for reporter ?

* ./pyEPR/core_quantum_analysis.py
  * LINE 130:	: remove all copies of same data
  * LINE 574:	: superseed by Convert.ZPF_from_EPR
  * LINE 607:	: avoide analyzing a previously analyzed variation
  * LINE 741:	: actually make into dataframe with mode labela and junction labels
  * LINE 782:	: ?
  * LINE 825:	: shouldmove these kwargs to the config

* ./pyEPR/project_info.py
  * LINE 134:	: introduce modal labels
