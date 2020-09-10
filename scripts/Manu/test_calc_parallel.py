# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 14:10:08 2020

@author: eflurin
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 10:03:32 2020

@author: eflurin
"""
from pyEPR import *
from pyEPR.ansys import Optimetrics, HfssDesign
import numpy as np
import time
import matplotlib.pyplot as plt
import scipy.optimize as sp
# Import sphere function as objective function
# Import PySwarms
import pyswarms as ps

def timestamp_name(name):
    secondsSinceEpoch = time.time()
    timeObj = time.localtime(secondsSinceEpoch)
    timestamp = '%d%d%d_%d%d%d' % (timeObj.tm_year,timeObj.tm_mon,timeObj.tm_mday,   timeObj.tm_hour, timeObj.tm_min, timeObj.tm_sec)
    return timestamp+'_'+name



################# 1.  Project and design. Open link to HFSS controls.
project_info = ProjectInfo(r'C:\HFSS_simu\\',
			     project_name = 'SMPD2', # Project file name (string). "None" will get the current active one.
			     design_name  = 'transmon_3D'       # Design name (string). "None" will get the current active one.
			    )




################# 2a. Junctions. Specify junctions in HFSS model
project_info.junctions['jtransmon'] = {'Lj_variable':'Jinduc', 'rect':'qubit_junction', 'line': 'qubit_junction_line', 'length':5e-6}
#

epr_hfss = DistributedAnalysis(project_info)

try:
    epr_hfss.save_calc_energy_electric()
    epr_hfss.save_calc_energy_magnetic()
except:
    print("expression already exist in the stack")


opti=Optimetrics( epr_hfss.design)
if 0:
    opti.import_setup("test2",r"C:\GitHub\pyEPR\scripts\Manu\202093_213743_parametric.txt")

energies=[]
for mode in range(2):
    epr_hfss.solutions.set_mode(mode+1, 0)
    opti.solve_setup("202099_102153_parametric")
    energies.append(opti.get_calc())
energies=np.array(energies)


#epr_hfss.setup.get_fields()
opti.solve_setup("202099_102153_parametric")
opti.get_calc()


self.solutions.set_mode(mode+1, 0)
        self.fields = self.setup.get_fields()


from pyEPR.ansys import Optimetrics, HfssDesign
