# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 10:03:32 2020

@author: eflurin
"""
from pyEPR import *
from pyEPR.ansys import Optimetrics, HfssDesign
import numpy as np
import time

# Import sphere function as objective function
# Import PySwarms
import pyswarms as ps
from pyswarms.utils.functions import single_obj as fx

def timestamp_name(name):
    secondsSinceEpoch = time.time()
    timeObj = time.localtime(secondsSinceEpoch)
    timestamp = '%d%d%d_%d%d%d' % (timeObj.tm_year,timeObj.tm_mon,timeObj.tm_mday,   timeObj.tm_hour, timeObj.tm_min, timeObj.tm_sec)
    return timestamp+'_'+name


################# 1.  Project and design. Open link to HFSS controls.
project_info = ProjectInfo(r'C:\HFSS_simu\\',
			     project_name = 'SMPD', # Project file name (string). "None" will get the current active one.
			     design_name  = 'transmon_3D'       # Design name (string). "None" will get the current active one.
			    )




################# 2a. Junctions. Specify junctions in HFSS model
project_info.junctions['jtransmon'] = {'Lj_variable':'Jinduc', 'rect':'qubit_junction', 'line': 'qubit_junction_line', 'length':5e-6}
#

################# 2c. Define ports.

#project_info.ports['Waste'] = {'rect':'WasP_connector_ohm', 'R': 50, 'line': 'WasP_connector_line'}
#project_info.ports['Buffer'] = {'rect':'BufP_connector_ohm', 'R': 50,  'line': 'BufP_connector_line'}
#project_info.ports['Qubit'] = {'rect':'qubit_connector_ohm', 'R': 50, 'line': 'qubit_connector_line'}
#project_info.ports['Flux'] = {'rect':'BufR_squid_loop_lumped', 'R': 50, 'line': 'BufR_squid_connector_line'}


################# 2b. Dissipative elements.
#project_info.dissipative['dielectrics_bulk']    = ['silicon_substrate']    # supply names here, there are more options in project_info.dissipative.
#project_info.dissipative['dielectric_surfaces'] = ['silicon_surface']




################# Define the loss function to be minimized
################# 1 - take an array x of variable values to inject in parametric sweep for HFSS
#################     'x' is of size (n,m) where 'n' is the number of variable and 'm' the number of variation to compute in parallel
################# 2 - run 'm' HFSS parametric variation in parallel
################# 3 - performs a pyEPR analysis on the new 'm' variations
################# 4 - identify the physical modes on physical criterions
################# 5 - compute the distance to the target for each variation
################# 6 - return loss as an array of size (m)

def loss_function(x):
    ###connect to HFSS
    epr_hfss = DistributedAnalysis(project_info)
    if 1:
        ###load the optimetrics module from the ansys package
        opti=Optimetrics( epr_hfss.design)
        
        ################# 1 - take an array x of variable values to inject in parametric sweep for HFSS
        ###load the optimetrics module from the ansys package
        parametric_name=timestamp_name('parametric')
        print(parametric_name)
        
        ###get the list of HFSS variable
        var=epr_hfss.design.get_variables()   
        
        ###list of variable to be optimized (should probably done outside)
        name=np.array(["connect_penetrationlength1","pad_length","pad_width","Jinduc","box_length"])
        
        ###dirty way to get the unit of the variable of interest (for some reason HFSS wants the same units than the project variable for the variation)
        units=[var[key][-2:] for key in name]
        
        ###dirty way to add the correct units to the 'x' array
        var_list=np.ones((x.shape[0],x.shape[1])).astype(str)
        for i in range(x.shape[0]):
            for j in range(x.shape[1]):
                var_list[i,j]=str(x[i,j])+units[j]
                
        ###create and save the array of string with the correct format to import to HFSS
        arr=np.vstack([name,var_list])
        np.savetxt(r"C:\GitHub\pyEPR\scripts\Manu\%s.txt"%parametric_name, arr, fmt='%s', delimiter='; ', newline='\n', header='', footer='', comments='# ', encoding=None)
        
        ###import the parametric setup with the list of variations (function I added to the ansys package)
        opti.import_setup(parametric_name,r"C:\GitHub\pyEPR\scripts\Manu\%s.txt"%parametric_name)
        
        ################# 2 - run 'm' HFSS parametric variation in parallel
        ###Solve the added parametric variation
        opti.solve_setup(parametric_name)
    
    ################# 3 - performs a pyEPR analysis on the new 'm' variations
    ###reload the list of the last variations
    epr_hfss = DistributedAnalysis(project_info)
    ###create the list containing the last 'm' variations
    var_list=list(epr_hfss.variations[-x.shape[0]:])
    epr_hfss.do_EPR_analysis(var_list)
    epr = QuantumAnalysis(epr_hfss.data_filename,var_list)
    epr.analyze_all_variations(var_list,cos_trunc = 5, fock_trunc = 4)
    
    ################# 4 - identify the physical modes on physical criterion
    freqs=np.array(epr.get_frequencies()).T   
    nb_mode=np.array(freqs).shape[1]
    nb_var=np.array(freqs).shape[0]
    chis=np.array(epr.get_chis()).reshape(nb_var,nb_mode,nb_mode)
    
    loss_allvar=[]
    for var in range(nb_var):
    
        ### get the frequencies of the current variation
        freq=freqs[var]
        ### get the anharmonicity of the current variation
        anharmonicity=np.abs(np.diag(chis[var]))
        ### get the Q of the current variation
        total_Q_from_HFSS = np.array(epr.Qs)[:,0]
        #total_Q_from_couplings = 1/(1/np.array(epr.Qm_coupling[str(0)])).sum(1)
        #Q_couplings_adjusted=np.array([total_Q_from_HFSS/total_Q_from_couplings]).T*np.array(epr.Qm_coupling[str(var)])
        
        

        ### sorting modes  
        ### define the qubit as the mode with the largest anharmanocity
        index={}
        index['qubit']=np.argsort(anharmonicity)[0]
        index['cav']=np.argsort(anharmonicity)[1]
        
        ### get the dispersive shifts of the current variation
        dispersiveshifts=np.array(epr.get_chis()).reshape(nb_var,nb_mode,nb_mode)[var,index['qubit']]
        
        
        ################# 5 - compute the distance to the target for each variation
        computed_val={}
        target_val={}
        weigth={}
        
        computed_val['qubit_anharmonicity']=anharmonicity[index['qubit']]
        computed_val['cav_DS'] = dispersiveshifts[index['cav']]
        computed_val['Freq_cav']=freq[index['cav']]
        computed_val['cav_Q'] = total_Q_from_HFSS[index['cav']]
        
#        def supfunc(x,delta):
#            return x*(x>delta)
#                
        ####definition of the targets (should probably done outside)
        target_val['qubit_anharmonicity']=180
        target_val['cav_DS'] = 3
        target_val['Freq_cav']= 7300
        target_val['cav_Q'] = 1e4
        
        weigth['qubit_anharmonicity']=1
        weigth['cav_DS'] = 1
        weigth['Freq_cav']= 10
        weigth['cav_Q'] = 1
        
        np.save(r"C:\GitHub\pyEPR\scripts\Manu\%s"%parametric_name,computed_val)
        
        loss=0 
        for key in target_val.keys():
            loss+=weigth[key]*((computed_val[key]-target_val[key])/target_val[key])**2
        loss_allvar.append(loss)
        
    ################# 6 - return loss as an array of size (m)
    return np.array(loss_allvar)
    


######### Main code
######### Defines the optimizer sequence


# Create bounds for each variable (to be determined on physical and geometrical criterion within HFSS)
min_bound = np.array([-1,0.1,0.1,5,20])
max_bound = np.array([3,1,2,15,60])
bounds = (min_bound, max_bound)


# Initialize swarm (I have no idea what this is)
options = {'c1': 0.5, 'c2': 0.3, 'w':0.9}

# Call instance of PSO with bounds argument
# The maximum number of particle corresponds to the number of variation HFSS can perform in parallel
optimizer = ps.single.GlobalBestPSO(n_particles=10, dimensions=5, options=options, bounds=bounds)

# Perform optimization
cost, pos = optimizer.optimize(loss_function, iters=30)

