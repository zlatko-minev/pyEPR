# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 10:03:32 2020

@author: eflurin
"""
from pyEPR import *
from pyEPR.ansys import Optimetrics, HfssDesign
import numpy as np
import time

import scipy.optimize as sp
# Import sphere function as objective function
# Import PySwarms
import pyswarms as ps

def timestamp_name(name):
    secondsSinceEpoch = time.time()
    timeObj = time.localtime(secondsSinceEpoch)
    timestamp = '%d%d%d_%d%d%d' % (timeObj.tm_year,timeObj.tm_mon,timeObj.tm_mday,   timeObj.tm_hour, timeObj.tm_min, timeObj.tm_sec)
    return timestamp+'_'+name

def jac(fxepsilon,fx,epsilon):
    return (fxepsilon-fx)/epsilon

################# 1.  Project and design. Open link to HFSS controls.
project_info = ProjectInfo(r'C:\HFSS_simu\\',
			     project_name = 'SMPD2', # Project file name (string). "None" will get the current active one.
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

################# 0 - define the variable position vector to be computed for evaluating the jacobian
################# 1 - take an array x of variable values to inject in parametric sweep for HFSS
#################     'x' is of size (n,m) where 'n' is the number of variable and 'm' the number of variation to compute in parallel
################# 2 - run 'm' HFSS parametric variation in parallel
################# 3 - performs a pyEPR analysis on the new 'm' variations
################# 4 - identify the physical modes on physical criterions
################# 5 - compute the distance to the target for each variation
################# 6 - compute and return the jacobian based on HFSS evals


# Create bounds for each variable (to be determined on physical and geometrical criterion within HFSS)
min_bound = np.array([-1.,0.1,0.1,5.,15.,0.05])
max_bound = np.array([3.,1.,2.,15.,30.,0.5])
bounds=[(i,j) for i,j in zip(min_bound,max_bound)]

def loss_f_and_g(x0):
    
    print('current_pos =', x0)
    
    ################# 0 - define the variable position vector to be computed for evaluating the jacobian
    ##### the epsilon vector is determined based on the bounds (to be refined), note that the gradient direction is chosen randomly
    bounds_span=max_bound-min_bound
    epsilon=bounds_span/100.*(2*np.random.randint(2,size=len(x0))-1)
    x_grad=x0+epsilon
    
    x=np.array([x0]*(len(x0)+1))
    
    for i in range(len(x0)):
        ##### check if the gradient positions are within the bounds, otherwize we take the oposite direction
        if x0[i]>max_bound[i]:
            x0[i]=max_bound[i]
            
        if x0[i]<min_bound[i]:
            x0[i]=min_bound[i]
            
        if x_grad[i]>max_bound[i] or x_grad[i]<min_bound[i]:
            epsilon[i]=-epsilon[i]
            x_grad[i]=x0[i]+epsilon[i]
            
        x[i+1,i]=x_grad[i]
    
    print('x',x)

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
        name=np.array(["connect_penetrationlength1","pad_length","pad_width","Jinduc","box_height","pad_spacing"])
        
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
    
    root_array=[]
    for var in range(nb_var):
    
        
        
        ### get the frequencies of the current variation
        chi_dic=epr.results.get_chi_O1()
        chis=np.abs(np.array(chi_dic[var_list[var]]))
        ### get the frequencies of the current variation
        freq_dic=epr.results.get_frequencies_O1()
        freq=np.abs(np.array(freq_dic[var_list[var]]))
        ### get the anharmonicity of the current variation
        anharmonicity=np.diag(chis)
        ### get the Q of the current variation
        total_Q_from_HFSS = np.array(epr.Qs)[:,var]
        #total_Q_from_couplings = 1/(1/np.array(epr.Qm_coupling[str(0)])).sum(1)
        #Q_couplings_adjusted=np.array([total_Q_from_HFSS/total_Q_from_couplings]).T*np.array(epr.Qm_coupling[str(var)])



        ### sorting modes  
        ### define the qubit as the mode with the largest anharmanocity
        index={}
        index['qubit']=np.argsort(anharmonicity)[-1]
        index['cav']=np.argsort(anharmonicity)[-2]
        
        ### get the dispersive shifts of the current variation
        dispersiveshifts=chis[index['qubit']]
        
        print('freq=',freq)
        print('anharmonicity=',anharmonicity)
        print('dispersiveshifts=',dispersiveshifts)
        print('total_Q_from_HFSS=',total_Q_from_HFSS)
        
        ################# 5 - compute the distance to the target for each variation
        computed_val={}
        target_val={}
        weigth={}
        
        computed_val['qubit_anharmonicity']=anharmonicity[index['qubit']]
        computed_val['cav_DS'] = dispersiveshifts[index['cav']]
        computed_val['Freq_cav']=freq[index['cav']]
        computed_val['cav_Q'] = total_Q_from_HFSS[index['cav']]
        computed_val['Freq_qubit']=freq[index['qubit']]
        
                
        ####definition of the targets (should probably done outside)
        target_val['qubit_anharmonicity']=180.
        target_val['cav_DS'] = 3.
        target_val['Freq_cav']= 7300.
        target_val['cav_Q'] = 1e4
        target_val['Freq_qubit']= 6300.

        
        weigth['qubit_anharmonicity']=1
        weigth['cav_DS'] = 1
        weigth['cav_Q'] = 1
        weigth['Freq_qubit']= 10
        weigth['Freq_cav']= 10

                
        np.save(r"C:\GitHub\pyEPR\scripts\Manu\%s_anh_DS_freq_Q"%parametric_name,computed_val)
        print(computed_val)
        
        k=0
        while len(computed_val)!=len(x0):
            keywordstring='dummy'+str(k)
            computed_val[keywordstring] = 1
            target_val[keywordstring] = 1
            weigth[keywordstring] = 1
            k+=1
            print(computed_val)
#        loss=0 
#        for key in target_val.keys():
#            print((computed_val[key]-target_val[key])/target_val[key])
#            loss+=weigth[key]*((computed_val[key]-target_val[key])/target_val[key])**2
#        
        root_vec=[weigth[key]*(computed_val[key]-target_val[key])/target_val[key] for key in target_val.keys()]
        root_array.append(root_vec)

    root_array=np.vstack(root_array)


    ################# 6 - compute and return the jacobian based on HFSS evals
    print('f_var=',root_array)

    root_xepsilon=root_array[1:]
    root_x=root_array[0]
    jac_root=jac(root_xepsilon,root_x,epsilon).T
    print('f=',root_x)
    print('jac_fx=',jac_root)
    np.save(r"C:\GitHub\pyEPR\scripts\Manu\%s_jac"%parametric_name,jac_root)
    np.save(r"C:\GitHub\pyEPR\scripts\Manu\%s_f"%parametric_name,root_x)


    
    return root_x, jac_root
        
  

######### Main code
######### Defines the optimizer sequence


######### position found by the Particle Swarm Optimizer
name=np.array(["connect_penetrationlength1","pad_length","pad_width","Jinduc","box_height","pad_spacing"])

x0=np.array([1.+np.random.rand()/100,  0.5,  0.5,  10.,  25., 0.15])



res=sp.root(loss_f_and_g, x0, jac=True)
print(res)


var=epr_hfss.design.get_variables()   
units=[var[key][-2:] for key in name]
for i in range(len(name)):
    epr_hfss.design.set_variable(name[i],str(res.x[i])+units[i])
