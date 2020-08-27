#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 11:33:44 2020

@author: emmanuel
"""


from bayes_opt import BayesianOptimization
from bayes_opt import SequentialDomainReductionTransformer
from bayes_opt import UtilityFunction
import numpy as np
import matplotlib.pyplot as plt

from pyEPR import *
from pyEPR.ansys import Optimetrics, HfssDesign
import time

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


def func_to_max_HFSS(x0):
    
    print('current_pos =', x0)
    

    x=np.array(x0)

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
    
    loss_allvar=[]
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
        loss=0 
        for key in target_val.keys():
            print((computed_val[key]-target_val[key])/target_val[key])
            loss+=(weigth[key]**2*(computed_val[key]-target_val[key])/target_val[key])**2
        loss_allvar.append(loss)
        
    f=-np.array(loss_allvar)




    
    print('f=',f)
    np.save(r"C:\GitHub\pyEPR\scripts\Manu\%s_f"%parametric_name,f)


    
    return f



def random_kappa(x0,sigma):
    kappa=sigma*np.random.randn()+x0
    while kappa<0.:
        kappa=sigma*np.random.randn()+x0
    return kappa


def maximize(func_to_max,guess=[],N_parallel=7,n_iter=30,x0=4,sigma=2,dropout_frac=0.2):
    points_list=[]
    targets_list=[]
    optimax=[]
    next_points=guess
    for k in range(n_iter):
        for i in range(N_parallel):
            
            optimizer = BayesianOptimization(f = func_to_max_HFSS, pbounds=pbounds, bounds_transformer=bounds_transformer)#,
            if points_list is not []:
                rand_index=np.random.permutation(len(points_list))[:int((1-dropout_frac)*len(points_list))]
                print(len(rand_index))
                for r in rand_index:
                    optimizer.register(params=points_list[r], target=targets_list[r])

            
            kappa=random_kappa(x0,sigma)
            print('random kappa =',kappa)
            utility = UtilityFunction(kind="ucb", kappa=kappa, xi=0.0)

            suggested_point=optimizer.suggest(utility)
            if suggested_point not in points_list:
                if suggested_point not in next_points:
                    next_points.append(suggested_point)


#        next_points_array=np.array([list(dic.values()) for dic in next_points])
        next_points_array=np.array([[dic[n] for n in name] for dic in next_points])
        targets=func_to_max(next_points_array)

#        for next_point,target in zip(next_points,targets):
#            optimizer.register(params=next_point, target=target)
        for target in targets:
            targets_list.append(target)
        for next_point in next_points:
            points_list.append(next_point)

        print()
        print(targets, next_points)
        next_points=[]

        print()
        print('"current_max =', optimizer.max)
        try:
            optimax.append(optimizer.max['target'])
        except:
            print('lol')

    print('result')
    print(optimizer.max)
    return optimax



min_bound = np.array([-1.,0.1,0.1,5.,15.,0.05])
max_bound = np.array([3.,1.,2.,15.,30.,0.5])
name=np.array(["connect_penetrationlength1","pad_length","pad_width","Jinduc","box_height","pad_spacing"])
pbounds = {n: (min_b,max_b) for n,min_b,max_b in zip(name,min_bound,max_bound)}


x0=np.array([1.+np.random.rand()/100,  0.5,  0.5,  10.,  25., 0.15])

guess=[{n:x for n,x in zip(name,x0)}]


bounds_transformer = SequentialDomainReductionTransformer()
## Reduces the bounds ###
optimizer = BayesianOptimization(f = func_to_max_HFSS, pbounds=pbounds)#,
                              #   bounds_transformer=bounds_transformer)



optimax=maximize(func_to_max_HFSS,guess)



















