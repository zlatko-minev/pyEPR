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
    print('calculations already in the stack')



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




def jac(fxepsilon_plus,fxepsilon_minus,fx,epsilon):
    return (fxepsilon_plus-fxepsilon_minus)/(2*epsilon)


def loss_f_and_g(x0):
    
    print('current_pos =', x0)
    
    ################# 0 - define the variable position vector to be computed for evaluating the jacobian
    ##### the epsilon vector is determined based on the bounds (to be refined), note that the gradient direction is chosen randomly
    bounds_span=max_bound-min_bound
    epsilon=bounds_span/30
    x_grad_plus=x0+epsilon
    x_grad_minus=x0-epsilon
    
    x=np.array([x0]*(2*len(x0)+1))
    
    for i in range(len(x0)):
        ##### check if the gradient positions are within the bounds, otherwize we take the oposite direction
#        if x_grad[i]>max_bound[i] or x_grad[i]<min_bound[i]:
#            epsilon[i]=-epsilon[i]
#            x_grad[i]=x0[i]+epsilon[i]
            
        x[2*i+1,i]=x_grad_plus[i]
        x[2*i+2,i]=x_grad_minus[i]
    
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
    computed_val_list= []
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
        precision={}
        
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
#        target_val['Freq_qubit']= 6300.

        
#        weigth['qubit_anharmonicity']=1
#        weigth['cav_DS'] = 1
#        weigth['cav_Q'] = 
#        weigth['Freq_qubit']= 0
#        weigth['Freq_cav']= 1
        
        precision['qubit_anharmonicity']=5
        precision['cav_DS'] = 0.2
        precision['cav_Q'] = 200
#        precision['Freq_qubit']= 10
        precision['Freq_cav']= 30

        loss=0
        m=0
        for key in target_val.keys():
            print("loss for", key)
            print((computed_val[key]-target_val[key])/precision[key])
            loss+=((computed_val[key]-target_val[key])/precision[key])**2
            m+=1
        loss_allvar.append(loss/m)
        
        computed_val_list.append(computed_val)

        print(computed_val_list)

#        loss=0 
#        for key in target_val.keys():
#            print((computed_val[key]-target_val[key])/target_val[key])
#            loss+=((weigth[key]*(computed_val[key]-target_val[key])/target_val[key])**2)
#        loss_allvar.append(loss)
        
    f=np.array(loss_allvar)


    ################# 6 - compute and return the jacobian based on HFSS evals

    fxepsilon_plus=f[1::2]
    fxepsilon_minus=f[2::2]
    fx=f[0]
    
    jac_fx=jac(fxepsilon_plus,fxepsilon_minus,fx,epsilon)
    print('f=',f)
    print('jac=',jac_fx)
    np.save(r"C:\GitHub\pyEPR\scripts\Manu\%s_summary"%parametric_name,{'x0':x0,'x':x,'score':f,'jac':jac_fx,'values':computed_val_list})


    
    return fx, jac_fx



def adam(loss_f_and_g,x0,bounds,eta,beta1=1-0.3,beta2=1-0.03**2,eps=0.01):
    fx_list=[]
    x_list=[]
    
    x=x0
    m=np.zeros(len(x0))
    v=np.zeros(len(x0))


    bounds_min=np.array(bounds)[:,0]
    bounds_max=np.array(bounds)[:,1]
    eta_unit=(bounds_max-bounds_min)*eta

    start_time=time.time()
    
    for i in range(41):
        start_steptime=time.time()
        print("step number #",i)
        fx,jac_fx=loss_f_and_g(x)
        fx_list.append(fx)
        x_list.append(x)
#        jac_fx+=np.random.randn(len(jac_fx))

        m=beta1*m+(1-beta1)*jac_fx
        v=beta2*v+(1-beta2)*jac_fx**2

        m_un=m/(1-beta1**(i+1))
        v_un=v/(1-beta2**(i+1))
        x=x-eta_unit*m_un/(np.sqrt(v_un)+eps)

        print('x=',x)
        for k in range(len(x)):
            if x[k]>bounds_max[k]:
                print("x hit the bound")
                x[k]=bounds_max[k]
                m[k]=0
                v[k]=0
            if x[k]<bounds_min[k]:
                print("x hit the bound")
                x[k]=bounds_min[k]
                m[k]=0
                v[k]=0
                
        now=time.time()
        print("enlapse time for step #%i = %i min"%(i,(now-start_steptime)/60))
        print("enlapse time since start = %i min"%((now-start_time)/60))
    
        plt.figure()
        plt.plot(fx_list)
        plt.ylabel('loss')
        plt.xlabel('step #')
        plt.show()

    return x_list,fx_list

######### Main code

# Create bounds for each variable (to be determined on physical and geometrical criterion within HFSS)

name=np.array(["connect_penetrationlength1","pad_length","pad_width","Jinduc","box_height","pad_spacing"])

min_bound = np.array([-1.,0.1,0.1,5.,15.,0.05])
max_bound = np.array([3.,1.,2.,15.,30.,0.5])
bounds=[(i,j) for i,j in zip(min_bound,max_bound)]

######### Defines the optimizer sequence

######### position found by the Particle Swarm Optimizer


x0=np.array( [ 0.91755446,  0.3519256,   0.40556735,  9.34107044, 24.97869501,  0.16953198])
x0+=np.random.rand(len(x0))/200

x,fx_list=adam(loss_f_and_g,x0,bounds,eta=0.01)










