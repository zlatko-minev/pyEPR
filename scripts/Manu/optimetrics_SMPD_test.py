## -*- coding: utf-8 -*-
#"""
#Created on Wed Jul 15 10:03:32 2020
#
#@author: eflurin
#"""
#from pyEPR import *
#from pyEPR.ansys import Optimetrics
#import numpy as np
#import time
#
## 1.  Project and design. Open link to HFSS controls.
#project_info = ProjectInfo(r'C:\HFSS_simu\\',
#			     project_name = 'SMPD', # Project file name (string). "None" will get the current active one.
#			     design_name  = 'fullcircuit_final1'       # Design name (string). "None" will get the current active one.
#			    )
#
#
#
#
## 2a. Junctions. Specify junctions in HFSS model
#project_info.junctions['jtransmon'] = {'Lj_variable':'Jinduc', 'rect':'qubit_junction', 'line': 'qubit_junction_line', 'length':5e-6}
##
#
#project_info.ports['Waste'] = {'rect':'WasP_connector_ohm', 'R': 50, 'line': 'WasP_connector_line'}
#project_info.ports['Buffer'] = {'rect':'BufP_connector_ohm', 'R': 50,  'line': 'BufP_connector_line'}
#project_info.ports['Qubit'] = {'rect':'qubit_connector_ohm', 'R': 50, 'line': 'qubit_connector_line'}
#project_info.ports['Flux'] = {'rect':'BufR_squid_loop_lumped', 'R': 50, 'line': 'BufR_squid_connector_line'}
#
#
#
### 2b. Dissipative elements.
#project_info.dissipative['dielectrics_bulk']    = ['silicon_substrate']    # supply names here, there are more options in project_info.dissipative.
#project_info.dissipative['dielectric_surfaces'] = ['silicon_surface']
#
## 3.  Run analysis
#epr_hfss = DistributedAnalysis(project_info)
#opti=Optimetrics(epr_hfss.design)
#
#
#secondsSinceEpoch = time.time()
#timeObj = time.localtime(secondsSinceEpoch)
#timestamp = '%d%d%d_%d%d%d' % (timeObj.tm_year,timeObj.tm_mon,timeObj.tm_mday,   timeObj.tm_hour, timeObj.tm_min, timeObj.tm_sec)
#parametric_name=timestamp+'_parametric'
#
#
#epr_hfss.design.get_variables()
#
#arr=np.array([["pos_end_cable","BufR_PB_coupling"],['0.1mil', '2mm'],['0.2mil', '3mm']])
#np.savetxt('%s.txt'%parametric_name, arr, fmt='%s', delimiter='; ', newline='\n', header='', footer='', comments='# ', encoding=None)
#
#opti.import_setup(parametric_name,r"C:\GitHub\pyEPR\scripts\Manu\%s.txt"%parametric_name)
#
##opti.edit_setup(parametric_name)
#
#opti.solve_setup(parametric_name)
#
#
## 4.  Hamiltonian analysis
#epr_hfss.do_EPR_analysis()
#epr = QuantumAnalysis(epr_hfss.data_filename)
#epr.analyze_all_variations(cos_trunc = 5, fock_trunc = 4)



freqs=np.array(epr.get_frequencies()).T

nb_mode=np.array(freqs).shape[1]
nb_var=np.array(freqs).shape[0]
chis=np.array(epr.get_chis()).reshape(nb_var,nb_mode,nb_mode)

var=1

freq=freqs[var]
anharmonicity=np.abs(np.diag(chis[var]))

total_Q_from_couplings = 1/(1/np.array(epr.Qm_coupling[str(0)])).sum(1)
total_Q_from_HFSS = np.array(epr.Qs)[:,0]
Q_couplings_adjusted=np.array([total_Q_from_HFSS/total_Q_from_couplings]).T*np.array(epr.Qm_coupling[str(var)])

Waste_port_index=0
Buffer_port_index=1
##sorting modes
index={}
index['qubit']=np.argmax(anharmonicity)
index['WasP']=np.argsort(Q_couplings_adjusted[Waste_port_index])[0]
index['WasR']=np.argsort(Q_couplings_adjusted[Waste_port_index])[1]
index['BufP']=np.argsort(Q_couplings_adjusted[Buffer_port_index])[0]
index['BufR']=np.argsort(Q_couplings_adjusted[Buffer_port_index])[1]

dispersiveshifts=np.array(epr.get_chis()).reshape(nb_var,nb_mode,nb_mode)[var,index['qubit']]

computed_val={}
target_val={}
weigth={}

computed_val['qubit_anharmonicity']=anharmonicity[index['qubit']]
computed_val['WasR_DS'] = dispersiveshifts[index['WasR']]
computed_val['BufR_DS'] = dispersiveshifts[index['BufR']]
computed_val['WasP_Q'] = total_Q_from_HFSS[index['WasP']]
computed_val['BufP_Q'] = total_Q_from_HFSS[index['BufP']]
computed_val['WasR_Q'] = total_Q_from_HFSS[index['WasR']]
computed_val['BufR_Q'] = total_Q_from_HFSS[index['BufR']]
computed_val['Qubit_Q'] = total_Q_from_HFSS[index['qubit']]
computed_val['Freq_WasP']=freq[index['WasP']]
computed_val['Freq_BufP']=freq[index['BufP']]
computed_val['Freq_WasR']=freq[index['WasR']]
computed_val['Freq_BufR']=freq[index['BufR']]

target_val['qubit_anharmonicity']=0.180
target_val['BufR_DS'] = 3e-3
target_val['WasR_DS'] = 4e-3
target_val['WasP_Q'] = 30
target_val['BufP_Q'] = 30
target_val['WasR_Q'] = 1e4
target_val['BufR_Q'] = 1e4
target_val['Qubit_Q'] = 1e7
target_val['Freq_WasP']=7300
target_val['Freq_WasR']=7300
target_val['Freq_BufP']=8000
target_val['Freq_BufR']=8000

weigth['qubit_anharmonicity']=1
weigth['BufR_DS'] = 1
weigth['WasR_DS'] = 1
weigth['WasP_Q'] = 1
weigth['BufP_Q'] = 1
weigth['WasR_Q'] = 1
weigth['BufR_Q'] = 1
weigth['Qubit_Q'] = 0
weigth['Freq_WasP']=10
weigth['Freq_WasR']=10
weigth['Freq_BufP']=10
weigth['Freq_BufR']=10
   
loss=0 
for key in target_val.keys():
    loss+=weigth[key]*((computed_val[key]-target_val[key])/target_val[key])**2
    

    
