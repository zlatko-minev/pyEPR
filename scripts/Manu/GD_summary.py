# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 11:09:34 2020

@author: eflurin
"""

import numpy as np
import os, sys, re
import matplotlib.pyplot as plt

def getfiles_h5(dirpath,extention_string='.h5'):
    a = [s for s in os.listdir(dirpath)
         if os.path.isfile(os.path.join(dirpath, s))]
    a.sort(key=lambda s: os.path.getmtime(os.path.join(dirpath, s)))
    return [i for i in a if i.endswith(extention_string)]

file_list=getfiles_h5(r'C:\GitHub\pyEPR\scripts\Manu','summary.npy')
file_list=getfiles_h5(r'C:\GitHub\pyEPR\scripts\Manu\GD_20200906','summary.npy')


f_max=[]
x0_max=[]
values_max=[]
jac_list=[]
for file in file_list:
    summary=np.load('C:\GitHub\pyEPR\scripts\Manu\GD_20200906\\'+file)[()]
    
    summary['score']
    
    
    
    f_max.append(summary['score'][0])
    x0_max.append(summary['x0'])
    values_max.append(summary['values'][0])
    
    jac_list.append(summary['jac'])
print(f_max)
print(x0_max)
print(values_max)

jac_list=np.array(jac_list)
x0_max=np.array(x0_max)

qubit_anharmonicity = [values['qubit_anharmonicity'] for values in values_max]
cav_DS = [values['cav_DS'] for values in values_max]
Freq_cav = [values['Freq_cav'] for values in values_max]
cav_Q = [values['cav_Q'] for values in values_max]
Freq_qubit = [values['Freq_qubit'] for values in values_max]


plt.figure()
plt.semilogy(f_max,label='f_max')
plt.legend()

plt.figure()
plt.plot(qubit_anharmonicity,label='qubit_anharmonicity')
plt.legend()

plt.figure()
plt.plot(cav_DS,label='cav_DS')
plt.legend()

plt.figure()
plt.plot(Freq_cav,label='Freq_cav')
plt.legend()

plt.figure()
plt.plot(cav_Q,label='cav_Q')
plt.legend()

plt.figure()
plt.plot(Freq_qubit,label='Freq_qubit')
plt.legend()


name=np.array(["connect_penetrationlength1","pad_length","pad_width","Jinduc","box_height","pad_spacing"])

for i in range(6):
    plt.figure()

    plt.plot(jac_list[:,i],label='jac_list_'+name[i])
    plt.legend()
    
    plt.figure()

    plt.plot(x0_max[:,i],label='x0_max_'+name[i])
    plt.legend()


plt.figure()
plt.plot(x0_max[:,1],label='x0_max'+name[1])
plt.plot(x0_max[:,2],label='x0_max'+name[2])



plt.figure()
plt.plot(jac_list[:,1],label='x0_max'+name[1])
plt.plot(jac_list[:,2],label='x0_max'+name[2])
