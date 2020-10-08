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
file_list=getfiles_h5(r'C:\GitHub\pyEPR\scripts\Manu\Pymoo_sumary','summary.npy')


f_max=[]
x0_max=[]
values_max=[]
for file in file_list[1:]:
    summary=np.load('C:\GitHub\pyEPR\scripts\Manu\Pymoo_sumary\\'+file)[()]
    
    summary['score']
    
    
    f_max.append(summary['score'])
    x0_max.append(summary['x0'])
    values_max.append(summary['values'])
    
print(f_max)
print(x0_max)
print(values_max)


plt.figure()
plt.semilogy(np.array(f_max)[:,:,0].flatten())
plt.semilogy(np.array(f_max)[:,:,1].flatten())
plt.semilogy(np.array(f_max)[:,:,2].flatten())
plt.semilogy(np.array(f_max)[:,:,3].flatten())


plt.figure()
plt.semilogy(np.array(f_max)[:,:,:].sum(-1).flatten())

x0_max=np.array(x0_max)

values_max=np.array(values_max).flatten()

qubit_anharmonicity = [values['qubit_anharmonicity'] for values in values_max]
cav_DS = [values['cav_DS'] for values in values_max]
Freq_cav = [values['Freq_cav'] for values in values_max]
cav_Q = [values['cav_Q'] for values in values_max]
Freq_qubit = [values['Freq_qubit'] for values in values_max]



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



