# -*- coding: utf-8 -*-
"""
Created on Wed Dec 25 20:35:42 2019

@author: Asaf
"""
from scipy.constants import hbar, e,pi,mega,nano
def Lj2Ej(Lj):
    """Lj in nH
    Ej in MHZ"""
    return 2*pi(((hbar/(2*e))**2)/(Lj/nano))/(hbar*mega)
    
    
