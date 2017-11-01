# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 11:21:01 2017

@author: Zlatko
"""

from pyEPR import *

if 1:
    # Specify the HFSS project to be analyzed
    project_info = Project_Info(r"D:\\hfss\\Catch_tomo\\tunage Q\\")
    project_info.project_name  = 'catch_tomo'  # Name of the project file (string). "None" will get the current active one.
    project_info.design_name   = 'First_Design'  # Name of the desgin file (string). "None" will get the current active one.
    project_info.setup_name    = None       # Name of the setup(string). "None" will get the current active one.

    ## Describe the junctions in the HFSS desgin
    project_info.junctions['Js'] = {'rect':'memory_jsjunc_junc',  'line': 'jsLine', 'Lj_variable':'$Lj', 'length': 0.00002}
    #project_info.junctions['jBob']   = {'rect':'qubitBob',    'line': 'bob_line',   'Lj_variable':'LJBob',   'length':0.0001}

    # Dissipative elments EPR
    project_info.dissipative.dielectric_surfaces = None         # supply names here, there are more options in  project_info.dissipative.

    # Run analysis
    epr_hfss    = pyEPR_HFSS(project_info)
    epr_hfss.do_EPR_analysis(modes=[1])

if 1: # Hamiltonian analysis
    filename = epr_hfss.data_filename
    #filename = r'C:\\Users\\rslqulab\\Desktop\\zkm\\2017_pyEPR_data\\\\/2017_08_Zlatko_Shyam_AutStab/2 pyEPR/2 pyEPR_20170825_170550.hdf5'
    epr      = pyEPR_Analysis(filename)

    result = epr.analyze_variation('0')
    #epr.analyze_all_variations(cos_trunc = 8, fock_trunc = 7)
    #epr.plot_Hresults()
