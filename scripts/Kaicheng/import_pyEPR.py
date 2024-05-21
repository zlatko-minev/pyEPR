# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 11:21:01 2017

@author: Zlatko
"""

from pyEPR import *

if 0:
    # Specify the HFSS project to be analyzed
    project_info = ProjectInfo(r"X:\Simulation\\hfss\\KC\\")
    project_info.project_name = "2013-12-03_9GHzCavity"  # Name of the project file (string). "None" will get the current active one.
    project_info.design_name = "9GHz_EM_center_SNAIL"  # Name of the design file (string). "None" will get the current active one.
    project_info.setup_name = (
        None  # Name of the setup(string). "None" will get the current active one.
    )

    ## Describe the junctions in the HFSS design
    project_info.junctions["snail"] = {
        "rect": "qubit",
        "line": "JunctionLine",
        "Lj_variable": "LJ",
        "length": 0.0001,
    }
    #    project_info.junctions['jBob']   = {'rect':'qubitBob',    'line': 'bob_line',   'Lj_variable':'LJBob',   'length':0.0001}

    # Dissipative elements EPR
    project_info.dissipative[
        "dielectric_surfaces"
    ] = None  # supply names here, there are more options in  project_info.dissipative.

    # Run analysis
    epr_hfss = DistributedAnalysis(project_info)
    epr_hfss.do_EPR_analysis()  # variations = ['1', '70']

if 1:  # Hamiltonian analysis
    #    filename = epr_hfss.data_filename
    filename = r"X:\Simulation\hfss\KC\pyEPR_results_2018\2013-12-03_9GHzCavity\9GHz_EM_center_SNAIL\9GHz_EM_center_SNAIL_20180726_170049.hdf5"
    # filename = r'C:\\Users\\rslqulab\\Desktop\\zkm\\2017_pyEPR_data\\\\/2017_08_Zlatko_Shyam_AutStab/2 pyEPR/2 pyEPR_20170825_170550.hdf5'
    epr = QuantumAnalysis(filename)

    # result = epr.analyze_variation('1', cos_trunc = 8, fock_trunc = 7)
    epr.analyze_all_variations(cos_trunc=None, fock_trunc=4)  # only quadratic part
    epr.plot_hamiltonian_results()

    if 1:
        f0 = epr.results.get_frequencies_HFSS()
        f1 = epr.results.get_frequencies_O1()
        chi = epr.results.get_chi_O1()
        mode_idx = list(f0.index)
        nmodes = len(mode_idx)
        cmap = cmap_discrete(nmodes)
