# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 11:21:01 2017

@author: Zlatko
"""

from pyEPR import *

if 1:
    # Specify the HFSS project to be analyzed
    project_info = ProjectInfo(
        r"C:\\Users\\rslqulab\Desktop\\Lysander\participation_ratio_project\\Shyam's autonomous stabilization simulations\\"
    )
    project_info.project_name = "2017_08_Zlatko_Shyam_AutStab"  # Name of the project file (string). "None" will get the current active one.
    project_info.design_name = "2 pyEPR"  # Name of the design file (string). "None" will get the current active one.
    project_info.setup_name = (
        None  # Name of the setup(string). "None" will get the current active one.
    )

    ## Describe the junctions in the HFSS design
    project_info.junctions["jAlice"] = {
        "rect": "qubitAlice",
        "line": "alice_line",
        "Lj_variable": "LJAlice",
        "length": 0.0001,
    }
    project_info.junctions["jBob"] = {
        "rect": "qubitBob",
        "line": "bob_line",
        "Lj_variable": "LJBob",
        "length": 0.0001,
    }

    # Dissipative elements EPR
    project_info.dissipative[
        "dielectric_surfaces"
    ] = None  # supply names here, there are more options in  project_info.dissipative.

    # Run analysis
    epr_hfss = DistributedAnalysis(project_info)
    epr_hfss.do_EPR_analysis()

if 1:  # Hamiltonian analysis
    filename = epr_hfss.data_filename
    # filename = r'C:\\Users\\rslqulab\\Desktop\\zkm\\2017_pyEPR_data\\\\/2017_08_Zlatko_Shyam_AutStab/2 pyEPR/2 pyEPR_20170825_170550.hdf5'
    epr = QuantumAnalysis(filename)

    # result = epr.analyze_variation('1', cos_trunc = 8, fock_trunc = 7)
    epr.analyze_all_variations(cos_trunc=8, fock_trunc=7)
    epr.plot_hamiltonian_results()
