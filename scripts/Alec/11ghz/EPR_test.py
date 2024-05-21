# -*- coding: utf-8 -*-
"""
Created on Wed Aug 23 10:34:00 2017

@author: alec-eickbusch
"""

from pyEPR import *

if 0:
    # Specify the HFSS project to be analyzed
    project_info = ProjectInfo(
        r"C:\Users\awe4\Documents\Backed\hfss_simulations\11ghz\\"
    )
    project_info.project_name = "11ghz_alec"  # Name of the project file (string). "None" will get the current active one.
    project_info.design_name = "11ghz_design1"  # Name of the design file (string). "None" will get the current active one.
    project_info.setup_name = (
        None  # Name of the setup(string). "None" will get the current active one.
    )

    project_info.junctions["bot_junc"] = {
        "rect": "bot_junction",
        "line": "bot_junc_line",
        "Lj_variable": "bot_lj",
        "length": 0.0001,
    }
    project_info.junctions["top_junc"] = {
        "rect": "top_junction",
        "line": "top_junc_line",
        "Lj_variable": "top_lj",
        "length": 0.0001,
    }

    project_info.dissipative[
        "dielectric_surfaces"
    ] = None  # supply names here, there are more options in  project_info.dissipative.

    epr_hfss = DistributedAnalysis(project_info)
    epr_hfss.do_EPR_analysis()

#%%
if 1:
    epr = QuantumAnalysis(epr_hfss.data_filename)  # Analysis results
    # result = epr.analyze_variation('1', cos_trunc = 8, fock_trunc = 7)
    epr.analyze_all_variations(cos_trunc=10, fock_trunc=7)

#%%
if 1:
    PM = OrderedDict()
    for n in range(3, 10):
        epr.analyze_all_variations(cos_trunc=10, fock_trunc=n)
        PM[n] = epr.results["0"]["chi_ND"]

    {k: v[0][0] for k, v in PM.items()}
