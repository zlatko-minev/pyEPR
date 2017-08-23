# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 11:21:01 2017

@author: rslqulab

LOGGING:
---------------------
import logging

logging.basicConfig(format='%(asctime)s,%(msecs)d %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s',
    datefmt='%d-%m-%Y:%H:%M:%S',
    level=logging.DEBUG)

logger = logging.getLogger('stackoverflow_rocks')
logger.debug("This is a debug log")
logger.info("This is an info log")
logger.critical("This is critical")
logger.error("An error occurred")
"""

from pyEPR import *

if 1:
    # Specify the HFSS project to be analyzed
    project_info = Project_Info(r"C:\\Users\\rslqulab\Desktop\\Lysander\participation_ratio_project\\Shyam's autonomous stabilization simulations\\")
    project_info.project_name  = '2017_08_Zlatko_Shyam_AutStab'  # Name of the project file (string). "None" will get the current active one.
    project_info.design_name   = 'TEST'  # Name of the desgin file (string). "None" will get the current active one.
    project_info.setup_name    = None    # Name of the setup(string). "None" will get the current active one.

    ## Describe the junctions in the HFSS desgin
    project_info.junctions['jAlice'] = {'rect':'qubitAlice',  'line': 'alice_line', 'Lj_variable':'LJAlice', 'length':0.0001}
    project_info.junctions['jBob']   = {'rect':'qubitBob',    'line': 'bob_line',   'Lj_variable':'LJBob',   'length':0.0001}

    # Dissipative elments EPR
    project_info.dissipative.dielectric_surfaces = None         # supply names here, there are more options in  project_info.dissipative.

    # Run analysis
    epr_hfss    = pyEPR_HFSS(project_info)
    epr_hfss.do_EPR_analysis()

if 1:
    epr = pyEPR_Analysis(epr_hfss.data_filename)   # Analysis results
    epr.analyze_variation('0', cos_trunc = 8, fock_trunc = 4)