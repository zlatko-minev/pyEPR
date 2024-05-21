# Zlatko

from pyEPR import *
import matplotlib.pyplot as plt

if 1:
    # Specify the HFSS project to be analyzed
    project_info = ProjectInfo(r"C:\Users\rslqulab\Desktop\zkm\2017_pyEPR_data\\")
    project_info.project_name = "2017-10 re-sim SM22-R3C1"
    project_info.design_name = "3. sweep both"
    project_info.setup_name = None

    ## Describe the junctions in the HFSS design
    project_info.junctions["jBright"] = {
        "rect": "juncV",
        "line": "juncH_line",
        "Lj_variable": "LJ1",
        "length": 0.0001,
    }
    project_info.junctions["jDark"] = {
        "rect": "juncH",
        "line": "juncV_line",
        "Lj_variable": "LJ2",
        "length": 0.0001,
    }

    # Dissipative elements EPR
    project_info.dissipative[
        "dielectric_surfaces"
    ] = None  # supply names here, there are more options in  project_info.dissipative.

    # Run analysis
    epr_hfss = DistributedAnalysis(project_info)
    epr_hfss.do_EPR_analysis()

if 1:  # Analysis result
    filename = epr_hfss.data_filename
    # filename = r'C:\Users\rslqulab\Desktop\zkm\2017_pyEPR_data\\/2017-10 re-sim SM22-R3C1/1. R3C1/1. R3C1_20171016_110756.hdf5'
    epr = QuantumAnalysis(filename)
    epr.plot_convergence_f_lin()

    epr._renorm_pj = True

    plt.close("all")
    epr.analyze_all_variations(cos_trunc=10, fock_trunc=8)
    epr.plot_hamiltonian_results()
    print(epr.data_filename)


#%%
if 1:
    import numpy as np
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt

    sns.reset_orig()
    # epr.hfss_variables.loc['_LJ2']

    kw_map = dict(
        vmin=-20, vmax=20, linewidths=0.5, annot=True, cmap="seismic"
    )  # RdYlGn_r

    target_f = pd.Series([4688, 5300, 9003], index=["D", "B", "C"])
    target_alpha = pd.Series([148, 174], index=["D", "B"])
    target_chi = pd.Series([85, 5, 0.33], index=["DB", "BC", "DC"])

    results = epr.results
    f_ND = results.get_frequencies_ND().rename(index={0: "D", 1: "B", 2: "C"})
    f_error = f_ND.apply(lambda x: 100 * (x.values - target_f) / x, axis="index")

    fig, axs = plt.subplots(1, 3, figsize=(15, 7.5))
    sns.heatmap(f_error.transpose(), ax=axs[0], **kw_map)

    chis = results.get_chi_ND()
    chis = xarray_unravel_levels(chis, ["variation", "m", "n"])
    alpha_ND = sort_df_col(chis.sel_points(m=[0, 1], n=[0, 1]).to_pandas())
    alpha_ND.index = target_alpha.index
    alpha_ND_err = alpha_ND.apply(
        lambda x: 100 * (x.values - target_alpha) / x, axis="index"
    )
    sns.heatmap(alpha_ND_err.transpose(), ax=axs[1], **kw_map)

    chi_ND = sort_df_col(chis.sel_points(m=[0, 1, 0], n=[1, 2, 2]).to_pandas())
    chi_ND.index = target_chi.index
    chi_ND_err = chi_ND.apply(lambda x: 100 * (x.values - target_chi) / x, axis="index")
    sns.heatmap(chi_ND_err.transpose(), ax=axs[2], **kw_map)
    axs[0].set_title("Freq.")
    axs[1].set_title("Anharmonicities")
    axs[2].set_title("cross-Kerrs")
