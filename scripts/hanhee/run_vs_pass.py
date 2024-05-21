# -*- coding: utf-8 -*-
"""
Example startup script to perform full quantization of a two qubit, one cavity Josephson circuit.
The results are saved, printed, and nicely plotted.

@author: Zlatko
"""

from pyEPR import ProjectInfo, DistributedAnalysis, QuantumAnalysis

# 1.  Project and design. Open link to HFSS controls.
project_info = ProjectInfo(
    "D:\LOGIQ-IBMQ\Cranes\HFSS simulation\\",
    project_name="2018-12-03 Zlatko pyEPR",  # Project file name (string). "None" will get the current active one.
    design_name="L-4 bus-EPR",  # Design name (string). "None" will get the current active one.
)

# 2a. Junctions. Specify junctions in HFSS model
for i in range(1, 3 + 1):  # specify N number of junctions
    i = str(i)
    project_info.junctions["j" + i] = {
        "Lj_variable": "Lj" + i,
        "rect": "Qubit" + i,
        "line": "Polyline" + i,
        "length": 30 * 10**-6,
    }

# 2b. Dissipative elements.
# project_info.dissipative['dielectrics_bulk']    = ['subs_Q1']    # supply names here, there are more options in project_info.dissipative.
# project_info.dissipative['dielectric_surfaces'] = ['interface']

#%%
# 3.  Run analysis
if 1:
    passes = range(1, 20, 1)
    epr_hfss = DistributedAnalysis(project_info)

    # CLEAR DATA
    # if not 'RES' in locals():
    from collections import OrderedDict

    RES = OrderedDict()

#%%%
setup_name = None
design = epr_hfss.design
setup_name = setup_name if not (setup_name is None) else design.get_setup_names()[0]
setup = design.get_setup(setup_name)
print(" HFSS setup name: %s" % setup_name)

#%%
from numpy import diag, sqrt, array
import pandas as pd
from pyEPR.toolbox import get_above_diagonal


def do_analysis(pass_, variation="0"):
    epr_hfss.do_EPR_analysis(variations=[variation])
    RES[pass_] = OrderedDict()

    epr = QuantumAnalysis(epr_hfss.data_filename)
    RES[pass_]["epr"] = epr
    # epr = RES[pass_]['epr']
    RES[pass_]["freq_hfss"] = epr.freqs_hfss[variation]
    RES[pass_]["Pmj_raw"] = epr.PM[variation]

    ## NORMED
    dum = epr.get_Pmj("0")
    RES[pass_]["Pmj_normed"] = dum["PJ"]  # DataFrame
    RES[pass_]["_Pmj_norm"] = dum["Pm_norm"]  # Series

    RES[pass_]["mats"] = epr.get_matrices(variation, print_=False)  # arrays
    # PJ, SJ, Om, EJ, PHI_zpf
    RES[pass_]["Hres"] = Hres = zkm_get_Hparams(*RES[pass_]["mats"])
    RES[pass_]["alpha"] = Hres["alpha"]  # easy access
    RES[pass_]["chi"] = Hres["chi"]

    ### Numerical diagonalization
    # RES[pass_]['ND']         = None
    if 1:
        print(" ND pass=%s variation=%s" % (pass_, variation))
        from pyEPR.core import pyEPR_ND

        f1_ND, CHI_ND = pyEPR_ND(
            epr.freqs_hfss[variation],
            epr.Ljs[variation],
            RES[pass_]["mats"][-1],  # PHI_zpf
            cos_trunc=10,
            fock_trunc=9,
        )
        RES[pass_]["ND"] = {}
        RES[pass_]["ND"]["f01"] = f1_ND
        RES[pass_]["ND"]["CHI"] = CHI_ND

    return Hres


def zkm_get_Hparams(PJ, SJ, Om, EJ, PHI):
    """
    Report all in MHz
    """
    M, J = PJ.shape
    res = {
        "alpha_p4R2_p6R1": [],  # {p=4,RWA=2} + {p=6,RWA2}
        "chi": {},  # {p=4,RWA=1}
        "omega_zx": {},
    }

    res["alpha_p4R2_p6R1"] = [
        1000
        * sum(
            [
                PJ[m, j] ** 2
                * (Om[m, m] ** 2)
                / (8.0 * EJ[j, j])
                * (
                    1.0
                    + PJ[m, j] * Om[m, m] / EJ[j, j] * (17.0 / 32.0 * PJ[m, j] - 0.25)
                )
                for j in range(J)
            ]
        )  # MHz
        for m in range(M)
    ]

    res["alpha_zpf"] = [
        1000
        * sum(
            [
                0.5 * EJ[j, j] * PHI[m, j] ** 4
                + 2.0 * 306.0 / Om[m, m] * (EJ[j, j] / 24.0) ** 2 * PHI[m, j] ** 8
                - 0.25 * EJ[j, j] * PHI[m, j] ** 6
                for j in range(J)
            ]
        )  # MHz
        for m in range(M)
    ]  # fully equivalent, checked

    ### Cross-Kerr
    for m in range(M):
        for m1 in range(m):
            res["chi"]["%d,%d" % (m1, m)] = 1000.0 * sum(
                [EJ[j, j] * PHI[m1, j] ** 2 * PHI[m, j] ** 2 for j in range(J)]
            )

    res["chi_zpf"] = {}  # debug only
    for m in range(M):
        for m1 in range(m):
            res["chi_zpf"]["%d,%d" % (m1, m)] = 1000.0 * sum(
                [
                    Om[m1, m1] * Om[m, m] * PJ[m1, j] * PJ[m, j] / (4.0 * EJ[j, j])
                    for j in range(J)
                ]
            )  # fully equivalent, checked

    res["chi2"] = {}  # higher order correction
    for m in range(M):
        for m1 in range(m):
            res["chi2"]["%d,%d" % (m1, m)] = 1000.0 * sum(
                [  # TODO: CHECK AND SIGN
                    (EJ[j, j] / 24.0) ** 2
                    * (
                        864.0
                        * PHI[m1, j] ** 6
                        * PHI[m, j] ** 2
                        / (+Om[m1, m1] - Om[m, m])
                        - 864.0
                        * PHI[m1, j] ** 2
                        * PHI[m, j] ** 6
                        / (-Om[m1, m1] + Om[m, m])
                        - 576.0 * PHI[m1, j] ** 6 * PHI[m, j] ** 2 / (Om[m1, m1])
                        - 576.0 * PHI[m1, j] ** 2 * PHI[m, j] ** 6 / (Om[m, m])
                        - 576.0
                        * PHI[m1, j] ** 4
                        * PHI[m, j] ** 4
                        / (Om[m1, m1] + Om[m, m])
                        + 288.0
                        * PHI[m1, j] ** 2
                        * PHI[m, j] ** 6
                        / (Om[m1, m1] - 3 * Om[m, m])
                        + 288.0
                        * PHI[m1, j] ** 6
                        * PHI[m, j] ** 2
                        / (Om[m, m] - 3 * Om[m1, m1])
                    )
                    for j in range(J)
                ]
            )

    ### CR Gate - analytical {p=4, RWA1}
    res["omega_zx"] = {}
    for m in range(M):
        for m1 in range(m):
            if m != m1:
                res["omega_zx"]["%d,%d" % (m, m1)] = 1000.0 * abs(
                    sum([EJ[j, j] * PHI[m, j] ** 3 * PHI[m1, j] for j in range(J)])
                )  # fully equivalent, checked

    res["alpha"] = res["alpha_p4R2_p6R1"]  # MHz
    res["f01"] = (
        1000 * diag(Om)
        - res["alpha"]
        - [
            sum(
                [
                    res["chi"]["%d,%d" % (min(m1, m), max(m1, m))]
                    for m1 in range(M)
                    if m1 != m
                ]
            )
            for m in range(M)
        ]
    )
    # MHz -kerrs .. todo # f_01 frequency
    return res


if 0:  # update all passes
    for pass_ in RES.keys():
        do_analysis(pass_)

# res = do_analysis(pass_)
# pd.Series(res['alpha'])


#%%
def do_plot(RES):
    """
    Make sure
        %matplotlib qt
    TODO: in future just setup once, and then update lines only
    """
    # live plot https://stackoverflow.com/questions/11874767/how-do-i-plot-in-real-time-in-a-while-loop-using-matplotlib
    # see also pylive
    import matplotlib

    matplotlib.use("Qt5Agg")
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np
    from pyEPR.toolbox_plotting import legend_translucent
    from pyEPR.toolbox import combinekw, xarray_unravel_levels, floor_10

    plt.ion()

    fig = plt.figure(1, figsize=(25, 10))
    fig.clf()
    fig, axs = plt.subplots(2, 3, subplot_kw=dict(), num=1, sharex=True)

    kw = dict(marker="o")
    leg_kw = dict(fontsize=9, ncol=1)

    # Frequency
    ax = axs[0, 0]
    df = pd.DataFrame({x: 1000.0 * RES[x]["freq_hfss"] for x in RES}).transpose()
    df.plot(ax=ax, **kw)
    ax.set_title("Linear mode, HFSS frequency $\omega_m/2\pi$ and dressed (MHz)")
    legend_translucent(
        ax, leg_kw=combinekw(leg_kw, dict(title="Mode #"))
    )  # ax.legend(title= 'Mode #')
    # Dressed frequency
    ax.set_prop_cycle(None)
    df = pd.DataFrame({x: RES[x]["Hres"]["f01"] for x in RES}).transpose()
    df.plot(ax=ax, legend=False, **combinekw(kw, dict(marker=None, alpha=0.5, ls="--")))

    # Pmj Norm
    ax = axs[1, 0]
    df = pd.DataFrame({x: RES[x]["_Pmj_norm"] for x in RES}).transpose()
    df.plot(ax=ax, **kw)
    ax.set_title("HFSS $p_{mj}$ norm")
    legend_translucent(ax, leg_kw=combinekw(leg_kw, dict(title="Mode #")))

    # Frequency
    ax = axs[0, 1]
    df = pd.DataFrame({x: RES[x]["alpha"] for x in RES}).transpose()
    df.plot(ax=ax, **kw)
    ax.set_title(r"Anharmonicity $\alpha_{mj}/2\pi$ (MHz)")
    legend_translucent(ax, leg_kw=combinekw(leg_kw, dict(title="Mode #")))
    if RES[pass_]["ND"] is not None:  # plot numerical solution
        ax.set_prop_cycle(None)
        df = pd.DataFrame({x: diag(RES[x]["ND"]["CHI"]) for x in RES}).transpose()
        df.plot(
            ax=ax, legend=False, **combinekw(kw, dict(marker=None, alpha=0.5, ls="--"))
        )

    ax = axs[0, 2]
    df = pd.DataFrame({x: RES[x]["chi"] for x in RES}).transpose()
    df.plot(ax=ax, **kw)
    ax.set_title(r"Cross-Kerr $\chi_{mm\prime}/2\pi$ (MHz)")
    legend_translucent(ax, leg_kw=combinekw(leg_kw, dict(title="Mode #")))
    #    if RES[pass_]['ND'] is not None: # plot numerical solution
    #        ax.set_prop_cycle(None)
    #        df = pd.DataFrame({x:-get_above_diagonal(RES[x]['ND']['CHI']) for x in RES}).transpose()
    ##        df = pd.DataFrame({x:RES[x]['Hres']['chi2'] for x in RES}).transpose()
    #        df.plot(ax=ax,legend=False,**combinekw(kw,dict(marker=None,alpha=0.5,ls='--')))

    ax = axs[1, 2]
    df = pd.DataFrame({x: RES[x]["Hres"]["omega_zx"] for x in RES}).transpose()
    df.plot(ax=ax, **kw)
    ax.set_title(r"Cross-Resonance $\omega_{ZX}/2\pi$ (MHz)")
    legend_translucent(ax, leg_kw=combinekw(leg_kw, dict(title="C,T")))

    # Pmj normed plot
    ax = axs[1, 1]
    da = xarray_unravel_levels(
        {x: RES[x]["Pmj_normed"] for x in RES}, names=["pass", "mode", "junction"]
    )
    for mode in da.coords["mode"]:
        for junc in da.coords["junction"]:
            junc_name = str(junc.values)[2:]
            ys = da.sel(mode=mode, junction=junc)
            ys.plot.line(ax=ax, label="%2s,%4s" % (str(mode.values), junc_name), **kw)

    min_ = floor_10(min(abs(np.min(da.values)), abs(np.max(da.values))))  # just in case
    ax.set_ylim(min_, 1.05)
    ax.set_yscale("log", nonposy="clip")
    legend_translucent(ax, leg_kw=combinekw(leg_kw, dict(title="$p_{mj}$")))
    ax.set_title("HFSS $p_{mj}$ normed")

    from matplotlib.widgets import Button

    class Index(object):
        ind = 0

        def __init__(self, button, ax):
            self.ax = ax
            self.button = button  # so it doesnt get erased

        def next(self, event):
            i = self.ind = (self.ind + 1) % 5
            ax = self.ax
            if i == 0:
                ax.set_ylim(min_, 1.05)
                ax.set_yscale("log", nonposy="clip")
            elif i == 1:
                ax.set_ylim(min_, 1.02)
                ax.set_yscale("linear", nonposy="clip")
            elif i == 2:
                ax.set_ylim(0.8, 1.02)
                ax.set_yscale("linear", nonposy="clip")
            elif i == 3:
                ax.set_ylim(10**-3, 10**-1)
                ax.set_yscale("log", nonposy="clip")
            elif i == 4:
                ax.set_ylim(5 * 10**-5, 2 * 10**-3)
                ax.set_yscale("log", nonposy="clip")
            self.button.label.set_text("Next %d" % (self.ind))
            fig.canvas.draw()
            fig.canvas.flush_events()

    pos1 = ax.get_position()
    pos2 = [pos1.x0 + 0.0, pos1.y0 + pos1.height + 0.002, 0.07, 0.04]
    axnext = plt.axes(pos2)
    bnext = Button(axnext, "Next")
    callback = Index(bnext, ax)
    bnext.on_clicked(callback.next)

    for ax in np.ndarray.flatten(axs):
        ax.set_xlabel("Pass number")
        # ax.autoscale(tight=True)

    fig.tight_layout()  # pad=0.4, w_pad=0.5, h_pad=1.0)
    fig.show()
    plt.pause(0.01)
    return df


# do_plot(RES)

# import threading
# t = threading.Thread(target=do_plot, args = (RES,))
# t.start()

#%%
import time

if 1:

    for pass_ in passes:
        print(" Running pass #%s" % (pass_), end="")
        setup.passes = str(pass_)
        try:
            ret = (
                setup.solve()
            )  # I tried to use a worker thread but this gets complicated with COM blocking interface
            if ret in [0, "0"]:
                print(".  Normal completion.")
                time.sleep(0.5)
                do_analysis(pass_)
                do_plot(RES)
            elif ret in ["-1", -1]:
                print(".  Simulation error.")
            print(ret)
        except KeyboardInterrupt:
            print("\n\n Keyboard interruption...")
            break
        # ABORT: -2147352567, 'Exception occurred.', (0, None, None, None, 0, -2147024349), None)
    do_plot(RES)

#%%
if 0:
    epr_hfss.do_EPR_analysis()

    # 4.  Hamiltonian analysis
    epr = QuantumAnalysis(epr_hfss.data_filename)
    epr.analyze_all_variations(cos_trunc=8, fock_trunc=7)
    epr.plot_hamiltonian_results()
