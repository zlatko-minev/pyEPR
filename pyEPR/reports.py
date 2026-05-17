"""
Module for reporting utility functions

@author: Zlatko K Minev
"""

import pandas as pd

from .toolbox.plotting import legend_translucent, plt


def _style_plot_convergence(
    ax, ylabel=None, xlabel="Pass number", ylabel_col="k", y_title=False
):
    ax.set_xlabel(xlabel)
    if ylabel:
        if y_title:
            ax.set_title(ylabel)
        else:
            ax.set_ylabel(ylabel, color=ylabel_col)
    ax.grid()
    ax.autoscale(tight=False)
    ax.set_axisbelow(True)  # Don't allow the axis to be on top of your data
    ax.minorticks_on()
    ax.grid(which="minor", linestyle=":", linewidth="0.5", color="black", alpha=0.2)
    ax.grid(which="major", alpha=0.5)


_style_plot_conv_kw = dict(marker="o", ms=4)


def plot_convergence_max_df(ax, s, kw={}, color="r"):
    """For a single pass"""
    s.plot(ax=ax, **{**dict(c="r"), **_style_plot_conv_kw, **kw})
    if (s > 0).any():
        ax.set_yscale("log")
    _style_plot_convergence(ax)
    fig = ax.figure
    fig.text(0.45, 0.95, s.name, ha="center", va="bottom", size="medium", color=color)
    ax.tick_params(axis="y", labelcolor=color)
    # ax.axhline(1.0, color='k', lw=1.5,alpha= 0.35)
    # ax.axhline(0.1, color='k', lw=1.5,alpha= 0.35)
    ax.grid(which="minor", linestyle=":", linewidth="0.5", color=color, alpha=0.25)
    ax.grid(which="major", color="#c4abab", alpha=0.5)
    ax.spines["left"].set_color(color)


def plot_convergence_solved_elem(ax, s, kw={}, color="b"):
    """For a single pass"""
    (s / 1000).plot(ax=ax, **{**dict(c="b"), **_style_plot_conv_kw, **kw})
    _style_plot_convergence(ax)
    # ax.set_ylim([100,None])
    # ax.set_yscale("log")
    ax.minorticks_off()
    ax.grid(False)
    ax.tick_params(axis="y", labelcolor=color)
    # ax.ticklabel_format(style='sci',scilimits=(0,0))
    fig = ax.figure
    fig.text(
        0.6,
        0.95,
        "Solved elements (1000s)",
        ha="center",
        va="bottom",
        size="medium",
        color=color,
    )
    ax.spines["left"].set_color("r")
    ax.spines["right"].set_color(color)


def plot_convergence_f_vspass(ax, s, kw={}):
    """For a single pass"""
    if s is not None:
        (s).plot(ax=ax, **{**_style_plot_conv_kw, **kw})
        _style_plot_convergence(ax, "Eigenmode f vs. pass [GHz]", y_title=True)
        legend_translucent(ax, leg_kw=dict(fontsize=6))


def plot_convergence_maxdf_vs_sol(ax, s, s2, kw={}):
    """
    ax, 'Max Δf %', 'Solved elements', kw for plot
    """
    s = s.copy()
    s.index = s2
    (s).plot(ax=ax, **{**_style_plot_conv_kw, **kw})
    _style_plot_convergence(ax, s.name, xlabel="Solved elements", y_title=True)
    if (s > 0).any():
        ax.set_yscale("log")
    if (s2 > 0).any():
        ax.set_xscale("log")


def plot_q3d_convergence_main(epr, RES):
    """Plot EPR convergence: loss (alpha) and frequency vs. pass number.

    Parameters
    ----------
    epr : QuantumAnalysis
        A pyEPR QuantumAnalysis object with HFSS convergence data available.
    RES : dict
        Results dictionary from ``QuantumAnalysis.analyze_all_variations``,
        containing ``"alpha"`` (loss rate, MHz) and ``"fQ"`` (frequency, GHz)
        as pandas Series indexed by pass number.

    Returns
    -------
    matplotlib.figure.Figure
        Figure with alpha (blue, left axis) and frequency in MHz (red, right axis)
        plotted vs. simulation pass number.
    """
    fig = epr.hfss_report_full_convergence(_display=False)

    ax = fig.axes[0]
    ax2 = ax.twinx()
    ax.cla()
    ax2.cla()
    RES["alpha"].plot(ax=ax, c="b")
    (RES["fQ"] * 1000).plot(ax=ax2, c="red")
    from matplotlib import pyplot as plt

    _style_plot_convergence(ax, "Alpha (blue),  Freq (red) [MHz]", y_title=True)
    ax2.set_ylabel("Frequency (MHz)", color="r")
    ax.set_ylabel("Alpha(MHz)", color="b")
    ax2.spines["right"].set_color("r")
    ax2.tick_params(axis="y", labelcolor="r")
    ax.tick_params(axis="y", labelcolor="b")
    ax.set_xlabel("Pass")
    fig.tight_layout()

    return fig


def plot_q3d_convergence_chi_f(RES):
    """Plot EPR convergence: χ (dispersive shift) and g (coupling) vs. pass number.

    Parameters
    ----------
    RES : dict
        Results dictionary from ``QuantumAnalysis.analyze_all_variations``,
        containing ``"chi_in_MHz"`` and ``"gbus"`` as columns of numeric values
        indexed by pass number.

    Returns
    -------
    matplotlib.figure.Figure
        Figure with two subplots: χ convergence (left) and g convergence (right),
        both in MHz.
    """
    df_chi = pd.DataFrame(RES["chi_in_MHz"].values.tolist())
    df_chi.index.name = "Pass"
    df_g = pd.DataFrame(RES["gbus"].values.tolist())
    df_g.index.name = "Pass"

    fig, axs = plt.subplots(1, 2, figsize=(9, 3.5))
    df_chi.plot(lw=2, ax=axs[0])
    df_g.plot(lw=2, ax=axs[1])
    _style_plot_convergence(axs[0])
    _style_plot_convergence(axs[1])
    axs[0].set_title(r"$\chi$ convergence (MHz)")
    axs[1].set_title(r"$g$ convergence (MHz)")

    return fig


# Backward-compatible private aliases (used by qiskit-metal and older scripts)
_plot_q3d_convergence_main   = plot_q3d_convergence_main
_plot_q3d_convergence_chi_f  = plot_q3d_convergence_chi_f
