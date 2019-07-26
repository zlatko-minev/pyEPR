from .toolbox_plotting import  legend_translucent, mpl, plt

def _style_plot_convergence(ax, ylabel=None, xlabel = 'Pass number', ylabel_col='k', y_title=False):
    ax.set_xlabel(xlabel)
    if ylabel:
        if y_title:
            ax.set_title(ylabel)
        else:
            ax.set_ylabel(ylabel,color=ylabel_col)
    ax.grid()
    ax.autoscale(tight=False)
    ax.set_axisbelow(True) # Don't allow the axis to be on top of your data
    ax.minorticks_on()
    ax.grid(which='minor', linestyle=':', linewidth='0.5', color='black', alpha=0.2)
    ax.grid(which='major', alpha=0.5)

    
_style_plot_conv_kw = dict(marker='o',ms=4)

def plot_convergence_max_df(ax, s, kw={}, color = 'r'):
    '''For a single pass'''
    s.plot(ax=ax, **{**dict(c='r'), **_style_plot_conv_kw, **kw})
    ax.set_yscale("log")
    _style_plot_convergence(ax)
    fig = ax.figure
    fig.text(0.45, 0.95, "Max Δf %", ha="center", va="bottom", size="medium",color=color)
    ax.tick_params(axis='y', labelcolor=color)
    #ax.axhline(1.0, color='k', lw=1.5,alpha= 0.35)
    #ax.axhline(0.1, color='k', lw=1.5,alpha= 0.35)
    ax.grid(which='minor', linestyle=':', linewidth='0.5', color=color, alpha=0.25)
    ax.grid(which='major', color='#c4abab', alpha=0.5)
    ax.spines['left'].set_color(color)

    
def plot_convergence_solved_elem(ax, s, kw={}, color = 'b'):
    '''For a single pass'''
    (s/1000).plot(ax=ax, **{**dict(c='b'), **_style_plot_conv_kw, **kw})
    _style_plot_convergence(ax )
    ax.set_ylim([100,None])
    #ax.set_yscale("log")
    ax.minorticks_off()
    ax.grid(False)
    ax.tick_params(axis='y', labelcolor=color)
    #ax.ticklabel_format(style='sci',scilimits=(0,0))
    fig = ax.figure
    fig.text(0.6, 0.95, 'Solved elements (1000s)', ha="center", va="bottom", size="medium",color=color)
    ax.spines['left'].set_color('r')
    ax.spines['right'].set_color(color)
    
    
def plot_convergence_f_vspass(ax, s, kw={}):
    '''For a single pass'''
    (s).plot(ax=ax, **{**_style_plot_conv_kw, **kw})
    _style_plot_convergence(ax, 'Eigenmode f vs. pass [GHz]', y_title=True)
    legend_translucent(ax, leg_kw =dict(fontsize=6))
    
def plot_convergence_maxdf_vs_sol(ax, s, s2, kw={}):
    '''
    ax, 'Max Δf %', 'Solved elements (1000s)', kw for plot
    '''
    s = s.copy()
    s.index = s2/1000
    (s).plot(ax=ax, **{**_style_plot_conv_kw, **kw})
    _style_plot_convergence(ax, 'Max Δf %', xlabel='Solved elements (1000s)', y_title=True)
    ax.set_yscale("log")
    #ax.set_xscale("log")