Examples
=========

Start-up quick example
-----------------------

The following code illustrates how to perform a complete analysis of a
simple two-qubit, one-cavity device in just a few lines of code with
``pyEPR``. In the HFSS file, before running the script, first specify
the non-linear junction rectangles and variables (see Sec. pyEPR Project
Setup in HFSS). All operations in the eigen analysis and Hamiltonian
computation are fully automated. The results are saved, printed, and
succinctly plotted.

.. code:: python

    # Load pyEPR. See the tutorial notebooks!
    import pyEPR as epr

    # 1. Connect to your Ansys, and load your design
    pinfo = epr.ProjectInfo(project_path = r'C:\sim_folder',
                            project_name = r'cavity_with_two_qubits',
                            design_name  = r'Alice_Bob')


    # 2a. Non-linear (Josephson) junctions
    pinfo.junctions['jAlice'] = {'Lj_variable':'Lj_alice', 'rect':'rect_alice', 'line': 'line_alice', 'Cj_variable':'Cj_alice'}
    pinfo.junctions['jBob']   = {'Lj_variable':'Lj_bob',   'rect':'rect_bob',   'line': 'line_bob', 'Cj_variable':'Cj_bob'}
    pinfo.validate_junction_info() # Check that valid names of variables and objects have been supplied.

    # 2b. Dissipative elements: specify
    pinfo.dissipative['dielectrics_bulk']    = ['si_substrate', 'dielectric_object2'] # supply names of hfss objects
    pinfo.dissipative['dielectric_surfaces'] = ['interface1', 'interface2']
    # Alternatively, these could be specified in ProjectInfo with
    # pinfo = epr.ProjectInfo(..., dielectrics_bulk = ['si_substrate', 'dielectric_object2'])

    # 3.  Perform microwave analysis on eigenmode solutions
    eprd = epr.DistributedAnalysis(pinfo)
    swp_var = 'Lj_alice' # Sweep variable from optimetric analysis that should be used on the x axis for the frequency plot
    eprd.quick_plot_frequencies(swp_var) # plot the solved frequencies before the analysis
    eprd.hfss_report_full_convergence() # report convergence
    eprd.do_EPR_analysis()

    # 4a.  Perform Hamiltonian spectrum post-analysis, building on mw solutions using EPR
    epra = epr.QuantumAnalysis(eprd.data_filename)
    epra.analyze_all_variations(cos_trunc = 8, fock_trunc = 7)

    # 4b. Report solved results
    swp_variable = 'Lj_alice' # suppose we swept an optimetric analysis vs. inductance Lj_alice
    epra.plot_hamiltonian_results(swp_variable=swp_variable)
    epra.report_results(swp_variable=swp_variable, numeric=True)
    epra.quick_plot_mode(0,0,1,numeric=True, swp_variable=swp_variable)

Video tutorials
------------------------------------------

.. raw:: html

    <div style="overflow:auto;">
    <table style="">
    <tr>
        <th>
        <a href="https://www.youtube.com/watch?v=fSRYvD-ITnQ&list=PLnak_fVcHp17tydgFosNtetDNjKbEaXtv&index=1">
        Tutorial 1 - Overview
        <br>
        <img src="https://img.youtube.com/vi/fSRYvD-ITnQ/0.jpg" alt="pyEPR Tutorial 1 - Overview" width=250>
        </a>
        </th>
        <th>
        <a href="https://www.youtube.com/watch?v=ZTi1pb6wSbE&list=PLnak_fVcHp17tydgFosNtetDNjKbEaXtv&index=2">
        Tutorial 2 - Setup of Conda & Git
        <br>
        <img src="https://img.youtube.com/vi/ZTi1pb6wSbE/0.jpg" alt="pyEPR Tutorial 2 - Setup of Conda & Git" width=250>
        </a>
        </th>
        <th>
        <a href="https://www.youtube.com/watch?v=L79nlXY2w4s&list=PLnak_fVcHp17tydgFosNtetDNjKbEaXtv&index=3">
        Tutorial 3 - Setup of Packages & Config
        <br>
        <img src="https://img.youtube.com/vi/L79nlXY2w4s/0.jpg" alt="pyEPR Tutorial 3 - Setup of Packages & Config" width=250>
        </a>
        </th>
    </tr>
    </table>
    </div>



Jupyter notebooks and example Ansys files
------------------------------------------

The most extensive way to learn pyEPR is to look over the pyEPR `example notebooks`_ and `example files`_ contained in the package repo.

The Jupyter notebooks can be viewed in a web browser `directly here`_.

.. _example notebooks: https://github.com/zlatko-minev/pyEPR/tree/master/_tutorial_notebooks
.. _example files: https://github.com/zlatko-minev/pyEPR/tree/master/_example_files
.. _directly here: https://nbviewer.jupyter.org/github/zlatko-minev/pyEPR/tree/master/_tutorial_notebooks/