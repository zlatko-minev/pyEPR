'''
Numerical diagonalization of quantum Hamiltonian and parameter
extraction.

@author: Phil Reinhold, Zlatko Minev, Lysander Christakis

Original code on black_box_hamiltonian and make_dispersive functions by Phil Reinhold
Revisions and updates by Zlatko Minev & Lysander Christakis
'''
# pylint: disable=invalid-name


from __future__ import print_function

from functools import reduce

import numpy as np

from .constants import Planck as h
from .constants import fluxQ, hbar
from .hamiltonian import MatrixOps

try:
    import qutip
    from qutip import basis, tensor
except (ImportError, ModuleNotFoundError):
    pass

__all__ = [ 'epr_numerical_diagonalization',
            'make_dispersive',
            'black_box_hamiltonian',
            'black_box_hamiltonian_nq']

dot = MatrixOps.dot
cos_approx = MatrixOps.cos_approx


# ==============================================================================
# ANALYSIS FUNCTIONS
# ==============================================================================

def epr_numerical_diagonalization(freqs, Ljs, ϕzpf,
             cos_trunc=8,
             fock_trunc=9,
             use_1st_order=False,
             return_H=False,
             non_linear_potential=None):
    '''
    Numerical diagonalization for pyEPR. Ask Zlatko for details.

    :param fs: (GHz, not radians) Linearized model, H_lin, normal mode frequencies in Hz, length M
    :param ljs: (Henries) junction linearized inductances in Henries, length J
    :param fzpfs: (reduced) Reduced Zero-point fluctuation of the junction fluxes for each mode
                across each junction, shape MxJ

    :return: Hamiltonian mode freq and dispersive shifts. Shifts are in MHz.
             Shifts have flipped sign so that down shift is positive.
    '''

    freqs, Ljs, ϕzpf = map(np.array, (freqs, Ljs, ϕzpf))
    assert(all(freqs < 1E6)
           ), "Please input the frequencies in GHz. \N{nauseated face}"
    assert(all(Ljs < 1E-3)
           ), "Please input the inductances in Henries. \N{nauseated face}"

    Hs = black_box_hamiltonian(freqs * 1E9, Ljs.astype(np.float), fluxQ*ϕzpf,
                               cos_trunc, fock_trunc, individual=use_1st_order,
                               non_linear_potential = non_linear_potential)
    f_ND, χ_ND, _, _ = make_dispersive(
        Hs, fock_trunc, ϕzpf, freqs, use_1st_order=use_1st_order)
    χ_ND = -1*χ_ND * 1E-6  # convert to MHz, and flip sign so that down shift is positive

    return (f_ND, χ_ND, Hs) if return_H else (f_ND, χ_ND)




def black_box_hamiltonian(fs, ljs, fzpfs, cos_trunc=5, fock_trunc=8, individual=False,
                          non_linear_potential = None):
    r"""
    :param fs: Linearized model, H_lin, normal mode frequencies in Hz, length N
    :param ljs: junction linearized inductances in Henries, length M
    :param fzpfs: Zero-point fluctuation of the junction fluxes for each mode across each junction,
                 shape MxJ
    :return: Hamiltonian in units of Hz (i.e H / h)
    All in SI units. The ZPF fed in are the generalized, not reduced, flux.

    Description:
     Takes the linear mode frequencies, :math:`\omega_m`, and the zero-point fluctuations, ZPFs, and
     builds the Hamiltonian matrix of :math:`H_{full}`, assuming cos potential.
    """
    n_modes = len(fs)
    njuncs = len(ljs)
    fs, ljs, fzpfs = map(np.array, (fs, ljs, fzpfs))
    ejs = fluxQ**2 / ljs
    fjs = ejs / h

    fzpfs = np.transpose(fzpfs)  # Take from MxJ  to JxM

    assert np.isnan(fzpfs).any(
    ) == False, "Phi ZPF has NAN, this is NOT allowed! Fix me. \n%s" % fzpfs
    assert np.isnan(ljs).any(
    ) == False, "Ljs has NAN, this is NOT allowed! Fix me."
    assert np.isnan(
        fs).any() == False, "freqs has NAN, this is NOT allowed! Fix me."
    assert fzpfs.shape == (njuncs, n_modes), "incorrect shape for zpf array, {} not {}".format(
        fzpfs.shape, (njuncs, n_modes))
    assert fs.shape == (n_modes,), "incorrect number of mode frequencies"
    assert ejs.shape == (njuncs,), "incorrect number of qubit frequencies"

    def tensor_out(op, loc):
        "Make operator <op> tensored with identities at locations other than <loc>"
        op_list = [qutip.qeye(fock_trunc) for i in range(n_modes)]
        op_list[loc] = op
        return reduce(qutip.tensor, op_list)

    a = qutip.destroy(fock_trunc)
    ad = a.dag()
    n = qutip.num(fock_trunc)
    mode_fields = [tensor_out(a + ad, i) for i in range(n_modes)]
    mode_ns = [tensor_out(n, i) for i in range(n_modes)]

    def cos(x):
        return cos_approx(x, cos_trunc=cos_trunc)
    
    if non_linear_potential is None:
        non_linear_potential = cos

    linear_part = dot(fs, mode_ns)
    cos_interiors = [dot(fzpf_row/fluxQ, mode_fields) for fzpf_row in fzpfs]
    nonlinear_part = dot(-fjs, map(non_linear_potential, cos_interiors))
    if individual:
        return linear_part, nonlinear_part
    else:
        return linear_part + nonlinear_part

bbq_hmt = black_box_hamiltonian

def make_dispersive(H, fock_trunc, fzpfs=None, f0s=None, chi_prime=False,
                    use_1st_order=False):
    r"""
    Input: Hamiltonian Matrix.
        Optional: phi_zpfs and normal mode frequencies, f0s.
        use_1st_order : deprecated
    Output:
        Return dressed mode frequencies, chis, chi prime, phi_zpf flux (not reduced), and linear frequencies
    Description:
        Takes the Hamiltonian matrix `H` from bbq_hmt. It them finds the eigenvalues/eigenvectors and  assigns quantum numbers to them --- i.e., mode excitations,  such as, for instance, for three mode, :math:`|0,0,0\rangle` or :math:`|0,0,1\rangle`, which correspond to no excitations in any of the modes or one excitation in the 3rd mode, resp.    The assignment is performed based on the maximum overlap between the eigenvectors of H_full and H_lin.   If this crude explanation is confusing, let me know, I will write a more detailed one |:slightly_smiling_face:|
        Based on the assignment of the excitations, the function returns the dressed mode frequencies :math:`\omega_m^\prime`, and the cross-Kerr matrix (including anharmonicities) extracted from the numerical diagonalization, as well as from 1st order perturbation theory.
        Note, the diagonal of the CHI matrix is directly the anharmonicity term.
    """
    if hasattr(H, '__len__'):  # is it an array / list?
        [H_lin, H_nl] = H
        H = H_lin + H_nl
    else:  # make sure its a quanutm object
        assert type(
            H) == qutip.qobj.Qobj, "Please pass in either a list of Qobjs or Qobj for the Hamiltonian"

    print("Starting the diagonalization")
    evals, evecs = H.eigenstates()
    print("Finished the diagonalization")
    evals -= evals[0]

    N = int(np.log(H.shape[0]) / np.log(fock_trunc))    # number of modes
    assert H.shape[0] == fock_trunc ** N

    def fock_state_on(d):
        ''' d={mode number: # of photons} '''
        return qutip.tensor(*[qutip.basis(fock_trunc, d.get(i, 0)) for i in range(N)])  # give me the value d[i]  or 0 if d[i] does not exist

    if use_1st_order:
        num_modes = N
        print("Using 1st O")

        def multi_index_2_vector(d, num_modes, fock_trunc):
            return tensor([basis(fock_trunc, d.get(i, 0)) for i in range(num_modes)])
            '''this function creates a vector representation a given fock state given the data for excitations per
                        mode of the form d={mode number: # of photons}'''

        def find_multi_indices(fock_trunc):
            multi_indices = [{ind: item for ind, item in enumerate([i, j, k])} for i in range(fock_trunc)
                             for j in range(fock_trunc)
                             for k in range(fock_trunc)]
            return multi_indices
            '''this function generates all possible multi-indices for three modes for a given fock_trunc'''

        def get_expect_number(left, middle, right):
            return (left.dag()*middle*right).data.toarray()[0, 0]
            '''this function calculates the expectation value of an operator called "middle" '''

        def get_basis0(fock_trunc, num_modes):
            multi_indices = find_multi_indices(fock_trunc)
            basis0 = [multi_index_2_vector(
                multi_indices[i], num_modes, fock_trunc) for i in range(len(multi_indices))]
            evalues0 = [get_expect_number(v0, H_lin, v0).real for v0 in basis0]
            return multi_indices, basis0, evalues0
            '''this function creates a basis of fock states and their corresponding eigenvalues'''

        def closest_state_to(vector0):

            def PT_on_vector(original_vector, original_basis, pertub, energy0, evalue):
                new_vector = 0 * original_vector
                for i in range(len(original_basis)):
                    if (energy0[i]-evalue) > 1e-3:
                        new_vector += ((original_basis[i].dag()*H_nl*original_vector).data.toarray()[
                                       0, 0])*original_basis[i]/(evalue-energy0[i])
                    else:
                        pass
                return (new_vector + original_vector)/(new_vector + original_vector).norm()
                '''this function calculates the normalized vector with the first order correction term
                   from the non-linear hamiltonian '''

            [multi_indices, basis0, evalues0] = get_basis0(
                fock_trunc, num_modes)
            evalue0 = get_expect_number(vector0, H_lin, vector0)
            vector1 = PT_on_vector(vector0, basis0, H_nl, evalues0, evalue0)

            index = np.argmax([(vector1.dag() * evec).norm()
                               for evec in evecs])
            return evals[index], evecs[index]

    else:
        def closest_state_to(s):
            def distance(s2):
                return (s.dag() * s2[1]).norm()
            return max(zip(evals, evecs), key=distance)

    f1s = [closest_state_to(fock_state_on({i: 1}))[0] for i in range(N)]
    chis = [[0]*N for _ in range(N)]
    chips = [[0]*N for _ in range(N)]
    for i in range(N):
        for j in range(i, N):
            d = {k: 0 for k in range(N)}       # put 0 photons in each mode (k)
            d[i] += 1
            d[j] += 1
            # load ith mode and jth mode with 1 photon
            fs = fock_state_on(d)
            ev, evec = closest_state_to(fs)
            chi = (ev - (f1s[i] + f1s[j]))
            chis[i][j] = chi
            chis[j][i] = chi

            if chi_prime:
                d[j] += 1
                fs = fock_state_on(d)
                ev, evec = closest_state_to(fs)
                chip = (ev - (f1s[i] + 2*f1s[j]) - 2 * chis[i][j])
                chips[i][j] = chip
                chips[j][i] = chip

    if chi_prime:
        return np.array(f1s), np.array(chis), np.array(chips), np.array(fzpfs), np.array(f0s)
    else:
        return np.array(f1s), np.array(chis), np.array(fzpfs), np.array(f0s)


def black_box_hamiltonian_nq(freqs, zmat, ljs, cos_trunc=6, fock_trunc=8, show_fit=False):
    """
    N-Qubit version of bbq, based on the full Z-matrix
    Currently reproduces 1-qubit data, untested on n-qubit data
    Assume: Solve the model without loss in HFSS.
    """
    nf = len(freqs)
    nj = len(ljs)
    assert zmat.shape == (nf, nj, nj)

    imY = (1/zmat[:, 0, 0]).imag
    # zeros where the sign changes from negative to positive

    (zeros,) = np.where((imY[:-1] <= 0) & (imY[1:] > 0))
    nz = len(zeros)

    imYs = np.array([1 / zmat[:, i, i] for i in range(nj)]).imag
    f0s = np.zeros(nz)
    slopes = np.zeros((nj, nz))
    import matplotlib.pyplot as plt
    # Fit a second order polynomial in the region around the zero
    # Extract the exact location of the zero and the associated slope
    # If you need better than second order fit, you're not sampling finely enough
    for i, z in enumerate(zeros):
        f0_guess = (freqs[z+1] + freqs[z]) / 2
        zero_polys = np.polyfit(
            freqs[z-1:z+3] - f0_guess, imYs[:, z-1:z+3].transpose(), 2)
        zero_polys = zero_polys.transpose()
        f0s[i] = f0 = min(np.roots(zero_polys[0]),
                          key=lambda r: abs(r)) + f0_guess
        for j, p in enumerate(zero_polys):
            slopes[j, i] = np.polyval(np.polyder(p), f0 - f0_guess)
        if show_fit:
            plt.plot(freqs[z-1:z+3] - f0_guess, imYs[:, z-1:z +
                                                     3].transpose(), lw=1, ls='--', marker='o', label=str(f0))
            p = np.poly1d(zero_polys[0, :])
            p2 = np.poly1d(zero_polys[1, :])
            plt.plot(freqs[z-1:z+3] - f0_guess, p(freqs[z-1:z+3] - f0_guess))
            plt.plot(freqs[z-1:z+3] - f0_guess, p2(freqs[z-1:z+3] - f0_guess))
            plt.legend(loc=0)

    zeffs = 2 / (slopes * f0s[np.newaxis, :])
    # Take signs with respect to first port
    zsigns = np.sign(zmat[zeros, 0, :])
    fzpfs = zsigns.transpose() * np.sqrt(hbar * abs(zeffs) / 2)

    H = black_box_hamiltonian(f0s, ljs, fzpfs, cos_trunc, fock_trunc)
    return make_dispersive(H, fock_trunc, fzpfs, f0s)

black_box_hamiltonian_nq = black_box_hamiltonian_nq
