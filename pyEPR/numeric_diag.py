'''
    Written by Phil Rheinhold
    Updated by Zlatko Minev & Lysander Christakis

    This file is tricky, use caution to modify.
'''

from __future__ import print_function    # Python 2.7 and 3 compatibility
import numpy as np

from qutip import basis, tensor
from functools import  reduce

from .toolbox import fluxQ, fact
from scipy.constants import pi, hbar, e, h


__all__ = ['bbq_hmt', 'make_dispersive', 'bbq_nq']

def dot(ais, bis):
    return sum(ai*bi for ai, bi in zip(ais, bis))

def cos_approx(x, cos_trunc=5):
    return sum((-1)**i * x**(2*i) / float(fact(2*i)) for i in range(2, cos_trunc + 1))

def bbq_hmt(fs, ljs, fzpfs, cos_trunc=5, fock_trunc=8, individual = False):
    """
    :param fs: Zero crossing frequencies in Hz, length N
    :param ljs: qubit junction inductances in Henries, length M
    :param fzpfs: Estimated zero-point fluxes for each mode across each junction, shape MxN
    :return: Hamiltonian in units of Hz (i.e H / h)
    All in SI units. The ZPF fed in are the generalized, not reduced flux.
    """
    import qutip
    nmodes = len(fs)
    nqubits = len(ljs)
    fs, ljs, fzpfs = map(np.array, (fs, ljs, fzpfs))
    ejs = fluxQ**2 / ljs
    fjs = ejs / h ;

    assert np.isnan(fzpfs).any() == False, "Phi ZPF has NAN, this is NOT allowed! Fix me. \n%s" %fzpfs
    assert np.isnan(ljs).any() == False, "Ljs has NAN, this is NOT allowed! Fix me."
    assert np.isnan(fs).any() == False, "freqs has NAN, this is NOT allowed! Fix me."
    assert fzpfs.shape == (nqubits, nmodes), "incorrect shape for zpf array, {} not {}".format(fzpfs.shape, (nqubits, nmodes))
    assert fs.shape == (nmodes,), "incorrect number of mode frequencies"
    assert ejs.shape == (nqubits,), "incorrect number of qubit frequencies"

    def tensor_out(op, loc):
        "Make operator <op> tensored with identities at locations other than <loc>"
        op_list = [qutip.qeye(fock_trunc) for i in range(nmodes)]
        op_list[loc] = op
        return reduce(qutip.tensor, op_list)

    a = qutip.destroy(fock_trunc)
    ad = a.dag()
    n = qutip.num(fock_trunc)
    mode_fields = [tensor_out(a + ad, i) for i in range(nmodes)]
    mode_ns = [tensor_out(n, i) for i in range(nmodes)]


    def cos(x):
        return cos_approx(x, cos_trunc=cos_trunc)

    linear_part = dot(fs, mode_ns)
    cos_interiors = [dot(fzpf_row/fluxQ, mode_fields) for fzpf_row in fzpfs]
    nonlinear_part = dot(-fjs, map(cos, cos_interiors))
    if individual:
        return linear_part, nonlinear_part
    else:
        return linear_part + nonlinear_part


def make_dispersive(H, fock_trunc, fzpfs= None, f0s = None, chi_prime = False,
                    use_1st_order = False):
    """
    Takes hamiltonian, and finds second order dispersive approximation
    I.e. returns final qubit freqs, chis, chi prime, phi_zpf flux (not reduced), and linear ferq
        Note, the diagonal chis are the alphas directly.
    """
    import qutip
    if hasattr(H, '__len__'): # is it an array / list?
        [H_lin, H_nl] = H
        H = H_lin + H_nl
    else:  # make sure its a quanutm object
        assert type(H) == qutip.qobj.Qobj, "Please pass in either  a list of Qobjs or Qobj for the Hamiltonian"

    evals, evecs = H.eigenstates()
    evals -= evals[0]

    N = int(np.log(H.shape[0]) / np.log(fock_trunc))    # number of modes
    assert H.shape[0] == fock_trunc ** N

    def fock_state_on(d):
        return qutip.tensor(*[qutip.basis(fock_trunc, d.get(i, 0)) for i in range(N)])  # give me the value d[i]  or 0 if d[i] does not exist


    if use_1st_order:
        num_modes = N
        print( "Using 1st O")

        def multi_index_2_vector(d, num_modes, fock_trunc):
            return tensor([basis(fock_trunc, d.get(i, 0)) for i in range(num_modes)])
            '''this function creates a vector representation a given fock state given the data for excitations per
                        mode of the form d={mode number: # of photons}'''

        def find_multi_indices(fock_trunc):
            multi_indices =[{ind:item for ind,item in enumerate([i,j,k])} for i in range(fock_trunc)
            for j in range(fock_trunc)
                 for k in range(fock_trunc)]
            return multi_indices
            '''this function generates all possible multi-indices for three modes for a given fock_trunc'''

        def get_expect_number(left, middle, right):
            return (left.dag()*middle*right).data.toarray()[0,0]
            '''this function calculates the expectation value of an operator called "middle" '''

        def get_basis0(fock_trunc, num_modes):
            multi_indices = find_multi_indices(fock_trunc)
            basis0        = [ multi_index_2_vector(multi_indices[i],num_modes,fock_trunc) for i in range(len(multi_indices)) ]
            evalues0      = [ get_expect_number(v0,H_lin, v0).real for v0 in basis0 ]
            return multi_indices, basis0, evalues0
            '''this function creates a basis of fock states and their corresponding eigenvalues'''

        def closest_state_to(vector0):

            def PT_on_vector(original_vector, original_basis, pertub, energy0, evalue):
                new_vector = 0 * original_vector;
                for i in range(len(original_basis)):
                    if (energy0[i]-evalue) > 1e-3:
                        new_vector += ((original_basis[i].dag()*H_nl*original_vector).data.toarray()[0,0])*original_basis[i]/(evalue-energy0[i])
                    else:
                        pass
                return (new_vector + original_vector)/(new_vector + original_vector).norm()
                '''this function calculates the normalized vector with the first order correction term
                   from the non-linear hamiltonian '''

            [multi_indices,basis0, evalues0] = get_basis0(fock_trunc, num_modes)
            evalue0 = get_expect_number(vector0, H_lin, vector0)
            vector1 = PT_on_vector(vector0, basis0, H_nl, evalues0,evalue0)

            index = np.argmax([(vector1.dag() * evec).norm() for evec in evecs])
            return evals[index], evecs[index]

    else:
        def closest_state_to(s):
            def distance(s2):
                return (s.dag() * s2[1]).norm()
            return max(zip(evals, evecs), key=distance)

    f1s = [closest_state_to(fock_state_on({i:1}))[0] for i in range(N)]
    chis = [[0]*N for _ in range(N)]
    if chi_prime: chips = [[0]*N for _ in range(N)]
    for i in range(N):
        for j in range(i, N):
            d = {k: 0 for k in range(N)}       # put 0 photons in each mode (k)
            d[i] += 1
            d[j] += 1
            fs = fock_state_on(d)              # load ith mode and jth mode with 1 photon
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


def bbq_nq(freqs, zmat, ljs, cos_trunc=6, fock_trunc=8, show_fit = False):
    """
    N-Qubit version of bbq, based on the full Z-matrix
    Currently reproduces 1-qubit data, untested on n-qubit data
    Assume: Solve the model without loss in HFSS.
    """
    nf = len(freqs)
    nj = len(ljs)
    assert zmat.shape == (nf, nj, nj)


    imY = (1/zmat[:,0,0]).imag
    # zeros where the sign changes from negative to positive

    (zeros,) = np.where((imY[:-1] <= 0) & (imY[1:] > 0))
    nz = len(zeros)

    imYs = np.array([1 / zmat[:, i, i] for i in range(nj)]).imag
    f0s = np.zeros(nz)
    slopes = np.zeros((nj, nz))
    import matplotlib.pyplot as plt
    # Fit a second order polynomial in the region around the zero
    # Extract the exact location of the zero and the assocated slope
    # If you need better than second order fit, you're not sampling finely enough
    for i, z in enumerate(zeros):
        f0_guess = (freqs[z+1] + freqs[z]) / 2
        zero_polys = np.polyfit(freqs[z-1:z+3] - f0_guess, imYs[:, z-1:z+3].transpose(), 2)
        zero_polys = zero_polys.transpose()
        f0s[i] = f0 = min(np.roots(zero_polys[0]), key=lambda r: abs(r)) + f0_guess
        for j, p in enumerate(zero_polys):
            slopes[j, i] = np.polyval(np.polyder(p), f0 - f0_guess)
        if show_fit:
            plt.plot(freqs[z-1:z+3] - f0_guess, imYs[:, z-1:z+3].transpose(), lw=1, ls ='--', marker = 'o', label =str(f0))
            p = np.poly1d(zero_polys[0,:]); p2 = np.poly1d(zero_polys[1,:])
            plt.plot(freqs[z-1:z+3] - f0_guess, p(freqs[z-1:z+3] - f0_guess)); plt.plot(freqs[z-1:z+3] - f0_guess, p2(freqs[z-1:z+3] - f0_guess));   plt.legend(loc = 0)


    zeffs = 2 / (slopes * f0s[np.newaxis, :])
    zsigns = np.sign(zmat[zeros, 0, :]) # Take signs with respect to first port
    fzpfs = zsigns.transpose() * np.sqrt(hbar * abs(zeffs) / 2)

    H = bbq_hmt(f0s, ljs, fzpfs, cos_trunc, fock_trunc)
    return make_dispersive(H, fock_trunc, fzpfs, f0s)
