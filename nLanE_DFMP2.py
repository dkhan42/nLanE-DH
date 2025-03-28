import numpy as np
from numpy import sqrt
from pyscf import dft
from pyscf.mp.dfmp2_native import DFRMP2
from pyscf.mp.dfump2_native import DFUMP2 
from scipy.integrate import quad
'''
THIS SCRIPT HAS BEEN WRITTEN BY GPT O3-MINI BASED ON THE nLanE.py SCRIPT. USER DISCRETION IS ADVISED.
'''

def eval_xc(mol, mf, xc):
    '''
    Evaluates the exchange-correlation energy of a given functional on a pre-converged mol object.
    '''
    mf_nscf = dft.KS(mol, xc=xc)
    mf_nscf.grids.level = 9 
    # Evaluate total energy using the provided density matrix
    e_nscf = mf_nscf.energy_tot(dm=mf.make_rdm1(ao_repr=True))
    return mf_nscf.scf_summary['exc']


def adiabatic_connection(lam, a, b, c):
    '''
    Analytical and convex adiabatic connection ansatz with correct lam = ∞ expansion.
    '''
    return a + ((b * sqrt(lam + 1)) / (c * lam + 1))


def solve_abc(Exx, Ec_MP2, W1):
    '''
    Returns the values of a, b, c in the above adiabatic connection ansatz by imposing the constraints:
      W(0) = Exx (HF exchange energy)
      W'(0) = 2*Ec_MP2 (Twice the MP2 correlation energy)
      W(1) = SCAN exchange energy + 2*SCAN correlation energy
    '''
    if np.abs(Ec_MP2) < 1e-6:  # For a one-electron system, set b = 0 and c = 0.
        a = Exx
        b, c = 0, 0
    else:
        alpha = (W1 - Exx) / (2 * Ec_MP2)
        c = sqrt(9 * alpha**2 - 16 * sqrt(2) * alpha + 12 * alpha + 4) - alpha + 2  # Analytical solution.
        c = c / (4 * alpha)
        b = 2 * Ec_MP2 / (0.5 - c)
        a = Exx - b
        if a > W1:
            c = -c
            b = 2 * Ec_MP2 / (0.5 - c)
            a = Exx - b
    return a, b, c


def return_energy(mol, xcscf='R2SCAN', W1x = 'SCAN', W1c = 'SCAN', integral='numerical'):
    '''
    Returns the energy of a given molecule using the double-hybrid (DH) functional evaluated on SCAN density and orbitals.
    MP2 correlation energy (from native DF-MP2) is used for the W'(0) slope constraint.
    
    In this implementation the SCF and MP2 steps use density fitting.
    For the MP2 part, we use the “native” DF-MP2 implementation (DFRMP2 for RHF references).
    Inputs :
    mol - mol object from pyscf GTO
    xcscf - Choice of functional to run the SCF. All semi-local functionals seem to provide similar results.
    W1x - Choice of exchange functional for the W1 limit. SCAN is recommended
    W1c - Choice of correlation functional for the W1 limit. SCAN is recommended
    integral - evaluate nLanE Exc analytically or numerically via scipy.quad
    '''
    # Use DFT KS object with density fitting.
    mf_scf = dft.KS(mol)
    mf_scf.xc = xcscf
    mf_scf.grids.atom_grid = (99, 590)
    if xcscf=="SCAN":
        mf_scf = mf_scf.newton() #SCAN SCFs are hard to converge
    # Enable density fitting for the SCF step.
    mf_scf = mf_scf.density_fit()
    energy = mf_scf.kernel()

    if not mf_scf.converged:
        raise ValueError("SCF not converged. Try reducing tolerance or use R2SCAN/PBE for SCF.")

    Exc = mf_scf.scf_summary['exc']
    # Evaluate HF exchange energy on the SCAN orbitals.
    Exx = eval_xc(mol, mf_scf, 'HF,')
    
    # Run DF-MP2 using the native implementation.
    if mf_scf.mol.spin == 0:
        dfmp2 = DFRMP2(mf_scf)
    else:
        dfmp2 = DFUMP2(mf_scf)

    dfmp2.kernel()
    Ec_MP2 = dfmp2.e_corr  # MP2 correlation energy.

    # Evaluate SCAN (or the chosen functional) correlation energy.
    Ec = eval_xc(mol,mf_scf,f',{W1c}') #SCAN correlation energy

    if xcscf=='SCAN':
        Ex = Exc - Ec #SCAN exchange energy
    else:
        Ex = eval_xc(mol,mf_scf,f'{W1x},') #SCAN exchange energy
    
    W1 = Ex + 2 * Ec #W1 limit approximation using SCAN
    a, b, c = solve_abc(Exx, Ec_MP2, W1)

    if integral == 'numerical':
        # Numerical integration of the adiabatic connection.
        Exc_ac, _ = quad(lambda lam: adiabatic_connection(lam, a, b, c), 0, 1, limit=1000)
    else:
        # Analytical integration (beware of potential numerical instabilities with np.arctan).
        Exc_ac = a + b * (((2 * np.sqrt(2) - 2) / c) +
                          (2 * np.sqrt(1 - c) / (c**1.5)) *
                          (np.arctan(np.sqrt(c / (1 - c))) - np.arctan(np.sqrt(2 * c / (1 - c)))))
    return energy - Exc + Exc_ac
