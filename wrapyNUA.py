"""
This is a pre- and postprocessing script for the Python wrapper around the 
Fortran numerical uncertainty core

Written by mrentschler, last modified 2024-12-02
"""

#%% Settings

path = './Tutorial_examples/steady/'
dataFileName = 'steady.dat'

# 0: Steady, 1: Unsteady
isUnsteady = 0

# -1: Relative grid step size, 1: Number of cells in 1D, 2: 2D, 3: 3D
gridStepSizeMethod = 1

# -1: Relative time step size, 1: Number of time steps
timeStepSizeMethod = -1

# Reference solution. 1: Value of finest discretization, 2: Second finest, ...
uncertaintySolution = 1

# 0: Plots only value of reference solution, 1: Plots all
showAllUncertainties = 1

# Text style. 0: Sans-serif, 1: Serif
serif = 0


#%% Examples
"""
Tutorial
--------
steady: isUnsteady = 0, gSSM = 1, tSSM = -1, uncertaintySolution = 1
unsteady: isUnsteady = 1, gSSM = 2, tSSM = -1, uncertaintySolution = 1

MARIN
-----
std-2d-structured_1: isUnsteady = 0, gSSM = 1, tSSM = -1, uncertaintySolution = 1
std-3d-structured_1: isUnsteady = 0, gSSM = 3, tSSM = -1, uncertaintySolution = 1
unstd-2d-structured_1: isUnsteady = 1, gSSM = 2, tSSM = 1, uncertaintySolution = 1
unstd-2d-unstructured_1: isUnsteady = 1, gSSM = 2, tSSM = 1, uncertaintySolution = 1
unstd-3d-unstructured_1: isUnsteady = 1, gSSM = -1, tSSM = -1, uncertaintySolution = 2
unstd-3d-unstructured_2: isUnsteady = 1, gSSM = 3, tSSM = 1, uncertaintySolution = 1
unstd-3d-unstructured_3: isUnsteady = 1, gSSM = 3, tSSM = 1, uncertaintySolution = 1
"""

#%% Imports

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from wrap import NUA


#%% Functions

def fit_p(r_i, phi_0, alpha, p):
    return (phi_0 + alpha*r_i**p)


def fit_p12(r_i, phi_0, alpha_1, alpha_2):
    return (phi_0 + alpha_1*r_i + alpha_2*r_i**2)


def fit_xt(r_i, tau_i, phi_0, alpha_x, p_x, alpha_t, p_t):    
    return (phi_0 + alpha_x*r_i**p_x + alpha_t*tau_i**p_t)


def fit_x12t(r_i, tau_i, phi_0, alpha_x1, alpha_x2, alpha_t, p_t):
    return (phi_0 + alpha_x1*r_i + alpha_x2*r_i**2 + alpha_t*tau_i**p_t)


def fit_xt12(r_i, tau_i, phi_0, alpha_x, p_x, alpha_t1, alpha_t2):    
    return (phi_0 + alpha_x*r_i**p_x + alpha_t1*tau_i + alpha_t2*tau_i**2)


def fit_x12t12(r_i, tau_i, phi_0, alpha_x1, alpha_x2, alpha_t1, alpha_t2):    
    return (phi_0 + alpha_x1*r_i + alpha_x2*r_i**2 + alpha_t1*tau_i + alpha_t2*tau_i**2)


def getData(path, dataFileName, isUnsteady, rMethod, tauMethod):
    # Read input file
    df = pd.read_csv(f'{path}/{dataFileName}', sep=r'\s+', comment='#')
    
    # Extract experimental measurements
    df_exp = df[df.iloc[:, 0] == 0]
    if len(df_exp) > 1:
        raise Exception('Only one line with experimental values can be specified!')
    
    # Extract simulation data (without experimental value and number of solutions)
    df_sim = df[df.iloc[:, 0] != 0].copy()
    
    # Calculate relative grid step size
    if rMethod == -1:
        r_i = df_sim.iloc[:, 0] / df_sim.iloc[:, 0].min()
    elif rMethod == 1 or rMethod == 2 or rMethod == 3:
        r_i = (df_sim.iloc[:, 0].max() / df_sim.iloc[:, 0])**(1 / rMethod)
    df_sim.insert(0, 'r_i', r_i)
    
    if isUnsteady:
        # Calculate relative time step size
        if tauMethod == -1:
            tau_i = df_sim.iloc[:, 2] / df_sim.iloc[:, 2].min()
        elif tauMethod == 1:
            tau_i = df_sim.iloc[:, 2].max() / df_sim.iloc[:, 2]
        df_sim.insert(1, 'tau_i', tau_i)
        
        # Sort in ascending order
        df_sim = df_sim.sort_values(['r_i', 'tau_i'])
        
        # Get variable names
        varNames = list(df_sim)[4:]
        
    else:
        # Sort in ascending order
        df_sim = df_sim.sort_values('r_i')
        
        # Get variable names
        varNames = list(df_sim)[2:]
    
    return varNames, df_sim, df_exp


def plotUncertainty_std(path, var, r_i, phi_i, uphi_i, p, alpha, phi_0, 
                        iFitType, exp, uSol, showAll, serif):
    # Text style
    if serif == 1:
        plt.rc('font', family='Serif', size=16)
        plt.rcParams['mathtext.fontset'] = 'stix'
    else:
        plt.rc('font', family='DejaVu Sans', size=16)
        plt.rcParams['mathtext.fontset'] = 'dejavusans'
    
    # Alphanumeric variable name
    var_an = ''.join(letter for letter in var if letter.isalnum())
    
    # Independent variable
    r_lin = np.linspace(0, np.max(r_i), num=1000)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(8, 6), dpi=300, tight_layout=True)
    
    # Data points
    ax.plot(r_i, phi_i, 'ko', fillstyle='none')
    # Error bars
    ax.errorbar(r_i, phi_i, uphi_i, color='black', fmt='none', capsize=5)
    # Error values as percentage
    uper_i = np.abs(uphi_i / phi_i[uSol - 1] * 100)
    if not showAll:
        r_i = [r_i[uSol - 1]]
        phi_i = [phi_i[uSol - 1]]
        uphi_i = [uphi_i[uSol - 1]]
        uper_i = [uper_i[uSol - 1]]
    for _r, _phi, _uphi, _uper in zip(r_i, phi_i, uphi_i, uper_i):
        ax.text(_r, _phi + _uphi, f'{_uper:.1f}%', ha='center', va='bottom')
    """ Alternative code to plot the percentages alternately above and below the bars
    for _r, _phi, _uphi, _uper, _sign, _va in zip(r_i, phi_i, uphi_i, uper_i, 
        np.resize([1, -1], len(r_i)), np.resize(['bottom', 'top'], len(r_i))):
        ax.text(_r, _phi + _sign * _uphi, f'{_uper:.1f}%', ha='center', va=_va) """
    
    # Experiment if available and other value than 0 (default placeholder value)
    if exp:
        ax.plot(0, exp, 'ks', clip_on=False, markerfacecolor='none', 
                markersize=10, label='Experiment')
    
    # Fit
    if iFitType == 3:
        # Two-term expansion
        ax.plot(r_lin, fit_p12(r_lin, phi_0, alpha, p), 'k-', label='Fit p=*1,2')
    else:
        # One-term expansion
        ax.plot(r_lin, fit_p(r_lin, phi_0, alpha, p), 'k-', label=f'Fit p={p:.3g}')
    
    ax.set_xlim(left=0)
    ax.set_ylabel(r'$%s$' % var)
    ax.set_xlabel(r'Relative grid step size')
    ax.legend()
    ax.grid()
    plt.savefig(f'{path}/Uncertainty_{var_an}.png')


def plotUncertainty_unstd(path, var, r_i, tau_i, phi_i, uphi_i, fp, iFitType, 
                          exp, uSol, showAll, serif):
    # Text style
    if serif == 1:
        plt.rc('font', family='Serif', size=16)
        plt.rcParams['mathtext.fontset'] = 'stix'
    else:
        plt.rc('font', family='DejaVu Sans', size=16)
        plt.rcParams['mathtext.fontset'] = 'dejavusans'
    
    # Alphanumeric variable name
    var_an = ''.join(letter for letter in var if letter.isalnum())
    
    # Independent variables
    r_lin = np.linspace(0.5, np.max(r_i), num=1000)
    tau_lin = np.linspace(0.5, np.max(tau_i), num=1000)
    r_lin, tau_lin = np.meshgrid(r_lin, tau_lin)
    
    # Create figure
    fig = plt.figure(figsize=(8, 6), dpi=300, tight_layout=True)
    ax = plt.axes(projection='3d')
    
    # Data points
    ax.scatter(r_i, tau_i , phi_i, facecolors='none', edgecolors='k')
    # Error bars
    ax.errorbar(r_i, tau_i, phi_i, uphi_i, ecolor='black', fmt='none', capsize=5)
    # Error values as percentage
    uper_i = np.abs(uphi_i / phi_i[uSol - 1] * 100)
    if not showAll:
        r_i = [r_i[uSol - 1]]
        tau_i = [tau_i[uSol - 1]]
        phi_i = [phi_i[uSol - 1]]
        uphi_i = [uphi_i[uSol - 1]]
        uper_i = [uper_i[uSol - 1]]
    for _r, _tau, _phi, _uphi, _uper in zip(r_i, tau_i, phi_i, uphi_i, uper_i):
        ax.text(_r, _tau, _phi + _uphi, f'{_uper:.1f}%', size=14, 
                ha='center', va='bottom')
    """ Alternative code to plot the percentages below the error bars
    for _r, _tau, _phi, _uphi, _uper in zip(r_i, tau_i, phi_i, uphi_i, uper_i):
        ax.text(_r, _tau, _phi - _uphi, f'{_uper:.1f}%', size=14, 
                ha='center', va='top') """
    
    # Fit
    if iFitType in [3, 11, 12]:
        # Two-term expansion in space, one-term in time
        phi_fit = fit_x12t(r_lin, tau_lin, fp[0], fp[1], fp[3], fp[2], fp[4])
        label = r'Fit with two terms in space and one in time, $p_t=%.3g$' % fp[4]
    elif iFitType in [6, 13, 14]:
        # One-term expansion in space, two-term in time
        phi_fit = fit_xt12(r_lin, tau_lin, fp[0], fp[1], fp[3], fp[2], fp[4])
        label = r'Fit with one term in space and two in time, $p_x=%.3g$' % fp[3]
    elif iFitType == 15:
        # Two-term expansion in space and time
        phi_fit = fit_x12t12(r_lin, tau_lin, fp[0], fp[1], fp[3], fp[2], fp[4])
        label = r'Fit with two terms in both space and time'
    else:
        # One-term expansion in space and time
        phi_fit = fit_xt(r_lin, tau_lin, fp[0], fp[1], fp[3], fp[2], fp[4])
        label = r'Fit with one term in both space and time, $p_x=%.3g$, $p_t=%.3g$' \
            % (fp[3], fp[4])
    surf = ax.contour3D(r_lin, tau_lin, phi_fit, 50, cmap='viridis', alpha=0.6, 
                        zorder=0)
    fig.colorbar(surf, shrink=0.5, aspect=20, location='right')
    
    ax.set_title(label, fontsize=14)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel('Relative grid step size')
    ax.set_ylabel('Relative time step size')
    ax.set_zlabel(r'$%s$' % var)
    ax.view_init(30, -135)
    plt.savefig(f'{path}/Uncertainty_{var_an}.png')


#%% Main

# Read data from file
varNames, df_sim, df_exp = getData(path, dataFileName, isUnsteady,
                                   gridStepSizeMethod, timeStepSizeMethod)    

# Initialize NUA object
nua = NUA()
print(f'Uncertainty tools version: {nua.toolversion()}\n')

for var in varNames:
    # Alphanumeric variable name
    var_an = ''.join(letter for letter in var if letter.isalnum())
    print(f'Fitting {var_an} ...')
    
    # Relative grid step size
    r_i = df_sim['r_i'].to_numpy()
    
    # Solutions of current variable
    phi_i = df_sim[var].to_numpy()
    
    # Experimental measurement
    try:
        exp = df_exp[var].to_numpy()[0]
    except:
        exp = 0
    
    if isUnsteady:
        # Relative time step size
        tau_i = df_sim['tau_i'].to_numpy()
        
        # Call unsteady uncertainty routine
        [uphi_i, fp, iFitType] = \
            nua.numerical_uncertainty_unstd(r_i, tau_i, phi_i)
        
        if iFitType == -1:
            raise Exception('... Error in uncertainty analysis')
            break
        elif iFitType in [3, 11, 12]:
            # Two-term expansion in space, one-term in time
            print(f'... iFitType={iFitType}, phi_0={fp[0]:.3g}, '
                  + f'alpha_x1={fp[1]:.3g}, alpha_x2={fp[3]:.3g}, '
                  + f'alpha_t={fp[2]:.3g}, p_t={fp[4]:.3g}')
        elif iFitType in [6, 13, 14]:
            # One-term expansion in space, two-term in time
            print(f'... iFitType={iFitType}, phi_0={fp[0]:.3g}, '
                  + f'alpha_x={fp[1]:.3g}, p_x={fp[3]:.3g}, '
                  + f'alpha_t1={fp[2]:.3g}, alpha_t2={fp[4]:.3g}')
        elif iFitType == 15:
            # Two-term expansion in space and time
            print(f'... iFitType={iFitType}, phi_0={fp[0]:.3g}, '
                  + f'alpha_x1={fp[1]:.3g}, alpha_x2={fp[3]:.3g}, '
                  + f'alpha_t1={fp[2]:.3g}, alpha_t2={fp[4]:.3g}')
        else:
            # One-term expansion in space and time
            print(f'... iFitType={iFitType}, phi_0={fp[0]:.3g}, '
                  + f'alpha_x={fp[1]:.3g}, alpha_t={fp[2]:.3g}, '
                  + f'p_x={fp[3]:.3g}, p_t={fp[4]:.3g}')
        
        # Uncertainty plot
        plotUncertainty_unstd(path, var, r_i, tau_i, phi_i, uphi_i, fp, 
                              iFitType, exp, uncertaintySolution, 
                              showAllUncertainties, serif)
    
    else:
        # Call steady uncertainty routine
        [uphi_i, p, alpha, phi_0, iFitType] = \
            nua.numerical_uncertainty_std(r_i, phi_i, len(r_i))
        
        if iFitType == -1:
            raise Exception('Error in uncertainty analysis, variable {var_an}')
        elif iFitType == 3:
            # Two-term expansion
            print(f'... iFitType={iFitType}, phi_0={phi_0:.3g}, '
                  + f'alpha_1={alpha:.3g}, alpha_2={p:.3g}')
        else:
            # One-term expansion
            print(f'... iFitType={iFitType}, phi_0={phi_0:.3g}, '
                  + f'alpha={alpha:.3g}, p={p:.3g}')
        
        # Uncertainty plot
        plotUncertainty_std(path, var, r_i, phi_i, uphi_i, p, alpha, phi_0, 
                            iFitType, exp, uncertaintySolution, 
                            showAllUncertainties, serif)

