#!/usr/local/bin/python3

# script which runs molConv varying T/Tc or T/TF for different molecules and distance measures 
# this is used to generate nice plots and to generate test cases for comparison
# - in terminal (cmd.exe) cd into folder with molConv executable
# - create a temporary folder in local directory, typically called "./tmp/", see "folder" variable. 
#   all output files and folders are generated in this folder. 
#   this works also on Windows, so no absolute path is needed or "\\".
# - run "python molConv_run.py" or "python3 molConv_run.py"
# last change 22/9/2026 by Andi

# TODO: not all figures have been updated after changes and might give missing keys errors.
#       use figure = None as reference which is working. 

import numpy as np
from mpmath import polylog, re, im
import subprocess
import os
from datetime import datetime

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as mticker
mpl.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 2
mpl.rcParams["savefig.directory"] = os.path.dirname(__file__) # save to current directy

#colors = cm.get_cmap('tab20')

# path to molConv executable
if os.name == 'nt':
    molConv_path = r"./molConv_VisualStudio/x64/Release/molConv"
else:
    molConv_path = r"./molConv"

# constants
hPlanck = 6.62607015e-34
hbar    = 1.054571817e-34
kB      = 1.380649e-23
mn      = 1.6749e-27
k_nK    = (1.054571817/1.380649)*1.0e-2   # get Tcrit and Tc in nK out of hbar*omega: k_nK = (hbar*1e9)/kB;
k_Epot  = (1.67492749804/1.380649)*1e-7   # get [Epot] = amu*(rad/s*mu)^2 -> nK; mn*(1e-6)^2/kB*1e9 = mn/kB*1e-3;
k_Ekin  = (1.67492749804/1.380649)*1e-1   # get [Ekin] = amu*(mm/s)^2     -> nK; mn*(1e-3)^2/kB*1e9 = mn/kB*1e3;
zeta    = float(re(polylog(3,1.0)))

# set output options for printing numpy arrays
# we print with 6 digits and %13 includes always space for optional negative sign
# use %13s for header
np.set_printoptions(precision=6, linewidth=150, formatter={'float': lambda x: '%13.6e'%x})

################################################################################################
# predefined molecules and parameters
################################################################################################
    
K40K40 = { 
    'label'     : 'K40K40',         # label used for result directory and plotting
    'species'   : ['K40','K40'],    # name per species
    'mass'      : [40   , 40  ],    # mass per species in amu
    'atom_stat' : ['FermiDirac','FermiDirac'], # statistics of both species
    'mol_stat'  : 'MaxwellBoltzmann',          # statistics of molecule
}

K39K39 = { 
    'label'     : 'K39K39',         # label used for result directory and plotting
    'species'   : ['K39','K39'],    # name per species
    'mass'      : [39   , 39  ],    # mass per species in amu
    'atom_stat' : ['BoseEinstein','BoseEinstein'], # statistics of both species
    'mol_stat'  : 'MaxwellBoltzmann',              # statistics of molecule
}

K40K39 = { 
    'label'     : 'K40K39',         # label used for result directory and plotting
    'species'   : ['K40','K39'],    # name per species
    'mass'      : [40   , 39  ],    # mass per species in amu
    'atom_stat' : ['FermiDirac','BoseEinstein'], # statistics of both species
    'mol_stat'  : 'MaxwellBoltzmann',            # statistics of molecule
}

K39K41 = { 
    'label'     : 'K39K41',         # label used for result directory and plotting
    'species'   : ['K39','K41'],    # name per species
    'mass'      : [39   , 41  ],    # mass per species in amu
    'atom_stat' : ['BoseEinstein','BoseEinstein'], # statistics of both species
    'mol_stat'  : 'MaxwellBoltzmann',            # statistics of molecule
}


#U_ratio = 0.5                       # trap depth ratio U_Cr/U_Li
Li6Cr53 = { 
    'label'     : 'Li6Cr53',        # label used for result directory and plotting
    'species'   : ['Li6','Cr53'],   # name per species
    'mass'      : [6    , 53   ],   # mass per species in amu
    'atom_stat' : ['FermiDirac','FermiDirac'], # statistics of both species
    'mol_stat'  : 'MaxwellBoltzmann',          # statistics of molecule
    #'N'         : [10000, 10000],   # atom number per species
    #'T'         : [150  , 150],     # temperature in nK (might be overwritten)
    #'f_rad'     : [100, 100*np.sqrt(U_ratio*6/53)], # radial (x) trap frequency in Hz
    #'f_vert'    : [100, 100*np.sqrt(U_ratio*6/53)], # vertical (y) trap frequencies in Hz
    #'f_ax'      : [ 15,  10],   # axial (z) trap frequency in Hz
}

Li6Cr52 = { 
    'label'     : 'Li6Cr52',        # label used for result directory and plotting
    'species'   : ['Li6','Cr52'],   # name per species
    'mass'      : [6    , 52   ],   # mass per species in amu
    'atom_stat' : ['FermiDirac','BoseEinstein'], # statistics of both species
    'mol_stat'  : 'MaxwellBoltzmann',            # statistics of molecule
    #'N'         : [10000, 10000],   # atom number per species
    #'T'         : [150  , 150],     # temperature in nK (might be overwritten)
    #'f_rad'     : [100, 100*np.sqrt(U_ratio*6/53)], # radial (x) trap frequency in Hz
    #'f_vert'    : [100, 100*np.sqrt(U_ratio*6/53)], # vertical (y) trap frequencies in Hz
    #'f_ax'      : [ 15,  10],       # axial (z) trap frequency in Hz
}

Li6Cr50 = {
    'label'     : 'Li6Cr50',        # label used for result directory and plotting
    'species'   : ['Li6','Cr50'],   # name per species
    'mass'      : [6    , 50   ],   # mass per species in amu
    'atom_stat' : ['FermiDirac','BoseEinstein'], # statistics per species
    'mol_stat'  : 'MaxwellBoltzmann',            # statistics of molecule
    #'N'         : [10000, 10000],   # atom number per species
    #'T'         : [150  , 150],     # temperature in nK (might be overwritten)
    #'f_rad'     : [100, 100*np.sqrt(U_ratio*6/53)], # radial (x) trap frequency in Hz
    #'f_vert'    : [100, 100*np.sqrt(U_ratio*6/53)], # vertical (y) trap frequencies in Hz
    #'f_ax'      : [ 15,  10],       # axial (z) trap frequency in Hz
}

# available distance measures and corresponding labels for plotting
distance_measure_all = {'delta_x*delta_p'           : r'$\Delta x \Delta p$',
                        'delta_x*delta_v'           : r'$\Delta x \Delta p$',
                        'delta_x'                   : r'$\Delta x$',
                        'delta_v'                   : r'$\Delta v$',
                        'cross(delta_x^,delta_p^)'  : r'$|\Delta \vec{x} \times \Delta \vec{p}$|', 
                        '<delta_x*delta_p>'         : r'$\langle \Delta x \Delta p \rangle$',
                        'max|delta_xi*delta_pi|'    : r'$max|\Delta x_i \Delta p_i|$',
                       }

################################################################################################
# user parameters
################################################################################################

def dict_update(d, u, invert=False):
    # retuns copy of updated dictionary d with u 
    # if invert==True inverts all lists in dict
    # if value in u is None key is removed from d. 
    # note: dict.update modifies dict and returns None
    c = {}
    for key,value in d.items():
        if isinstance(value, (list, tuple, np.ndarray)):
            c[key] = value[::-1] if invert else value
        else:
            c[key] = value
    for key,value in u.items():
        if value is None:
            try:
                del c[key]
            except KeyError:
                pass
        else:
            c[key] = value
    return c

# select figure to generate
figure = ['Greene2a', 'Greene2b', 'Greene2c', 'Greene2c-inv', 
          'timing', 'threads', 
          'LiCr_v1.2', 'test_data', 
          None][-1]

if figure == 'Greene2a':
    # this takes about 10' per molecule on my laptop
    
    # plot title
    title = 'K40K40 molecule conversion efficiency vs. distance measures (Greene Fig. 2a)'

    # atom number and trapping frequencies
    N       = [30000, 30000]    # atom number per species
    f_rad   = [470, 470]        # radial (x) trap frequency in Hz
    f_vert  = [470, 470]        # vertical (y) trap frequencies in Hz
    f_ax    = [6.7, 6.7]        # axial (z) trap frequency in Hz

    # variation
    # vary T/TF of first species, T[1] = T[0]
    vary        = 'T0=T1'
    vary_values = np.linspace(0.1, 1.6, 8) # scaling of T/TF 

    # random number and distribution generators
    rng = 'Lehmer128'
    dng = 'Metropolis'
    
    # other settings
    threads = 8
    reps    = 5

    # dict of molecule settings.
    # first species is varied, second is fixed or varied accordingly (see vary options).
    molecules = {
        'K40K40a': dict_update(K40K40, {
            'label'             : 'K40K40',        
            'rng'               : rng,
            'dng'               : dng,
            'distance measure'  : 'delta_x*delta_p',
            'gamma'             : 0.19*2,
            'data_args'         : [{'color': 'Red'}],
            'vary'              : vary,
            'vary_values'       : vary_values,
            'N'                 : N,                            # atom number per species
            'f_rad'             : f_rad,                        # radial (x) trap frequency in Hz
            'f_vert'            : f_vert,                       # vertical (y) trap frequencies in Hz
            'f_ax'              : f_ax,                         # axial (z) trap frequency in Hz
            'threads'           : threads,                      # number of threads
            'repetitions'       : reps,                         # repetitions per variation
            'calc_size'         : [1.0,1.0],                    # calculation size scaling
        }),
    }
    if False: 
        molecules = {
        'K40K39a': dict_update(K40K39, {
            'label'             : 'K40K39-test',
            'rng'               : rng,
            'dng'               : dng,
            'distance measure'  : 'delta_x*delta_p',
            'gamma'             : 0.38,
            'data_args'         : [{'color': 'Blue'}],
            'vary'              : vary,
            'vary_values'       : vary_values,
            'N'                 : N,                            # atom number per species
            'f_rad'             : f_rad,                        # radial (x) trap frequency in Hz
            'f_vert'            : f_vert,                       # vertical (y) trap frequencies in Hz
            'f_ax'              : f_ax,                         # axial (z) trap frequency in Hz
            'threads'           : threads,                      # number of threads
            'repetitions'       : reps,                         # repetitions per variation
            'calc_size'         : [1.0,1.0],                    # calculation size scaling
        }),
        'K39K40a': dict_update(K40K39, {
            'label'             : 'K39K40',
            'rng'               : rng,
            'dng'               : dng,
            'distance measure'  : 'delta_x*delta_p',
            'gamma'             : 0.38,
            'data_args'         : [{'color': 'Violet'}],
            'vary'              : vary,
            'vary_values'       : vary_values,
            'N'                 : N,                            # atom number per species
            'f_rad'             : f_rad,                        # radial (x) trap frequency in Hz
            'f_vert'            : f_vert,                       # vertical (y) trap frequencies in Hz
            'f_ax'              : f_ax,                         # axial (z) trap frequency in Hz
            'threads'           : threads,                      # number of threads
            'repetitions'       : reps,                         # repetitions per variation
            'calc_size'         : [1.0,1.0],                    # calculation size scaling
        }, invert=True),
        'K39K39a': dict_update(K39K39, {
            'rng'               : rng,
            'dng'               : dng,
            'distance measure'  : 'delta_x*delta_p',
            'gamma'             : 0.19*2,
            'data_args'         : [{'color': 'Green'}],
            'vary'              : vary,
            'vary_values'       : vary_values,
            'N'                 : N,                            # atom number per species
            'f_rad'             : f_rad,                        # radial (x) trap frequency in Hz
            'f_vert'            : f_vert,                       # vertical (y) trap frequencies in Hz
            'f_ax'              : f_ax,                         # axial (z) trap frequency in Hz
            'threads'           : threads,                      # number of threads
            'repetitions'       : reps,                         # repetitions per variation
            'calc_size'         : [1.0,1.0],                    # calculation size scaling
        }),
        'K39K41a': dict_update(K39K41, {
            'rng'               : rng,
            'dng'               : dng,
            'distance measure'  : 'delta_x*delta_p',
            'gamma'             : 0.19*2,
            'data_args'         : [{'color': 'Orange'}],
            'vary'              : vary,
            'vary_values'       : vary_values,
            'N'                 : N,                            # atom number per species
            'f_rad'             : f_rad,                        # radial (x) trap frequency in Hz
            'f_vert'            : f_vert,                       # vertical (y) trap frequencies in Hz
            'f_ax'              : f_ax,                         # axial (z) trap frequency in Hz
            'threads'           : threads,                      # number of threads
            'repetitions'       : reps,                         # repetitions per variation
            'calc_size'         : [1.0,1.0],                    # calculation size scaling
        }),
    }
    if False: 
        molecules = {
        'K40K40b': dict_update(K40K40, {
            'distance measure'  : 'cross(delta_x^,delta_p^)',
            'gamma'             : 0.26*3,
            'color'             : 'Green',
            'N'                 : N,                            # atom number per species
            #'T'                 : T,                            # temperature in nK (not used)
            'f_rad'             : f_rad,                        # radial (x) trap frequency in Hz
            'f_vert'            : f_vert,                       # vertical (y) trap frequencies in Hz
            'f_ax'              : f_ax,                         # axial (z) trap frequency in Hz
            'vary'              : vary,
            'vary_values'       : vary_values,    
            'vary_fixed'        : vary_fixed,
        }), 
        'K40K40c': dict_update(K40K40, {
            'distance measure'  : '<delta_x*delta_p>',
            'gamma'             : 0.0085*2,
            'color'             : 'Blue',
            'N'                 : N,                            # atom number per species
            #'T'                 : T,                            # temperature in nK (not used)
            'f_rad'             : f_rad,                        # radial (x) trap frequency in Hz
            'f_vert'            : f_vert,                       # vertical (y) trap frequencies in Hz
            'f_ax'              : f_ax,                         # axial (z) trap frequency in Hz
            'vary'              : vary,
            'vary_values'       : vary_values,
            'vary_fixed'        : vary_fixed,
        }), 
        'K40K40d': dict_update(K40K40, {
            'distance measure'  : 'max|delta_xi*delta_pi|',
            'gamma'             : 0.050*2,
            'color'             : 'Violet',
            'N'                 : N,                            # atom number per species
            #'T'                 : T,                            # temperature in nK (not used)
            'f_rad'             : f_rad,                        # radial (x) trap frequency in Hz
            'f_vert'            : f_vert,                       # vertical (y) trap frequencies in Hz
            'f_ax'              : f_ax,                         # axial (z) trap frequency in Hz
            'vary'              : vary,
            'vary_values'       : vary_values,
            'vary_fixed'        : vary_fixed,
        }), 
        'K40K40(MB)': dict_update(K40K40, {
            'label'             : 'K40K40(MB)',
            'atom_stat'         : ['MaxwellBoltzmann','MaxwellBoltzmann'],
            'mol_stat'          : 'MaxwellBoltzmann',
            'distance measure'  : 'delta_x*delta_p',
            'gamma'             : 0.19*2,
            'color'             : 'Black',
            'N'                 : N,                            # atom number per species
            'T'                 : [310, 310],                   # use K40 TF for scaling
            'f_rad'             : f_rad,                        # radial (x) trap frequency in Hz
            'f_vert'            : f_vert,                       # vertical (y) trap frequencies in Hz
            'f_ax'              : f_ax,                         # axial (z) trap frequency in Hz
            'vary'              : vary,
            'vary_values'       : vary_values,
            'vary_fixed'        : vary_fixed,
        }), 
        }

elif figure == 'Greene2b':

    # plot title
    title = 'K39K39 molecule conversion efficiency vs. distance measures (Greene Fig. 2b)'

    # variation
    # vary T/Tc of first species, T[1] = T[0]
    vary        = 'T0=T1'
    vary_values = np.linspace(0.1, 1.6, 16) # scaling of T/Tc 
    vary_fixed  = None # not used

    # dict of molecule settings.
    # first species is varied, second is fixed or varied accordingly (see vary options).
    molecules = {
        'K39K39a': dict_update(K39K39, {
            'distance measure'  : 'delta_x*delta_p',
            'gamma'             : 0.19*2,
            'color'             : 'Red',
            'vary'              : vary,
            'vary_values'       : vary_values,
            'vary_fixed'        : vary_fixed,
        }), 
        'K39K39b': dict_update(K39K39, {
            'distance measure'  : 'cross(delta_x^,delta_p^)',
            'gamma'             : 0.26*3,
            'color'             : 'Green',
            'vary'              : vary,
            'vary_values'       : vary_values,
            'vary_fixed'        : vary_fixed,
        }), 
        'K39K39c': dict_update(K39K39, {
            'distance measure'  : '<delta_x*delta_p>',
            'gamma'             : 0.0085*2,
            'color'             : 'Blue',
            'vary'              : vary,
            'vary_values'       : vary_values,
            'vary_fixed'        : vary_fixed,
        }), 
        'K39K39d': dict_update(K39K39, {
            'distance measure'  : 'max|delta_xi*delta_pi|',
            'gamma'             : 0.050*2,
            'color'             : 'Violet',
            'vary'              : vary,
            'vary_values'       : vary_values,
            'vary_fixed'        : vary_fixed,
        }), 
        'K39K39(MB)': dict_update(K39K39, {
            'label'             : 'K39K39(MB)',
            'T'                 : [160, 160], # use K39 Tc for scaling
            'atom_stat'         : ['MaxwellBoltzmann','MaxwellBoltzmann'],
            'mol_stat'          : 'MaxwellBoltzmann',
            'distance measure'  : 'delta_x*delta_p',
            'gamma'             : 0.19*2,
            'color'             : 'Black',
            'vary'              : vary,
            'vary_values'       : vary_values,
            'vary_fixed'        : vary_fixed,
        }), 
    }
elif figure == 'Greene2c':
    # plot title
    title = 'K39K40 molecule conversion efficiency vs. distance measures (Greene Fig. 2c)'

    # variation
    # vary T/Tc of first species, T[1] = T[0]
    vary        = 'T0=T1'
    vary_values = np.linspace(0.05, 1.6, 32) # scaling of T/Tc 
    vary_fixed  = None # not used

    # dict of molecule settings.
    # first species is varied, second is fixed or varied accordingly (see vary options).
    molecules = {
        'K39K40a': dict_update(K40K39, {
            'distance measure'  : 'delta_x*delta_p',
            'gamma'             : 0.19*2,
            'color'             : 'Red',
            'vary'              : vary,
            'vary_values'       : vary_values,
            'vary_fixed'        : vary_fixed,
        }, invert=True), 
        'K39K40b': dict_update(K40K39, {
            'distance measure'  : 'cross(delta_x^,delta_p^)',
            'gamma'             : 0.26*3,
            'color'             : 'Green',
            'vary'              : vary,
            'vary_values'       : vary_values,
            'vary_fixed'        : vary_fixed,
        }, invert=True), 
        'K39K40c': dict_update(K40K39, {
            'distance measure'  : '<delta_x*delta_p>',
            'gamma'             : 0.0085*2,
            'color'             : 'Blue',
            'vary'              : vary,
            'vary_values'       : vary_values,
            'vary_fixed'        : vary_fixed,
        }, invert=True), 
        'K39K40d': dict_update(K40K39, {
            'distance measure'  : 'max|delta_xi*delta_pi|',
            'gamma'             : 0.050*2,
            'color'             : 'Violet',
            'vary'              : vary,
            'vary_values'       : vary_values,
            'vary_fixed'        : vary_fixed,
        }, invert=True), 
        'K39K40(MB)': dict_update(K40K39, {
            'label'             : 'K39K40(MB)',
            'T'                 : [160, 160], # use K39 Tc for scaling
            'atom_stat'         : ['MaxwellBoltzmann','MaxwellBoltzmann'],
            'mol_stat'          : 'MaxwellBoltzmann',
            'distance measure'  : 'delta_x*delta_p',
            'gamma'             : 0.19*2,
            'color'             : 'Black',
            'vary'              : vary,
            'vary_values'       : vary_values,
            'vary_fixed'        : vary_fixed,
        }, invert=True), 
    }

elif figure == 'Greene2c-inv':
    # plot title
    title = 'K40K39 molecule conversion efficiency vs. distance measures (Greene Fig. 2c, vary K40 instead of K39)'

    # variation
    # vary T/Tc of first species, T[1] = T[0]
    vary        = 'T0=T1'
    vary_values = np.linspace(0.1, 1.6, 16) # scaling of T/Tc 
    vary_fixed  = None # not used

    # dict of molecule settings.
    # first species is varied, second is fixed or varied accordingly (see vary options).
    molecules = {
        'K40K39a': dict_update(K40K39, {
            'distance measure'  : 'delta_x*delta_p',
            'gamma'             : 0.19*2,
            'color'             : 'Red',
            'vary'              : vary,
            'vary_values'       : vary_values,
            'vary_fixed'        : vary_fixed,
        }), 
        'K40K39b': dict_update(K40K39, {
            'distance measure'  : 'cross(delta_x^,delta_p^)',
            'gamma'             : 0.26*3,
            'color'             : 'Green',
            'vary'              : vary,
            'vary_values'       : vary_values,
            'vary_fixed'        : vary_fixed,
        }), 
        'K40K39c': dict_update(K40K39, {
            'distance measure'  : '<delta_x*delta_p>',
            'gamma'             : 0.0085*2,
            'color'             : 'Blue',
            'vary'              : vary,
            'vary_values'       : vary_values,
            'vary_fixed'        : vary_fixed,
        }), 
        'K40K39d': dict_update(K40K39, {
            'distance measure'  : 'max|delta_xi*delta_pi|',
            'gamma'             : 0.050*2,
            'color'             : 'Violet',
            'vary'              : vary,
            'vary_values'       : vary_values,
            'vary_fixed'        : vary_fixed,
        }), 
        'K40K39(MB)': dict_update(K40K39, {
            'label'             : 'K40K39(MB)',
            'T'                 : [310, 310], # use K40 TF for scaling
            'atom_stat'         : ['MaxwellBoltzmann','MaxwellBoltzmann'],
            'mol_stat'          : 'MaxwellBoltzmann',
            'distance measure'  : 'delta_x*delta_p',
            'gamma'             : 0.19*2,
            'color'             : 'Black',
            'vary'              : vary,
            'vary_values'       : vary_values,
            'vary_fixed'        : vary_fixed,
        }), 
    }

elif figure == 'timing':
    # vary N of both species in steps of 2 for fixed Cr52 T/Tc and TLi = TCr and measure calculation time
    # note: on my laptop (Lenovo Thinkpad E14) this takes 13-14' (9') per test on Ubuntu 22.04 LTS (Windows 10)
    #       on Linux there is only a small difference between Lehmer64 and Lehmer128, on Windows its larger but still faster.
    #       the actual time depends on the PC, performance settings, and if the laptop is on Battery or charging, 
    #       and how the OS performs CPU throttling (see also 'threads' test below) and thermal management. 
    #       A difference of 50% even on similar hardware is not unusual.
    title = 'timing test'
    num = 10 # number of variations
    Nratio = 0.5 # Nc/N. larger values give longer calculation time since overlap between BEC and FG is smaller than MB and FG
    test0 = {
        'label'             : 'Cr52Li6_timing_Lehmer64',
        'distance measure'  : 'delta_x*delta_p',
        'gamma'             : 0.19*2,
        'rng'               : 'Lehmer64',
        'color'             : 'Gray',
        'data_args'         : [{'color':'Orange'},{'color':'Green'},{'color':'Blue'},{'color':'Red'}],
        'vary'              : 'N0&N1',
        'N'                 : [128000, 128000], # maximum N of both species
        'vary_values'       : np.array([2**(i-num+1) for i in range(num)]), # scaling of N of both species
        'vary_fixed'        : (1-Nratio)**(1/3), # fixed T/TF or T/Tc of first species. second species T1=T0 
        'repetitions'       : 5,    # repetitions per variation
    }
    test1 = test0.copy()
    test1.update({
        'label'             : 'Cr52Li6_timing_Lehmer128',
        'color'             : 'Black',
        'data_args'         : [{'color':'Orange', 'edgecolor':'Orange', 'facecolor':'White'},
                               {'color':'Green' , 'edgecolor':'Green' , 'facecolor':'White'},
                               {'color':'Blue'  , 'edgecolor':'Blue'  , 'facecolor':'White'},
                               {'color':'Red'   , 'edgecolor':'Red'   , 'facecolor':'White'}],
        'rng'               : 'Lehmer128'}
    )
    molecules = {
        'Cr52Li6_timing_0': dict_update(Li6Cr52, test0, invert=True),
        'Cr52Li6_timing_1': dict_update(Li6Cr52, test1, invert=True),
    }
elif figure == 'threads':
    # perform the same calculation with different number of threads
    # this takes about 45' on my laptop with N=75k
    # note: this checks how much one gains by increasing number of threads and how much overhead this introduces.
    #       for independent threads the calculation time should scale by 1/number of threads.
    #       but for the calculation the threads must access shared data which requires synchronization of data access
    #       and causes temporary blocking of threads and additional overhead leading to a scaling worse than 1/number of threads.
    #       additionally, the CPU reduces the clock frequency when more cores are active which makes scaling even worse.
    #       this depends on CPU performance setting thermal management of the hardware and the OS. 
    title = 'thread performance test'
    num = 8 # maximum number of threads
    Nratio = 0.5 # Nc/N. larger values give longer calculation time since overlap between BEC and FG is smaller than MB and FG
    molecules = {
        'Cr52Li6_%i_threads'%(1+i) : dict_update(Li6Cr52, {
            'label'             : 'Cr52Li6_%i_threads'%(1+i),
            'distance measure'  : 'delta_x*delta_p',
            'gamma'             : 0.19*2,
            'repetitions'       : 5,    # repetitions per variation
            'vary'              : 'N0&N1',
            'N'                 : [50000, 50000], # maximum N of both species
            'vary_values'       : np.array([1]), # scaling of N of both species
            'vary_fixed'        : (1-Nratio)**(1/3), # fixed T/TF or T/Tc of first species. second species T1=T0 
            'color'             : 'Blue',
            'data_args'         : [{'color':'Orange'},{'color':'Green'},{'color':'Blue'},{'color':'Red'}],
            'threads'           : 1+i,
        }, invert=True)
    for i in range(num) }
elif figure == 'LiCr_v1.2':
    # Li6Cr53 test case for comparison with v1.2 (22/4/2024)
    title = 'LiCr_v1.2_test_case'
    molecules = {
        'Li6Cr53_v1.2' : dict_update(Li6Cr53, {
            'label'             : 'Li6Cr53',
            'atom_stat'         : ['FermiDirac','FermiDirac'], # statistics of both species            
            'mol_stat'          : 'BoseEinstein',              # statistics of molecule
            'distance measure'  : 'delta_x*delta_p',           # former 'momentum-space'
            'rng'               : 'Lehmer64',
            'gamma'             : 0.19*2,
            'vary'              : 'T0=T1', # vary T/TF of Li6, TCr = TLi
            'N'                 : [100000, 100000], # N of both species
            'T'                 : [6.322, 6.322],   # T of both species
            'f_rad'             : [100, 33.646329], # radial   trapping frequency in Hz
            'f_vert'            : [100, 33.646329], # vertical trapping frequency in Hz
            'f_ax'              : [100, 33.646329], # axial    trapping frequency in Hz
            'vary_values'       : [6.322/404.784],  # T/TF of Li6
            'vary_fixed'        : 0.1,              # not used
            'color'             : 'Blue',
            'threads'           : 1,
            'repetitions'       : 1,                # repetitions per variation
            'calc_size'         : [1.5,2.5],
        })
    }
elif figure == 'test_data':
    # comparison with test data 2023/12/11 from Alessio
    # see Results.dat and ExpConfig.dat in ./LiCr/Alessio_PSD/20240222
    # trapping frequencies and TLi are in TrapFreq.csv and Temperature.csv in ./LiCr/Alessio_PSD/20231206
    # Tvar[0][0] is not given and adjusted manually
    # MolPSD_8 does not exist in Temperature.csv or is mis-labelled as 2nd MolPSD_2 (assume this is the case)
    # there are some inconsistencies, take the ones from Results.dat:
    # - MolPSD_2 NCr = 40k vs. 20k 
    # - MolPSD_5 NCr = 36k vs. 40k
    # - MolPSD_7 NCr = 38k vs. 40k
    mol = Li6Cr53
    TTcrit = [None]*2
    Tvar   = [None]*2
    Nvar   = [None]*2
    f_rad  = [None]*2
    f_ax   = [None]*2
    f_vert = [None]*2
    td_label   =          ['MolPSD_3','MolPSD_1','MolPSD_0','MolPSD_2','MolPSD_8','MolPSD_4','MolPSD_5','MolPSD_6','MolPSD_7']
    TTcrit[0]  = np.array([ 0.2      ,   0.35   ,   0.35   ,   0.12   ,   0.20   ,   0.12   ,   0.12   ,   0.12   ,   0.12   ])
    TTcrit[1]  = np.array([ 1.0      ,   1.0    ,   2.0    ,   0.50   ,   0.45   ,   0.5    ,   0.50   ,   0.50   ,   0.50   ])
    Tvar  [0]  = np.array([  60      ,   120    ,   166    ,   30     ,   52     ,   30     ,   30     ,   30     ,   30     ])
    Tvar  [1]  = np.array([ 680      ,   270    ,   560    ,   130    ,   135    ,   130    ,   130    ,   130    ,   130    ])
    Nvar  [0]  = np.array([ 450e3    ,   120e3  ,   100e3  ,   60e3   ,   80e3   ,   60e3   ,   60e3   ,   60e3   ,   60e3   ])
    Nvar  [1]  = np.array([  90e3    ,    90e3  ,   130e3  ,   40e3   ,   40e3   ,   40e3   ,   36e3   ,   40e3   ,   38e3   ])
    f_rad [0]  = np.array([ 300.0    ,   340.0  ,   340.0  ,   304.0  ,   755.0  ,   304.0  ,   304.0  ,   304.0  ,   304.0  ])
    f_vert[0]  = f_rad[0].copy() # ensure to copy, otherwise scaling below is wrong!
    f_ax  [0]  = np.array([ 17.2     ,   17.2   ,   17.2   ,   17.1   ,   19.3   ,   17.1   ,   17.1   ,   17.1   ,   17.1   ])
    f_rad [1]  = np.array([ 182.0    ,   132.0  ,   132.0  ,   97.0   ,   164.0  ,   97.0   ,   97.0   ,   97.0   ,   97.0   ])
    f_vert[1]  = f_rad[1].copy() # ensure to copy, otherwise scaling below is wrong!
    f_ax  [1]  = np.array([ 14.0     ,   14.0   ,   14.0   ,   14.0   ,   14.2   ,   14.0   ,   14.0   ,   14.0   ,   14.0   ])

    # expected results
    td_N_mol   = np.array([ 12e3     ,   17.2e3 ,   10e3   ,   7.6e3  ,   7.6e3  ,   19.76e3,   15.86e3,   24.96e3,   27.04e3])
    td_f_mol   = np.array([ 0.1333   ,   0.1911 ,   0.0769 ,   0.19   ,   0.19   ,   0.494  ,   0.4406 ,   0.624  ,   0.71158])
    td_T_mol   = np.array([ 308      ,   228    ,   388    ,   87     ,   87     ,   140    ,   115    ,   110    ,   120    ])
    td_PSD_mol = np.array([ 0.012    ,   0.062  ,   0.0062 ,   0.11   ,   0.11   ,   0.0936 ,   0.1144 ,   0.1131 ,   0.0988 ])

    # check conversion efficiency
    # it seems the td_f_mol in list was calculated relative to NCr and not vs. min(NLi, NCr)
    max_error = 1e-3
    error = np.abs(td_f_mol-td_N_mol/np.min(Nvar, axis=0))
    count = np.count_nonzero(error >= max_error)
    if count > 0:
        print('warning: conversion efficiency inconsistent for', count, 'entries!')
        for i in range(len(error)):
            if error[i] >= max_error:
                print('\t', td_label[i], 'eta = ', td_N_mol[i], ' / ', min(Nvar[0][i], Nvar[1][1]), ' = %.3f' % (td_N_mol[i]/min(Nvar[0][i], Nvar[1][1])), '!= %.3f' % td_f_mol[i], '(error %.1e)' % error[i])
        #exit()

    if False:
        # make spheric trap for comparison
        f_ax[0] = f_rad[0].copy()
        f_ax[1] = f_rad[1].copy()
        label = 'MolPSD_spheric'
    else:
        # use original cigar-shaped trap
        label = 'MolPSD_cigar'

    # Li and Cr atom number does not match T and T/TF coming from fit. 
    # rescale the frequencies to make these matching
    rescale = [True, True]
    for i in range(2):
        if rescale[i]:
            TF = Tvar[i]/TTcrit[i]
            omega_bar = 2*np.pi*(f_rad[i]*f_ax[i]*f_vert[i])**(1/3)
            sc = TF / (omega_bar*((6.0*Nvar[i])**(1.0/3.0))*k_nK)
            print('note: scale', mol['species'][i], 'trap frequencies to match T/TF of fit!')
            print('scaling factor:', sc)
            f_rad [i] *= sc
            f_ax  [i] *= sc
            f_vert[i] *= sc

    title = 'LiCr 2023/12/11 test data'
    molecules = {
        'MolPSD' : dict_update(mol, {
            'label'             : label,
            'atom_stat'         : ['FermiDirac','FermiDirac'],  # statistics of both species            
            'mol_stat'          : 'BoseEinstein',               # statistics of molecule
            'distance measure'  : 'delta_x*delta_p',            # former 'momentum-space'
            'rng'               : 'Lehmer128',
            'gamma'             : 0.38,                          # increased to match data
            'vary'              : 'list',                       # direct input of data
            'color'             : 'Blue',
            'data_args'         : [{'color':'Blue'},{'color':'Orange'},{'color':'Green'}], # simulated, measured Nmol/max_pairs, measured data
            'N'                 : None,                         # clear N of both species
            'T'                 : None,                         # clear T of both species
            'f_rad'             : None,                         # clear radial   trapping frequency in Hz
            'f_vert'            : None,                         # clear vertical trapping frequency in Hz
            'f_ax'              : None,                         # clear axial    trapping frequency in Hz
            'threads'           : 8,                            # number of threads
            'repetitions'       : 5,                            # repetitions per variation
            'calc_size'         : [2.0,2.0],
        })
    }
    
elif False:
    # custom data set

    title = 'test'

    vary        = 'T0=T1'
    vary_values = np.linspace(0.1, 1.6, 4) # scaling of T/TF 
    vary_fixed  = None # not used

    U_ratio = 0.5
    molecules = {
        'test': dict_update(Li6Cr52, {
            'label'             : 'test',
            'atom_stat'         : ['BoseEinstein','MaxwellBoltzmann'],  # statistics of both species            
            'N'                 : [50000, 50000],
            'f_rad'             : [100, 100*np.sqrt(U_ratio*6/53)], # radial (x) trap frequency in Hz
            'f_vert'            : [ 50,  50*np.sqrt(U_ratio*6/53)], # vertical (y) trap frequencies in Hz
            'f_ax'              : [ 15,  10],   # axial (z) trap frequency in Hz
            'rng'               : 'Lehmer128',
            'dng'               : 'Metropolis',
            'distance measure'  : 'delta_x*delta_p',
            'gamma'             : 0.19*2,
            'data_args'         : [{'color':'Red'}],
            'vary'              : 'T0=T1', # vary-values = Li T/TF, Cr T = Li T
            'vary_values'       : np.linspace(0.7, 1.0, 1), # Li T/TF
            'threads'           : 8,                            # number of threads
            'repetitions'       : 5,                            # repetitions per variation
            'calc_size'         : [1.0,1.0],                    # calculation size/energy scaling
        }, invert=True), 
    }

else:
    # custom data set

    title = 'test'

    molecules = {
        'K40K40': dict_update(K40K40, {
            'label'             : 'K40K40',
            'rng'               : 'Lehmer128',
            'dng'               : 'Metropolis',
            'distance measure'  : 'delta_x*delta_p',
            'gamma'             : 0.38,
            'data_args'         : [{'color': 'Red'}],
            'N'                 : [30000, 30000],   # atom number per species
            'T'                 : [100  , 100],     # temperature in nK (might be overwritten)
            'f_rad'             : [470, 470],       # radial (x) trap frequency in Hz
            'f_vert'            : [470, 470],       # vertical (y) trap frequencies in Hz
            'f_ax'              : [6.7, 6.7],       # axial (z) trap frequency in Hz
            'vary'              : 'T0=T1',
            'vary_values'       : np.linspace(0.1, 1.6, 16), # scaling of T/TF
            'vary_fixed'        : None,
            'threads'           : 8,                            # number of threads
            'repetitions'       : 5,                            # repetitions per variation
            'calc_size'         : [1.0,1.0],                    # calculation size/energy scaling
        }), 
    }

# if True recalculate existing results, otherwise just plot results
recalc = True

# general output folder
# note: relative path and '/' works also on Windows. 
folder = './tmp/'
#folder = './test_cases/20260827_v1.6/Greene/'

# if not None generate histogram and atoms/molecule files with this filename in same folder as result file.
# attention: files might be large and generation might take some time!
histogram_file  = [None, '_hist.dat'    ][0]
export_atom_0   = [None, '_atom0.csv'   ][0]
export_atom_1   = [None, '_atom1.scv'   ][0]
export_molecule = [None, '_molecule.csv'][0]

# number of bins for histogram
num_bins = 300

# parameters for calculation of chemical potential
error_allowed = 1e-11
step_size     = 0.1

# offset of 2nd species
offset = [0,0,0]

# if 0 stops at first matching pair (default, faster), otherwise searches nearest pair (slower)
find_nearest = 0

# show T/Tc or T/TF of second species as twin axis
show_second_species = True

# acceptable timing error in s
t_err = 1e-3

# seed values if not None. (must be within {} and ',' as separator. no spaces allowed. several values allowed, >2 not really needed )
seed = [None, "{0x53b2c9e8,0x869c86ac}"][0]

################################################################################################
# results to be collected from result_file
################################################################################################

# note that parameters might appear several times per iteration in result file.
# to distinguish we need additional parameters "atomic species" and "molecule"
# on the right side give format + units within (). units must match result file,
# format must be "%s", "%f" or "%i"
# each parameter must be at beginning of line and must not have additional characters before '='
# note: for similar starting parameters give longest first to avoid warning messages

results = {}
results['repetition       [summary]']   = ["%i/"] # TODO: read repetitions
results["repetition"]                   = ["%i/"] # TODO: read repetitions
results["atomic species"]               = ["%s"]
results["N"]                            = ["%i"]
results["T"]                            = ["%f(nK)"]
results["TF"]                           = ["%f(nK)"]
results["T/TF"]                         = ["%f"]
results["Tc"]                           = ["%f(nK)"]
results["T/Tc"]                         = ["%f"]
results['conv. efficiency [summary]']   = ["%f+/-(%)"] # TODO: read error
results["conv. efficiency"]             = ["%f(%)"]
results["molecule"]                     = ["%s"]
results["N_mol"]                        = ["%i"]
results["chemical potential"]           = ["%f(nK)"]
results["calculation time [summary]"]   = ["%f+/-(s)"]
results["calculation time [total]"]     = ["%f(s)"]
results["calculation time"]             = ["%f(s)"]

################################################################################################
# functions
################################################################################################

def plot(data, data_labels, curves, curve_labels, data_args=None, curve_args=None, title=None, 
            fig_size=(7*4/3,7), fig_pos=[0.08, 0.07, 0.91, 0.87],
            xticks=None, xlabel=None, ylabel=None, 
            x_range=None, y_range=None, log_scale=None,
            label_pos='upper right', label_cols=1,
            bitmap=None, bm_args={},
            twin=None):
    """
    plot data and curves with given labels.
    data         = list of [[x0,y0],[x1,y1],...] plotted as data points
    data_labels  = labels of each data
    data_args    = dictionary with plot arguments like color, style, etc. if None, default values are taken.
    curves       = list of [[x0,y0],[x1,y1],...] plotted as lines between the points
    curve_labels = labels of each curve
    curve_args   = dictionary with plot arguments like color, style, etc. if None, default values are taken.
    title        = if not None title of figure
    xlabel       = if not None x label of figure
    ylabel       = if not None y label of figure
    xticks       = list of [[ticks],[labels]] for x-axis
    bitmap       = 2d bitmap of pixels
    bm_args      = dictionary with bitmap arguments
    twin         = if not None give axis returned by previoys plot and all data is plotted on twin axis
    """
    if twin is not None:
        ax = twin.twinx()
    else:
        fig = plt.figure(figsize=fig_size) #size in inches (width,height)
        ax = fig.add_axes(fig_pos)

    if title is not None: ax.set_title(title, fontsize=14)
    if xlabel is not None: ax.set_xlabel(xlabel, fontsize=14, labelpad=5)
    if ylabel is not None: ax.set_ylabel(ylabel, fontsize=14, labelpad=5)

    ax.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', top  ='on', labelsize=12)
    ax.xaxis.set_tick_params(which='minor', size= 7, width=2, direction='in', top  ='on', labelsize=12)
    ax.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', right='on', labelsize=12)
    ax.yaxis.set_tick_params(which='minor', size= 7, width=2, direction='in', right='on', labelsize=12)
    
    if x_range is not None:
        if len(x_range) == 2:
            ax.set_xlim(x_range[0], x_range[1])
        elif len(x_range) == 3:
            ax.set_xlim(x_range[0], x_range[1])
            if x_range[2] is not None:
                ax.set_xticks(x_range[2])
        else:
            print('x_range must be [min,max] or [[min,max],[ticks]]')
            exit()
            
    if y_range is not None:
        if len(y_range) == 2:
            ax.set_ylim(y_range[0], y_range[1])
        elif len(y_range) == 3:
            ax.set_ylim(y_range[0], y_range[1])
            if y_range[2] is not None:
                ax.set_yticks(y_range[2])
        else:
            print('y_range must be [min,max] or [[min,max],[ticks]]')
            exit()

    if log_scale is not None:
        if log_scale[0]: ax.set_xscale('log')
        if log_scale[1]: ax.set_yscale('log')
        
    if xticks is not None:
        #ax.set_xticks(xticks[0], xticks[1]) # does not work!?
        ax.set_xticks(xticks[0])
        ax.set_xticklabels(xticks[1])
        ax.tick_params(axis='x', rotation=70)

    if False:
        scale = 0.05
        y_min = np.min([np.min(d[0]) for d in data])
        y_max = np.max([np.max(d[0]) for d in data])
        dy = (y_max-y_min)*scale
        ax.set_xlim(y_min-dy, y_max+dy)

        y_min = np.min([np.min(d[1]) for d in data])
        y_max = np.max([np.max(d[1]) for d in data])
        dy = (y_max-y_min)*scale
        ax.set_ylim(y_min-dy, y_max+dy)

    #plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    #plt.yscale('log')

    if bitmap is not None:    
        plt.imshow(bitmap, **bm_args)

    for i,d in enumerate(data):
        if d is not None:
            if data_args is None or data_args[i] is None:
                args = {}
            else:
                args = data_args[i]
                #if not ('color' in args):
                #    args['color'] = colors(((i*2)+int(i//colors.N))%colors.N)
            if len(d) == 3: # plot error bar
                plot_error_bar_exclude = ['s', 'edgecolor', 'facecolor'] # remove options not good for error bars
                e_args = {k:v for k,v in args.items() if k not in plot_error_bar_exclude}
                if not 'linewidth' in e_args:
                    e_args['linewidth'] = 2
                if not 'zorder' in e_args:
                    e_args['zorder'] = 1
                for j in range(len(d[0])):
                    ax.plot([d[0][j]]*2, [d[1][j]-d[2][j], d[1][j]+d[2][j]], label=None, **e_args)
            # plot data point
            ax.scatter(d[0], d[1], label=data_labels[i], **args)

    for i,d in enumerate(curves):
        if d is not None:
            if (curve_args is None) or (curve_args[i] is None):
                args = {}
            else:
                args = curve_args[i].copy() # without copy can get unexpected results when lists were multiplied
            if not ('color' in args):
                args['color'] = colors(((i*2)+int(i//colors.N))%colors.N)
            ax.plot(d[0], d[1], label=curve_labels[i], **args)

    if label_pos is not None:
        label_xy = {'upper right' : (0.97, 0.97), 
                    'upper left'  : (0.03, 0.97),
                    'upper center': (0.50, 0.97),
                    'center right': (0.97, 0.50),  
                    'lower left'  : (0.03, 0.05), 
                    'lower right' : (0.97, 0.03),
                    'lower center': (0.50, 0.03)}[label_pos]
        ax.legend(bbox_to_anchor=label_xy, loc=label_pos, frameon=True, fontsize=10, ncol=label_cols, framealpha=0.5)
    
    #plt.savefig(title + '.png', dpi=300, transparent=False, bbox_inches='tight')
    return ax
    
def split_unit(text, sep=['(',')']):
    # returns [format,unit] from text=format+'('+unit+')'
    # if no unit given returns unit = ''
    text_ = text.split(sep[0])
    if len(text_) == 1: 
        return [text, '']
    if len(text_) > 2:
        text_ = [text_[0], sep[0].join(text_[1:])]
    if len(sep) == 1:
        return text_
    elif len(sep) == 2:
        _text = text_[1].split(sep[1])
        if len(_text) != 2: # no or several closing brackets!
            print('illegal format given:', text)
            exit()
        else:
            return [text_[0],_text[0]]

def read_number(text, fmt, force_unit):
    # read number as given format string
    # remove units from text and format assumed either with space of within round brackets
    # if force_unit = True and number given then text must contain unit given in fmt. if not returns None
    text_ = split_unit(text, sep=[' '])
    fmt_  = split_unit(fmt)
    if len(text_) != 2:
        return None
    if fmt_[0] == "%s":
        return text.strip()
    else:
        if fmt_[0] == "%i": 
            if force_unit and text_[1] != fmt_[1]:
                return None
            try:
                return int(text_[0])
            except ValueError:
                print("could not convert '%s' as integer!" % text)
                exit()
        elif fmt_[0] == "%i/": # %i / %i TODO: return also 2nd integer
            text_ = text.split(' ')
            if force_unit:
                if len(fmt_[1]) > 0 and (len(text_) != 4 or text_[1] != '/' or text_[3] != fmt_[1]):
                    print('format %i/ with unit error!', text_, fmt_)
                    return None
                elif len(fmt_[1]) == 0 and (len(text_) != 3 or text_[1] != '/'):
                    print('format %i/ without unit error!', text_, fmt_)
                    return None
            elif len(text_) != 4 and len(text_) != 3:
                print('format %i/ invalid (length)!', text_, fmt_)
                return None
            elif text[1] != '/':                
                print('format %i/ invalid (/)!', text_, fmt_)
                return None
            try:
                value = [int(text_[0]), int(text_[2])]
                return value[0]
            except ValueError:
                print("could not convert '%s' as integer!" % text)
                exit()
        elif fmt_[0] == "%f":
            if force_unit and text_[1] != fmt_[1]:
                return None
            try:
                return float(text_[0])
            except ValueError:
                print("could not convert '%s' as float!" % text)
                exit()
        elif fmt_[0] == "%f+/-": # %f +/- %f TODO: return also second number
            text_ = text.split(' ')
            if force_unit:
                if len(fmt_[1]) > 0 and (len(text_) != 4 or text_[1] != '+/-' or text_[3] != fmt_[1]):
                    print('format %f+/- with unit error!', text_, fmt_)
                    return None
                elif len(fmt_[1]) == 0 and (len(text_) != 3 or text_[1] != '+/-'):
                    print('format %f+/- without unit error!', text_, fmt_)
                    return None
            elif len(text_) != 4 and len(text_) != 3:
                print('format %f+/- invalid (length)!', text_, fmt_)
                return None
            elif text[1] != '+/-':                
                print('format %f+/- invalid (+/-)!', text_, fmt_)
                return None
            try:
                value = [float(text_[0]), float(text_[2])]
                return value[0]
            except ValueError:
                print("could not convert '%s' as float +/- float!" % text)
                exit()
        else:
            print("format '%s' unknown!" % (fmt))
            exit()

def read_file(filename, params, allow_multiple=False, show=False, force_unit=True):
    # read parameter or result file and replaces parameter_name with parameter_fmt
    # returns [replaced_text, replace_order, replace_params, number_lines]
    # replaced_text = non-commented valid text of file with fmt replaced for old params
    # replace_order = list with params keys as they appeared in file
    # replace_param = list with old params
    # lines         = number of lines in replaced_text
    # if allow_multiple = False: function will give an error if one parameter appears multiple times.
    # if allow_multiple = True:  function returns also multiple instances of the same paramamter.
    lines = 0
    replace_text  = "" # content of file
    replace_order = [] # list with order of parameters to be inserted
    replace_param = [] # list with old parameters
    with open(filename, 'r') as f:
        while(True):
            line = f.readline()
            if len(line) == 0: break
            line_str = line.strip()
            if len(line_str) > 0:
                if line_str[0:2] != '//':
                    lines += 1
 
                    # find parameter in line
                    # to be sure we find not a sub-string we add " " to finding parameter
                    # we also check then that nothing is between parameter and "="
                    found = False
                    for p,value in params.items():
                        start = line_str.find(p + " ")
                        if start >= 0: 
                            # parameter 'p' found
                            if start > 0:
                                if show:
                                    print("note: line %i: parameter '%s' found at nonzero offset! (ignore)" % (lines, p))
                                    print(line_str)
                                continue                            
                            if not allow_multiple and p in replace_order:
                                print("error line %i: parameter '%s' found twice!" % (lines, p))
                                exit()
                            # find '=' after parameter
                            stop = line_str.find("=", start + len(p))
                            if stop > 0:
                                if len(line_str[start+len(p):stop-1].strip()) != 0:
                                    # non white space between parameter p and '='
                                    print('%3i:'%lines, "extended parameter '%s' > '%s'! (ignore)" % (p, line_str[start:stop-1].strip()))
                                    continue
                                old = line_str[stop+1:].strip()
                                number = read_number(old, value[0], force_unit)
                                if number is None:
                                    # parameter with wrong unit
                                    txt = split_unit(old, sep=[' '])
                                    fmt = split_unit(value[0])
                                    print('%3i:'%lines, "parameter '%s' = '%s' with wrong unit '%s' != '%s'! (ignore)" % (p, old, txt[1], fmt[1]))
                                    continue             
                                if show: print('%3i:'%lines, line_str[start:stop+1], "'%s' -> '%s'" % (old, params[p][0]))
                                replace_order.append(p)
                                replace_param.append(number)
                                replace_text += line_str[start:stop+1] + " " + params[p][0] + "\n"
                                found = True
                                break
                            else:
                                print("error line %i: cannot find '=' in '%s'" % (lines, line_str))
                                print(line)
                                print(p)
                                print([start, len(p), stop])
                                exit()
                    if not found: # normal line
                        if show: print('%3i:'%lines, line_str)
                        replace_text += line
    return [replace_text, replace_order, replace_param, lines]    
        
def find_all(lst, value):
    "returns all indices of value in list lst. empty list if not found."
    index = []
    l = lst.copy()
    removed = 0
    while len(l) > 0:
        try:
            i = l.index(value)
        except ValueError as e:
            break
        index.append(removed + i)
        l = l[i+1:]
        removed += i+1 
    return index

def print_data(data, header, left=' ', ffmt=' %13s'):
    fmt = (left + (ffmt*len(header)))
    print(fmt % tuple(header))
    print(np.transpose(data))
        
def format_coord(M, x_values, label, x, y):
    xi = round(x)
    if xi < 0: xi = 0
    if xi > len(x_values[0]): xi = len(x_values[0])-1
    yi = round(y)
    if yi < 0: yi = 0
    if yi > len(x_values[1]): yi = len(x_values[1])-1
    zi  = M[yi,xi]
    return "%.1f, %.1f, %s = %.3f"%(x_values[0][xi],x_values[1][yi],label,zi)
    
################################################################################################
# calculation starts here
################################################################################################

if __name__ == '__main__':

    # loop over molecules
    data        = []
    data_labels = []
    data_args   = []
    label_x     = [None, None]
    twin_data   = []
    twin_labels = []
    twin_args   = []
    twin_y      = None

    for mol_key, mol_dict in molecules.items():

        # specific values must be defined in molecules:
        print(mol_key)
        print(mol_dict)
        species                 = mol_dict['species']
        mass                    = mol_dict['mass']
        statistics              = mol_dict['atom_stat'] + [mol_dict['mol_stat']]
        distance_measure        = mol_dict['distance measure']
        gamma                   = mol_dict['gamma']
        label                   = mol_dict['label']
        #color                   = mol_dict['color']
        data_args_mol           = mol_dict['data_args']
        f_rad                   = mol_dict['f_rad']
        f_vert                  = mol_dict['f_vert']
        f_ax                    = mol_dict['f_ax']
        vary                    = mol_dict['vary']
        vary_values             = mol_dict['vary_values']
        repetitions             = mol_dict['repetitions']
        random_number_generator = mol_dict['rng']
        distribution_generator  = mol_dict['dng']
        threads                 = mol_dict['threads']
        calc_size               = mol_dict['calc_size']

        # depending on variation these might not be needed. if used will give an error later down.
        try:
            N                   = mol_dict['N']
        except KeyError:
            N                   = None
        try:
            T                   = mol_dict['T']
        except KeyError:
            T                   = None
        try:            
            vary_fixed          = mol_dict['vary_fixed']
        except KeyError:
            vary_fixed          = None

        if label is None:
            # automatically generate a label for folder and plotting
            label = species[0]+species[1]
        
        # output folder
        mol_folder = folder + label + '_' + distance_measure + '_%.1e'%gamma+'/'
        if os.name == 'nt': 
            mol_folder = mol_folder.replace('*','_')

        # single species (True) or two species (False)
        single_species = (species[0] == species[1])
        
        # number of species 1 or 2
        num_species = 1 if single_species else 2
        
        # species names as expected in result_file
        # for Bosons thermal and condensed part are generated separately
        species_names = [[]]*num_species
        for i in range(num_species):
            if statistics[i].startswith('BoseEinstein'):
                species_names[i] = [species[i] + '(thermal)', species[i] + '(condensed)']
            else:
                species_names[i] = [species[i]]
        #species_names_flat = [nn for n in species_names for nn in n]

        # molecule names as expected in result_file
        # for Bosons thermal and condensed part are generated separately
        mol_names = []
        for i,sp0 in enumerate(species_names[0]):
            for j,sp1 in enumerate(species_names[num_species-1]):
                if single_species:
                    if i <= j: # Fermions=0+0, Bosons=0+0,0+1,1+1
                        if i == j:
                            mol_names.append([sp0+sp1+' (single species)', i, j])
                        else:
                            mol_names.append([sp0+sp1, i, j])
                else: # F+F=0+0, B+B=0+0,0+1,1+0,1+1, F+B=0+0,0+1, B+F=0+0,1+0
                    mol_names.append([sp0+sp1, i, j])
                
        # parameter file used as template
        default_parameter_file = mol_folder + label + "_params.txt"

        # parameter file name. give '%i' for each variation.
        parameter_file = mol_folder + label + "_params_%i.txt"

        # result file name. this contains all results generated by molConv for each variation
        result_file = mol_folder + label + "_result.dat"

        # result summary file (csv) for all variations
        csv_file = mol_folder + label + "_result.csv"

        # list of parameters for each variation 
        omega_x     = [np.zeros(shape=(1,), dtype=float) for _ in range(num_species)]
        omega_y     = [np.zeros(shape=(1,), dtype=float) for _ in range(num_species)]
        omega_z     = [np.zeros(shape=(1,), dtype=float) for _ in range(num_species)]
        omega_bar   = [np.zeros(shape=(1,), dtype=float) for _ in range(num_species)]

        if vary == 'list':
            # list of T/Tcrit, N, omega for each species and check == T given 
            Tcrit = [np.zeros(shape=(1,), dtype=float) for _ in range(num_species)]
            max_error = 1e-2
            for i in range(num_species):
                omega_x  [i] = 2*np.pi*f_rad [i]
                omega_y  [i] = 2*np.pi*f_vert[i]
                omega_z  [i] = 2*np.pi*f_ax  [i]
                omega_bar[i] = (omega_x[i]*omega_y[i]*omega_z[i])**(1.0/3.0)
                if statistics[i] == 'BoseEinstein':
                    Tcrit[i] = hbar * omega_bar[i] * ( Nvar[i] / zeta )**(1.0/3.0) *1e9/kB # in nK
                elif statistics[i] == 'FermiDirac':
                    Tcrit[i] = omega_bar[i]*(6*Nvar[i])**(1.0/3.0)*k_nK # nK
                else:
                    print('species', species_names[i], "statistics", statistics[i], "not implemented!")
                    exit()
                error = np.abs(Tvar[i] - TTcrit[i]*Tcrit[i])
                count = np.count_nonzero(error > max_error)
                print('species', species_names[i], 'vary=list T error = %.1e' % np.max(error), 'nK', '(ok)' if count == 0 else '(error)')
                if count > 0:
                    print('error for', count, 'entries: T error > %.1e' % max_error)
                    print('T              ', Tvar[i].astype(float))
                    print('T/Tcrit*Tcrit  ', TTcrit[i]*Tcrit[i])
                    print('delta T        ', error)
                    print('N              ', Nvar[i])
                    print('N(Tcrit, omega)', np.round(((Tcrit[i] / (omega_bar[i]*k_nK))**3)/6)) 
                    print('delta N        ', np.abs(Nvar[i] - np.round(((Tcrit[i] / (omega_bar[i]*k_nK))**3)/6)))
                    exit()
        else:
            Tvar        = [np.zeros(shape=(1,), dtype=float) for _ in range(num_species)]
            Nvar        = [np.zeros(shape=(1,), dtype=float) for _ in range(num_species)]
            Tcrit       = [np.zeros(shape=(1,), dtype=float) for _ in range(num_species)]
            TTcrit      = [np.zeros(shape=(1,), dtype=float) for _ in range(num_species)]

            for i in range(num_species):
                omega_x  [i][0] = 2*np.pi*f_rad [i]
                omega_y  [i][0] = 2*np.pi*f_vert[i]
                omega_z  [i][0] = 2*np.pi*f_ax  [i]
                omega_bar[i]    = (omega_x[i]*omega_y[i]*omega_z[i])**(1.0/3.0)

            if vary == 'T0' or vary == 'T0&T1' or vary == 'T0=T1':
                # vary temperature of first species to get T/Tc or T/TF = vary_values
                # for MB statistics we take T = vary_values*T[0] nK
                if statistics[0] == 'BoseEinstein':
                    TTcrit[0]    = vary_values
                    Nvar  [0][0] = N[0]
                    Tcrit [0][0] = hbar * omega_bar[0] * ( N[0] / zeta )**(1.0/3.0) *1e9/kB # in nK
                    Tvar  [0]    = TTcrit[0]*Tcrit[0] # nK
                elif statistics[0] == 'FermiDirac':
                    TTcrit[0]    = vary_values
                    Nvar  [0][0] = N[0]
                    Tcrit [0][0] = omega_bar[0]*(6*N[0])**(1.0/3.0)*k_nK # nK
                    Tvar  [0]    = TTcrit[0]*Tcrit[0] # nK
                elif statistics[0] == 'MaxwellBoltzmann':
                    TTcrit[0]    = vary_values # x-axis values
                    Nvar  [0][0] = N[0]
                    Tcrit [0][0] = 0
                    Tvar  [0]    = vary_values*T[0] # nK

                if num_species == 2:
                    if vary == 'T0':
                        # second species with fixed T/Tc or T/TF = vary_fixed
                        # for MB statistics we take T = T[1] in nK
                        if statistics[1] == 'BoseEinstein':
                            TTcrit[1][0]  = vary_fixed
                            Nvar  [1][0] = N[1]
                            Tcrit [1][0] = hbar * omega_bar[1] * ( N[1] / zeta )**(1.0/3.0) *1e9/kB # in nK
                            Tvar  [1][0] = TTcrit[1]*Tcrit[1] # nK
                        elif statistics[0] == 'FermiDirac':
                            TTcrit[1][0] = vary_fixed
                            Nvar  [1][0] = N[1]
                            Tcrit [1][0] = omega_bar[1]*(6*N[1])**(1.0/3.0)*k_nK # nK
                            Tvar  [1][0] = TTcrit[1]*Tcrit[1] # nK
                        elif statistics[1] == 'MaxwellBoltzmann':
                            TTcrit[1]    = vary_values # x-axis values
                            Nvar  [1][0] = N[1]
                            Tcrit [1][0] = 0
                            Tvar  [1][0] = T[1] # nK
                    elif vary == 'T0&T1':
                        # vary temperature also of second species to get T/Tc or T/TF = vary_values
                        # for MB statistics we take T = vary_values*T[1] nK
                        if statistics[1] == 'BoseEinstein':
                            TTcrit[1]    = vary_values
                            Nvar  [1][0] = N[1]
                            Tcrit [1][0] = hbar * omega_bar[1] * ( N[1] / zeta )**(1.0/3.0) *1e9/kB # in nK
                            Tvar  [1]    = TTcrit[1]*Tcrit[1] # nK
                        elif statistics[0] == 'FermiDirac':
                            TTcrit[1]    = vary_values
                            Nvar  [1][0] = N[1]
                            Tcrit [1][0] = omega_bar[1]*(6*N[1])**(1.0/3.0)*k_nK # nK
                            Tvar  [1]    = TTcrit[1]*Tcrit[1] # nK
                        elif statistics[1] == 'MaxwellBoltzmann':
                            TTcrit[1]    = vary_values # x-axis values
                            Nvar  [1][0] = N[1]
                            Tcrit [1][0] = 0
                            Tvar  [1]    = vary_values*T[1] # nK
                    elif vary == 'T0=T1':
                        # second species with T = T of first species
                        if statistics[1] == 'BoseEinstein':
                            Tvar  [1]    = Tvar[0] # nK
                            Tcrit [1][0] = hbar * omega_bar[1] * ( N[1] / zeta )**(1.0/3.0) *1e9/kB # in nK
                            TTcrit[1]    = Tvar[0]/Tcrit[1][0]
                            Nvar  [1][0] = N[1]
                        elif statistics[1] == 'FermiDirac':
                            Tvar  [1]    = Tvar[0] # nK
                            Tcrit [1][0] = omega_bar[1]*(6*N[1])**(1.0/3.0)*k_nK # nK
                            TTcrit[1]    = Tvar[1]/Tcrit[1][0]
                            Nvar  [1][0] = N[1]
                        elif statistics[1] == 'MaxwellBoltzmann':
                            TTcrit[1]    = vary_values # x-axis values
                            Nvar  [1][0] = N[1]
                            Tcrit [1][0] = 0 # in nK
                            Tvar  [1]    = Tvar[0] # nK
            elif vary == 'N0&N1':
                # vary N of both species with vary_values*N[i]
                # 1st species fixed T/Tc or T/TF = vary_fixed. 2nd species T = T of 1st species.
                # for MB statistics we take T = vary_values*T[0] nK
                if statistics[0] == 'BoseEinstein':
                    Nvar  [0]    = (N[0]*vary_values).astype(int)
                    TTcrit[0][0] = vary_fixed
                    Tcrit [0]    = hbar * omega_bar[0] * ( Nvar[0] / zeta )**(1.0/3.0) *1e9/kB # in nK
                    Tvar  [0]    = TTcrit[0][0]*Tcrit[0] # nK
                elif statistics[0] == 'FermiDirac':
                    Nvar  [0]    = (N[0]*vary_values).astype(int)
                    TTcrit[0][0] = vary_fixed
                    Tcrit [0]    = omega_bar[0]*(6*Nvar[0])**(1.0/3.0)*k_nK # nK
                    Tvar  [0]    = TTcrit[0][0]*Tcrit[0] # nK
                elif statistics[0] == 'MaxwellBoltzmann':
                    Nvar  [0]    = (N[0]*vary_values).astype(int)
                    TTcrit[0]    = vary_values # x-axis
                    Tcrit [0][0] = 0
                    Tvar  [0]    = vary_values*T[0] # nK
                if num_species == 2:
                    if statistics[1] == 'BoseEinstein':
                        Nvar  [1]    = (N[1]*vary_values).astype(int)
                        Tcrit [1]    = hbar * omega_bar[1] * ( Nvar[1] / zeta )**(1.0/3.0) *1e9/kB # in nK
                        Tvar  [1]    = Tvar[0] # nK
                        TTcrit[1]    = Tvar[1]/Tcrit[1]
                    elif statistics[1] == 'FermiDirac':
                        Nvar  [1]    = (N[1]*vary_values).astype(int)
                        Tcrit [1]    = omega_bar[1]*(6*Nvar[1])**(1.0/3.0)*k_nK # nK
                        Tvar  [1]    = Tvar[0] # nK
                        TTcrit[1]    = Tvar[1]/Tcrit[1]
                    elif statistics[1] == 'MaxwellBoltzmann':
                        Nvar  [1]    = (N[1]*vary_values).astype(int)
                        Tcrit [1][0] = 0
                        Tvar  [1]    = Tvar[0] # nK
                        TTcrit[1]    = vary_values # x-axis
            else:
                print('vary', vary, 'not implemented!')
                exit()

        # dict with parameters
        params = {}
        print("\nspecies '%s'%s" % (label, ' (single species)' if single_species else ''))
        print("vary '%s'" % (vary))
                
        # expand all dimensions to have same shape of all parameters
        num_var = 0
        for i in range(num_species):
            num_var = np.max([len(omega_x[i]), len(omega_y[i]), len(omega_z  [i]), len(omega_bar[i]),
                              len(Nvar   [i]), len(Tvar   [i]),
                              len(Tcrit  [i]), len(TTcrit [i]),
                              num_var])                          
        if num_var > 1:
            for i in range(num_species):
                if len(omega_x  [i]) == 1: omega_x  [i] = np.array([omega_x  [i][0]]*num_var)
                if len(omega_y  [i]) == 1: omega_y  [i] = np.array([omega_y  [i][0]]*num_var)
                if len(omega_z  [i]) == 1: omega_z  [i] = np.array([omega_z  [i][0]]*num_var)
                if len(omega_bar[i]) == 1: omega_bar[i] = np.array([omega_bar[i][0]]*num_var)
                if len(Nvar     [i]) == 1: Nvar     [i] = np.array([Nvar     [i][0]]*num_var)
                if len(Tvar     [i]) == 1: Tvar     [i] = np.array([Tvar     [i][0]]*num_var)
                if len(Tcrit    [i]) == 1: Tcrit    [i] = np.array([Tcrit    [i][0]]*num_var)
                if len(TTcrit   [i]) == 1: TTcrit   [i] = np.array([TTcrit   [i][0]]*num_var)
            
        # maximum possible number of pairs - used to calculate conversion efficiency
        # note: this can be used only when N is not varied!
        if single_species:
            max_pairs = Nvar[0]/2
        else:
            max_pairs = np.where(Nvar[0] <= Nvar[1], Nvar[0], Nvar[1])

        # save parameters
        
        params["num_species"]               = ["%i", num_species]
        params["species_0"]                 = ["%s", species[0]]
        params["statistics_0"]              = ["%s", statistics[0]]
        params["mass_0"]                    = ["%i", mass[0]]
        params["N_0"]                       = ["%i", Nvar[0]]
        params["T_0"]                       = ["%f", Tvar[0]]
        params["fx_0"]                      = ["%f", omega_x[0]/(2*np.pi)]
        params["fy_0"]                      = ["%f", omega_y[0]/(2*np.pi)]
        params["fz_0"]                      = ["%f", omega_z[0]/(2*np.pi)]
        params["size_0"]                    = ["%f", calc_size[0]]

        if not single_species:
            params["species_1"]             = ["%s", species[1]]
            params["statistics_1"]          = ["%s", statistics[1]]
            params["mass_1"]                = ["%i", mass[1]]
            params["N_1"]                   = ["%i", Nvar[1]]
            params["T_1"]                   = ["%f", Tvar[1]]
            params["fx_1"]                  = ["%f", omega_x[1]/(2*np.pi)]
            params["fy_1"]                  = ["%f", omega_y[1]/(2*np.pi)]
            params["fz_1"]                  = ["%f", omega_z[1]/(2*np.pi)]
            params["cx"]                    = ["%f", offset[0]]
            params["cy"]                    = ["%f", offset[1]]
            params["cz"]                    = ["%f", offset[2]]
            params["size_1"]                = ["%f", calc_size[1]]

        params["statistics_mol"]            = ["%s", statistics[2]]
        params["repeat"]                    = ["%i", repetitions]
        params["error_allowed"]             = ["%e", error_allowed]
        params["step_size"]                 = ["%f", step_size]
        params["threads"]                   = ["%i", threads]
        params["num_bins"]                  = ["%i", num_bins]
        params["find_nearest"]              = ["%i", find_nearest]
        params["gamma"]                     = ["%f", gamma]
        params['distance_measure']          = ['%s', distance_measure]
        params['distribution_generator']    = ['%s', distribution_generator]
        params['random_number_generator']   = ['%s', random_number_generator]
        params["result_file"]               = ["%s", result_file]
        
        if histogram_file is not None:
            params["histogram_file"]        = ["%s", mol_folder + label + histogram_file ]
        if export_atom_0 is not None:
            params["export_atom_0"]         = ["%s", mol_folder + label + export_atom_0  ]
        if not single_species and export_atom_1 is not None:
            params["export_atom_1"]         = ["%s", mol_folder + label + export_atom_1  ]
        if export_molecule is not None:
            params["export_molecule"]       = ["%s", mol_folder + label + export_molecule]
            
        if seed is not None:
            params["seed"]                  = ["%s", seed]

        # print all parameters for all species
        if True:
            for i in range(len(species)):
                print("\nspecies '%s', mass %i amu:" % (species[i], mass[i]))
                #print([len(d) for d in [Nvar[i],Tvar[i],Tcrit[i],TTcrit[i],omega_x[i]/(2*np.pi),omega_y[i]/(2*np.pi),omega_z[i]/(2*np.pi),omega_bar[i]/(2*np.pi)]])
                print_data([Nvar[i],Tvar[i],Tcrit[i],TTcrit[i],omega_x[i]/(2*np.pi),omega_y[i]/(2*np.pi),omega_z[i]/(2*np.pi),omega_bar[i]/(2*np.pi)],
                           ['N','T(nK)','Tcrit(nK)','T/Tcrit','fx(Hz)','fy(Hz)','fz(Hz)','<f>(Hz)'])
                if single_species: break
        
        if recalc:

            # create result folder if does not exist
            os.makedirs(os.path.dirname(mol_folder), exist_ok=True)

            # delete old result file (molConv appends them)
            try:
                os.remove(result_file)
            except FileNotFoundError:
                pass

            # variation of parameters with molConv
            print()
            for var in range(num_var):

                text  = "// parameter file for molConv Monte Carlo\n"
                text += "// command line: molConv.exe -p LiCr_params.txt\n"
                text += "// automatically generated by molConv_run.py\n"
                text += ("// %s\n" % str(datetime.now()))
                text += "// variation %i/%i with %i repetitions\n" % (var+1, num_var, repetitions)
                print(parameter_file)
                text += "// parameter file: " + (parameter_file % var) + "\n\n"
                print(text)
                    
                # save parameter file
                for p,value in params.items():
                    if isinstance(value[1], (list, np.ndarray)): 
                        xi = value[1][var]
                    else:
                        xi = value[1]
                    xi = (value[0] % xi)
                    print("parameter '%s' = %s" % (p, xi))
                    text += (p + " = " + xi + "\n")

                with open(parameter_file % (var), 'w') as f:
                    f.write(text)
                        
                # run molConv
                # note: this gives result = 255 for molConv returning -1
                result = subprocess.run([molConv_path, "-p", "%s" % (parameter_file % (var))]).returncode
                
                if result != 0:
                    print('molConv error', result, '\n')
                    break
                else:
                    print('molConv ok\n')
        else:
            result = 0

        if result != 0:
            break
        else:
            print('parse result file:', result_file)
        
            # read result file
            [result_text, result_order, result_param, lines] = read_file(result_file, results, allow_multiple=True, show=False, force_unit=True)

            if False:            
                print(result_order)
                print(result_param)

            # parameters in order they appear in file
            rep_  = [result_param[i] for i in find_all(result_order, 'repetition')]
            sp_   = [result_param[i] for i in find_all(result_order, 'atomic species')]
            N_    = [result_param[i] for i in find_all(result_order, 'N')]
            T_    = [result_param[i] for i in find_all(result_order, 'T')]
            Tc_   = [result_param[i] for i in sorted(find_all(result_order, 'Tc'  ) + find_all(result_order, 'TF'  ))]
            TTc_  = [result_param[i] for i in sorted(find_all(result_order, 'T/Tc') + find_all(result_order, 'T/TF'))]
            mol_  = [result_param[i] for i in find_all(result_order, 'molecule')]
            Nmol_ = [result_param[i] for i in find_all(result_order, 'N_mol')]
            emol_ = [result_param[i] for i in find_all(result_order, 'conv. efficiency')]
            tc_   = [result_param[i] for i in find_all(result_order, 'calculation time')]
            tt_   = [result_param[i] for i in find_all(result_order, 'calculation time [total]')]
            ts_   = [result_param[i] for i in find_all(result_order, 'calculation time [summary]')]
                        
            # parameters per species, variation and repetition
            N_rep     = [np.zeros(shape=(num_var,repetitions), dtype=float) for _ in range(num_species)]
            T_rep     = [np.zeros(shape=(num_var,repetitions), dtype=float) for _ in range(num_species)]
            Tc_rep    = [np.zeros(shape=(num_var,repetitions), dtype=float) for _ in range(num_species)]
            TTc_rep   = [np.zeros(shape=(num_var,repetitions), dtype=float) for _ in range(num_species)]
            t_s_rep   = [np.zeros(shape=(num_var,repetitions), dtype=float) for _ in range(num_species)]
            t_m_rep   =  np.zeros(shape=(num_var,repetitions), dtype=float)
            t_tot_rep =  np.zeros(shape=(num_var,repetitions), dtype=float)            
            Nmol_rep  =  np.zeros(shape=(num_var,repetitions), dtype=float)
            emol_rep  =  np.zeros(shape=(num_var,repetitions), dtype=float)

            # walk through all variations and repetitions, check consistency and assign parameters
            cr    = 0   # repetition + variation counter (rep)
            cs    = 0   # species counter (T, N)
            css   = 0   # summary of species counter (Tc, TTc)
            cm    = 0   # molecule counter (Nmol, emol)
            cc    = 0   # caclulation time counter (t_s, t_m, t_tot)
            for v in range(num_var):
            
                for r in range(repetitions):
            
                    # check repetitions
                    if cr >= len(rep_) or rep_[cr] != r+1:
                        print(len(rep_), cr, repetitions, num_species, num_var, v)
                        print('%i/%i error repetition (index %i):'%(v,r,cr), rep_[cr], '!=', r+1)
                        print(rep_)
                        exit()

                    # check species
                    # note: for bosons condensed part might be missing when not condended.
                    species_N = [[0]*len(n) for n in species_names]
                    species_t = [[0]*len(n) for n in species_names]
                    for s in range(num_species):
                        num = 0
                        for si,name in enumerate(species_names[s]):
                            if cs >= len(sp_) or sp_[cs] != name:
                                if statistics[s] == 'BoseEinstein' and name.endswith('condensed)'):
                                    print('%i/%i species:'%(v,r), name, 'missing (assume Nc=0)')
                                    continue
                                else:
                                    print('%i/%i error species:'%(v,r), sp_[cs], '!=', name)
                                    exit()
                            elif cc >= len(tc_):
                                print('%i/%i error species %s calculation time index out of range:'%(v,r,name), cc, '>=', len(tc_))
                                exit()
                            else:
                                species_N[s][si] = N_[cs]
                                species_t[s][si] = tc_[cc] # creation time for individual sub-species (thermal, condensed)
                                num += 1
                                cs += 1
                                cc += 1
                        
                        # check species summary
                        # here we check that N is the sum of same species before 
                        # and T is the same as given before
                        if cs >= len(sp_) or sp_[cs] != species[s]+' [summary]':
                            print('%i/%i error species %i:'%(v,r,s), sp_[cs], '!=', species[s])
                            exit()
                        elif cs >= len(N_) or N_[cs] != sum(species_N[s]):
                            print('%i/%i error species %i N:'%(v,r,s), N_[cs], '!=', sum(species_N[s]))
                            exit()
                        elif cs >= len(T_) or np.any(np.abs(np.array(T_[cs-num: cs]) - T_[cs]) != 0):
                            print('%i/%i error species %i T:'%(v,r,s), T_[cs], '!=', T_[cs-num: cs])
                            print('error:', np.array(T_[cs-num: cs]) - T_[cs])
                            exit()
                        elif cc >= len(tc_) or (tc_[cc]+t_err) < sum(species_t[s]):
                            print('%i/%i error species %i sum of calculation time:'%(v,r,s), tc_[cc], '<', sum(species_t[s]))
                            print(species_t, cc, len(tc_))
                            exit()
                        elif statistics[s] != 'MaxwellBoltzmann':
                            if Tc_[css] != round(Tcrit[s][v],3):
                                print('%i/%i error species %i Tc or TF (index %i):'%(v,r,s,css), Tc_[css], '!=', round(Tcrit[s][v],3))
                                exit()
                            elif TTc_[css] != round(T_[cs]/Tc_[css],3):
                                print('%i/%i error species %i T/Tc or T/TF (calc.):'%(v,r,s), TTc_[css], '!=', round(T_[cs]/Tc_[css],3))
                                exit()
                            elif TTc_[css] != round(TTcrit[s][v],3):
                                print('%i/%i error species %i T/Tc or T/TF (value):'%(v,r,s), TTc_[css], '!=', round(TTcrit[s][v],3))
                                exit()
                            Tc_rep [s][v, r] = Tc_ [css]
                            TTc_rep[s][v, r] = TTc_[css]
                            css += 1
                        else:
                            # for MB distribution there is no Tc or TF in summary
                            # take the T/Tc or T/TF from variation for x-axis of plot
                            Tc_rep [s][v, r] = Tcrit [s][v]
                            TTc_rep[s][v, r] = TTcrit[s][v]
                            
                        N_rep  [s][v, r] = N_ [cs]
                        T_rep  [s][v, r] = T_ [cs]
                        t_s_rep[s][v, r] = tc_[cc] # total creation time of species
                        cs += 1
                        cc += 1

                    # check molecule entries
                    # we check also conversion efficiency calculation from N of each species combination i,j
                    num = 0
                    mol_t = [0]*len(mol_names)
                    for mm, (name,i,j) in enumerate(mol_names):
                        if single_species:
                            pairs = int(species_N[0][i]//2) if (i==j) else min(species_N[0][i], species_N[0][j])
                        else:
                            pairs = min(species_N[0][i], species_N[1][j])
                        if cm >= len(mol_) or mol_[cm] != name:
                            if statistics[0] == 'BoseEinstein' and species_names[0][i].endswith('condensed)') and species_N[0][i] == 0:
                                print('%i/%i molecule:'%(v,r), name, 'missing (species 0 non-condensed, ok)')
                            elif statistics[1] == 'BoseEinstein' and species_names[num_species-1][j].endswith('condensed)') and species_N[num_species-1][j] == 0:
                                print('%i/%i molecule:'%(v,r), name, 'missing (species 1 non-condensed, ok)')
                            else:
                                print('%i/%i error molecule:'%(v,r), mol_[cm], '!=', name)
                                print('molecule names are:', mol_names)
                                print(species_N, i, j)
                                exit()
                        elif cm >= len(emol_) or cm >= len(Nmol_) or emol_[cm] != round(Nmol_[cm]/pairs*100,3):
                            print('%i/%i error molecule %s efficiency:'%(v,r,name), emol_[cm], '!=', round(Nmol_[cm]/pairs*100,3))
                            print('efficiency = N_mol/pairs = ', Nmol_[cm], '/', pairs, '=', Nmol_[cm]/pairs)
                            exit()
                        elif cc >= len(tc_):
                            print('%i/%i error molecule %s calculation time index out of range:'%(v,r,name), cc, '>=', len(tc_))
                            exit()
                        else:
                            mol_t[mm] = tc_[cc] # creation time for each molecule
                            num += 1
                            cm += 1
                            cc += 1

                    # check molecule summary
                    # here we check that N_mol is the sum of individual molecules (thermal+BEC) given before
                    # and we check conversion efficiency calculation = N_mol / max_pairs
                    name = species[0] + species[1] + ' [summary]'
                    pairs = sum(species_N[0])/2 if (single_species) else min(sum(species_N[0]), sum(species_N[1]))
                    if cm >= len(mol_) or mol_[cm] != name:
                        print('%i/%i error molecule (summary):'%(v,r), mol_[cm], '!=', name)
                        print('molecule names are:', mol_names)
                        exit()
                    elif vary != 'N0&N1' and pairs != max_pairs[v]:
                        print('%i/%i error max. number pairs:'%(v,r), pairs, '!=', max_pairs[v], '\nwas N varied?')
                        print(species_N)
                        exit()                        
                    elif cm >= len(Nmol_) or Nmol_[cm] != np.sum(Nmol_[cm-num:cm]):
                        print('%i/%i error molecule N_mol:'%(v,r), Nmol_[cm], '!=', np.sum(Nmol_[cm-num:cm]))
                        exit()
                    elif cm >= len(emol_) or emol_[cm] != round(Nmol_[cm]/pairs*100,3):
                        print('%i/%i error molecule efficiency:'%(v,r), emol_[cm], '!=', round(Nmol_[cm]/pairs*100,3))
                        print('efficiency = N_mol/pairs = ', Nmol_[cm], '/', pairs, '=', Nmol_[cm]/pairs)
                        exit()
                    elif cc >= len(tc_) or (tc_[cc]+t_err) < sum(mol_t):
                        print('%i/%i error molecule %s sum of calculation time:'%(v,r,name), tc_[cc], '<', sum(mol_t))
                        print(mol_t, cc, len(tc_))
                        exit()
                    Nmol_rep[v, r] = Nmol_[cm]
                    emol_rep[v, r] = emol_[cm]
                    t_m_rep [v, r] = tc_  [cc]
                    cm += 1
                    cc += 1

                    # total calculation time = sum of atoms + molecules calculation time
                    num = sum([sum(s) for s in species_t]) + sum(mol_t)
                    if cr >= len(tt_) or (tt_[cr]+t_err) < num:
                        print('%i/%i error total calculation time:'%(v,r), tt_[cr], '<', num)
                        print(species_t, mol_t, cr, len(tt_))
                        exit()
                    t_tot_rep[v, r] = tt_[cr]

                    # next iteration/variation
                    cr += 1

            # check total length. 
            # this ensures that there is not more data in file than expected.
            if len(rep_) != cr:
                print('%i/%i error repetition length:'%(v,r), len(rep_), '!=', cr)
                exit()
            elif len(sp_) != cs:
                print('%i/%i error species length:'%(v,r), len(sp_), '!=', cs)
                exit()
            elif len(N_) != cs:
                print('%i/%i error N length:'%(v,r), len(N_), '!=', cs)
                exit()
            elif len(T_) != cs:
                print('%i/%i error T length:'%(v,r), len(T_), '!=', cs)
                exit()
            elif len(Tc_) != css:
                print('%i/%i error Tc length:'%(v,r), len(Tc_), '!=', css)
                exit()
            elif len(TTc_) != css:
                print('%i/%i error T/Tc length:'%(v,r), len(TTc_), '!=', css)
                exit()
            elif len(mol_) != cm:
                print('%i/%i error molecule length:'%(v,r), len(mol_), '!=', cm)
                exit()
            elif len(Nmol_) != cm:
                print('%i/%i error molecule N length:'%(v,r), len(Nmol_), '!=', cm)
                exit()
            elif len(emol_) != cm:
                print('%i/%i error molecule eff length:'%(v,r), len(emol_), '!=', cm)
                exit()
            elif len(tc_) != cc:
                print('%i/%i error calculation time length:'%(v,r), len(tc_), '!=', cc)
                exit()
            elif len(tt_) != cr:
                print('%i/%i error calculation time total length:'%(v,r), len(tt_), '!=', cr)
                exit()
            elif len(ts_) != cr/repetitions:
                print('%i/%i error calculation time summary length:'%(v,r), len(ts_), '!=', cr/repetitions)
                exit()
                
            print('result file ok')
            
            # get mean and standard deviation for each variation
            N_avg      = [np.mean(N_rep     [s], axis=1) for s in range(num_species)]
            T_avg      = [np.mean(T_rep     [s], axis=1) for s in range(num_species)]
            Tcrit_avg  = [np.mean(Tc_rep    [s], axis=1) for s in range(num_species)]
            TTcrit_avg = [np.mean(TTc_rep   [s], axis=1) for s in range(num_species)]
            t_s_avg    = [np.mean(t_s_rep   [s], axis=1) for s in range(num_species)]
            t_m_avg    =  np.mean(t_m_rep      , axis=1)
            t_tot_avg  =  np.mean(t_tot_rep    , axis=1)
            Nmol_avg   =  np.mean(Nmol_rep     , axis=1)
            emol_avg   =  np.mean(emol_rep     , axis=1)
            N_err      = [np.std (N_rep     [s], axis=1) for s in range(num_species)]
            T_err      = [np.std (T_rep     [s], axis=1) for s in range(num_species)]
            Tcrit_err  = [np.std (Tc_rep    [s], axis=1) for s in range(num_species)]
            TTcrit_err = [np.std (TTc_rep   [s], axis=1) for s in range(num_species)]
            t_s_err    = [np.std (t_s_rep   [s], axis=1) for s in range(num_species)]
            t_m_err    =  np.std (t_m_rep      , axis=1)
            t_tot_err  =  np.std (t_tot_rep    , axis=1)
            Nmol_err   =  np.std (Nmol_rep     , axis=1)
            emol_err   =  np.std (emol_rep     , axis=1)
            
            # export result file
            if num_species == 1:
                csv_file_header = 'N,T/nK,Tcrit/nK,T/Tcrit,N_mol,N_mol_error,eff_mol,eff_mol_error,t_atom/s,t_atom_error/s,t_mol/s,t_mol_error/s,t_tot/s,t_tot_error/s\n'
                csv_data = [N_avg[0], 
                            T_avg[0], 
                            Tcrit_avg[0], 
                            TTcrit_avg[0], 
                            Nmol_avg, Nmol_err, 
                            emol_avg, emol_err, 
                            t_s_avg[0], t_s_err[0], 
                            t_m_avg, t_m_err, 
                            t_tot_avg, t_tot_err]
            else:
                csv_file_header = 'N_0,N_1,T_0/nK,T_1/nK,Tcrit_0/nK,Tcrit_1/nK,T/Tcrit_0,T/Tcrit_1,N_mol,N_mol_error,eff_mol,eff_mol_error,t_0/s,t_0_error/s,t_1/s,t_1_error/s,t_mol/s,t_mol_error/s,t_tot/s,t_tot_error/s\n'
                csv_data = [N_avg[0], N_avg[1], 
                            T_avg[0], T_avg[1], 
                            Tcrit_avg[0], Tcrit_avg[1], 
                            TTcrit_avg[0], TTcrit_avg[1], 
                            Nmol_avg, Nmol_err, 
                            emol_avg, emol_err, 
                            t_s_avg[0], t_s_err[0], 
                            t_s_avg[1], t_s_err[1], 
                            t_m_avg, t_m_err, 
                            t_tot_avg, t_tot_err]        

            print(label, 'average result:')
            print(csv_file_header)
            print(np.transpose(csv_data))

            fmt = '%.6f'
            csv = csv_file_header
            for i in range(num_var):
                csv += (','.join([fmt%d[i] for d in csv_data]))+'\n'
            with open(csv_file, 'w') as f:
                f.write(csv)
            print("result written to '%s'" % csv_file)

            # collect plotting data
            lm = r'%s, '%label
            ld = distance_measure_all[distance_measure]
            lg = r', $\gamma$=%.1e'%gamma
            if figure == 'timing':
                x = N_avg[0]
                d            = [x, t_s_avg[0], t_s_err[0], t_s_avg[1], t_s_err[1], t_m_avg, t_m_err, t_tot_avg, t_tot_err]
                data        += [[x, t_s_avg[0], t_s_err[0]],
                                [x, t_s_avg[1], t_s_err[1]], 
                                [x, t_m_avg   , t_m_err   ],
                                [x, t_tot_avg , t_tot_err]]
                data_labels += [lm+species[0], lm+species[1], lm+species[0] + species[1], lm+'total']
                data_args   += data_args_mol 

            elif figure == 'threads':
                if len(data_labels) == 0:
                    data_labels += [species[0], species[1], species[0] + species[1], 'total']
                    data_args   += data_args_mol
                    data         = [[[],[],[]], [[],[],[]], [[],[],[]], [[],[],[]]]
                data[0][0] += [threads]
                data[0][1] += [t_s_avg[0][0]]
                data[0][2] += [t_s_err[0][0]]
                data[1][0] += [threads]
                data[1][1] += [t_s_avg[1][0]] 
                data[1][2] += [t_s_err[1][0]]
                data[2][0] += [threads]
                data[2][1] += [t_m_avg   [0]]
                data[2][2] += [t_m_err   [0]]
                data[3][0] += [threads]
                data[3][1] += [t_tot_avg [0]]
                data[3][2] += [t_tot_err [0]]

            elif figure == 'test_data':
                # test data list, x-axis = series number, y-axis = molecule conversion efficiency comparison calculated vs. measured
                x = np.arange(len(emol_avg))
                data        += [[x, emol_avg, emol_err],[x, np.array(td_N_mol)/max_pairs*100.0], [x, np.array(td_f_mol)*100.0]]
                data_labels += [lm+ld+lg+' simulated', 'measured Nmol/max_pairs', 'measured']
                data_args   += data_args_mol 

            else:
                d            = [TTcrit_avg[0], emol_avg, emol_err]
                data        += [d]
                data_labels += [lm+ld+lg]
                data_args   += data_args_mol 

                if num_species == 2:
                    twin_data   += [[TTcrit_avg[0], TTcrit_avg[1]]]
                    twin_labels += [species[1]]
                    twin_args   += [{'color':args['color'], 'edgecolor':args['color'], 'facecolor':'White'} for args in data_args_mol]

    if result == 0:

        if figure == 'timing':
            # timing figure, x-axis = number of atoms, y-axis = calculation time
            ax = plot(
                    title        = title, 
                    data         = data, 
                    data_labels  = data_labels,
                    data_args    = data_args, 
                    curves       = [],   
                    curve_labels = None,
                    curve_args   = None, 
                    xlabel       = 'atom number',
                    x_range      = None,
                    ylabel       = 'calculation time (s)',
                    y_range      = None,
                    log_scale    = [True,True],
                    label_pos    = 'upper left', 
                    label_cols   = 1,
                    fig_size     = (10*4/3,7), # figure (width,height)
                    fig_pos      = [0.07, 0.08, 0.92, 0.87] # sub plot (left,bottom,width,height)
                    ) 

        elif figure == 'threads':
            # threads performance figure, x-axis = number of threads, y-axis = calculation time and gain
            print(data)
            ax = plot(
                    title        = title, 
                    data         = data, 
                    data_labels  = data_labels,
                    data_args    = data_args, 
                    curves       = [],   
                    curve_labels = None,
                    curve_args   = None, 
                    xlabel       = 'threads',
                    x_range      = None,
                    ylabel       = 'calculation time (s)',
                    y_range      = None,
                    log_scale    = [False,True],
                    label_pos    = 'upper right', 
                    label_cols   = 1,
                    fig_size     = (10*4/3,7), # figure (width,height)
                    fig_pos      = [0.07, 0.08, 0.87, 0.87] # sub plot (left,bottom,width,height)
                    ) 
            # get thread efficiency
            twin_labels = data_labels
            twin_args   = [{'color':da['color'], 'edgecolor':da['color'], 'facecolor':'White'} for da in data_args]
            twin_data   = [[np.array(d[0]), (d[1][0]/np.array(d[0]))/np.array(d[1])*100.0] for d in data]
            print(twin_args)
            print(twin_data)
            plot(   twin         = ax,
                    data         = twin_data, 
                    data_labels  = twin_labels,
                    data_args    = twin_args, 
                    curves       = [],   
                    curve_labels = None,
                    curve_args   = None, 
                    ylabel       = 'thread performance (%)',
                    y_range      = None,
                    label_pos    = None,
                    ) 
        
        elif figure == 'test_data':
            # test data list, x-axis = series number, y-axis = molecule conversion efficiency comparison calculated vs. measured
            ax = plot(
                    title        = title, 
                    data         = data, 
                    data_labels  = data_labels,
                    data_args    = data_args, 
                    curves       = [],   
                    curve_labels = None,
                    curve_args   = None, 
                    xlabel       = 'test data #',
                    x_range      = None,
                    ylabel       = 'conversion efficiency (%)',
                    y_range      = None,
                    log_scale    = [False,True],
                    label_pos    = 'upper left', 
                    label_cols   = 1,
                    fig_size     = (10*4/3,7), # figure (width,height)
                    fig_pos      = [0.07, 0.08, 0.92, 0.87] # sub plot (left,bottom,width,height)
                    ) 
        else:
            # other figures: x-axis = T/TF or T/Tc, y-axis: molecule creation efficiency
            ax = plot(
                    title        = title, 
                    data         = data, 
                    data_labels  = data_labels,
                    data_args    = data_args, 
                    curves       = [],   
                    curve_labels = None,
                    curve_args   = None, 
                    xlabel       = 'T/Tc or T/TF (first species)', #label_x,
                    x_range      = [-0.05,1.65],
                    ylabel       = 'efficiency (%)', #'efficiency (%' + (', solid symbols)' if show_second_species and num_species == 2 else ')'),
                    y_range      = [-5,105],
                    label_pos    = 'upper right', 
                    label_cols   = 1,
                    fig_size     = (10*4/3,7), # figure (width,height)
                    fig_pos      = [0.07, 0.08, 0.87, 0.87] # sub plot (left,bottom,width,height)
                    ) 
            if num_species == 2:
                plot(   twin         = ax,
                        data         = twin_data, 
                        data_labels  = twin_labels,
                        data_args    = twin_args, 
                        curves       = [],   
                        curve_labels = None,
                        curve_args   = None, 
                        ylabel       = 'T/Tc or T/TF (second species, open symbols)',
                        y_range      = [-0.05,1.05, [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]],
                        label_pos    = None,
                        ) 

        plt.show()
            
        
################################################################################################
# calculation ends here
################################################################################################

