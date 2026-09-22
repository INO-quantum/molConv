#!/usr/local/bin/python3

# checks distribution in histograms and atoms/molecule files
# run molConv_run.py and export histograms and atoms/molecules files
# give here parameter file name and repetition number (>=1) to be analyzed
# plots distribution of all atoms and molecule obtained from histogram files or/and atoms files
# and compares with expected distribution using parameter file.
# if load_atoms == True generates additionally histogram directly from atoms/molecule files.
# last changed 22/9/2026 by andi

# TODO: at the moment reads only histogram files and does not load atoms csv files

param_file  = './tmp/test_delta_x*delta_p_3.8e-01/test_params_0.txt'
repetition  = 1
load_atoms  = False # has no effect at the moment!

# figure identifier
fig_x = 'position'
fig_v = 'velocity'
fig_E = 'energy'

# figures and species to be plotted. None = all
fig_sel     = [fig_x, fig_v]
fig_species = ['Li6', 'Cr52']

from molConv_run import plot
import numpy as np
#from numpy.random import uniform, randint, normal, multivariate_normal
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from datetime import datetime
from mpmath import polylog, re, im

mpl.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 2
mpl.rcParams["savefig.directory"] = os.path.dirname(__file__) # save to current directy

# constants. updated 3/8/2025 mn (old value was 1.67492749804e-27)
# TODO: I use here mn but also amu could be used which is slightly (<1%) different. not clear which is more appropriate? 
hPlanck = 6.62607015e-34                         # Planck constant in Js
hbar    = 1.054571817e-34                        # reduced Planck constant in Js
kB      = 1.380649e-23                           # Boltzmann constant in J/K
mn      = 1.67492750056e-27                      # neutron mass in kg
a_Bohr  = 5.29177210544e-11                      # Bohr radius in m
k_Epot  = (1.67492750056/1.380649     )*1e-7     # get [Epot] = amu*(rad/s*mu)^2 -> nK; mn*(1e-6)^2/kB*1e9 = mn/kB*1e-3;
k_Ekin  = (1.67492750056/1.380649     )*1e-1     # get [Ekin] = amu*(mm/s)^2     -> nK; mn*(1e-3)^2/kB*1e9 = mn/kB*1e3;
k_nK    = (1.054571817  /1.380649     )*1e-2     # get TF and Tc in nK out of hbar*omega: k_nK = (hbar*1e9)/kB;
k_h     = (6.62607015   /1.67492750056)*1e2      # scale phase space density (gamma): k_h = h/(mn*1e-9)
k_Rc    = np.sqrt(1.054571817/1.67492750056*1e5) # get BEC size in mu = 1e6*np.sqrt(hbar/mn)

# test polylog. value checked with Mathematica 10.0 (3/8/2025)
zeta = float(re(polylog(3,1)))
error = abs(zeta-1.2020569031595942854)
if error > 1e-16:
    print('zeta(3) = Li_3(1) =', zeta, ', error = %.3e > 1e-16!' % error)
    exit()
else:
    print('zeta(3) = Li_3(1) =', zeta, ', error = %.3e (ok)' % error)

def get_int(text):
    try:
        return int(text)
    except ValueError:
        print("cannot interpret '", text, ", as integer!")
        exit()

# histogram header entry, unit, column index, column name, color, figure
hist_header = [
        ['// stat'    , ''    , None, None  , None   , None  ],
        ['// rep'     , ''    , None, None  , None   , None  ],
        ['// N'       , ''    , None, None  , None   , None  ],
        ['// T'       , 'nK'  , None, None  , None   , None  ],
        ['// mu'      , 'nK'  , None, None  , None   , None  ],
        [None         , None  , 0   , '#'   , 'Black', None  ],
        ['// max x'   , 'um'  , 1   , 'x'   , 'Red'  , fig_x ],
        ['// max y'   , 'um'  , 2   , 'y'   , 'Blue' , fig_x ],
        ['// max z'   , 'um'  , 3   , 'z'   , 'Green', fig_x ],
        ['// max vx'  , 'mm/s', 4   , 'vx'  , 'Red'  , fig_v ],
        ['// max vy'  , 'mm/s', 5   , 'vy'  , 'Blue' , fig_v ],
        ['// max vz'  , 'mm/s', 6   , 'vz'  , 'Green', fig_v ],
        ['// max Etot', 'nK'  , 7   , 'Etot', 'Red'  , fig_E ],
        ['// max Etot', 'nK'  , 8   , 'Epot', 'Blue' , fig_E ],
        ['// max Etot', 'nK'  , 9   , 'Ekin', 'Green', fig_E ],
    ]
hist_header_entry   = [h[0] for h in hist_header]
hist_units          = [h[1] for h in hist_header]
hist_colum_index          = [h[2] for h in hist_header]
hist_colum_name     = [h[3] for h in hist_header]
hist_colors         = [h[4] for h in hist_header]
hist_figure         = [h[5] for h in hist_header]
hist_header_len     = len(set([h for h in hist_header_entry if h is not None])) + 3
hist_num_cols       = len([h for h in hist_colum_index if h is not None])

def load_params(filename, skip=['//', '#'], max_lines=None):
    "loads paramaters as dictionary from file skipping all lines starting with skip until max_lines reached"
    params = {}
    index = 0
    with open(filename, 'r') as f:
        while True:
            line = f.readline()
            if len(line) == 0: break
            else:
                line = line.strip()
                if len(line) > 0:
                    ok = True
                    for s in skip:
                        if line.startswith(s): 
                            ok = False
                            break
                    if ok:
                        try:
                            pos   = line.index('=')
                            key   = line[:pos].strip()
                            value = line[pos+1:].strip()
                            s = value.split(' ')
                            if len(s) == 1:
                                unit = ''
                            elif len(s) == 2:
                                value = s[0]
                                unit  = s[1]
                            else:
                                print('line%i'% index, key, 'unexpected value:', value)
                                exit()
                            try:
                                value = float(value)
                            except ValueError:
                                pass
                            params[key] = [value, unit]
                            print('%-25s'%key, value, unit)
                        except ValueError:
                            pass
            index += 1
            if max_lines is not None and index >= max_lines: break
    return params

def FermiDirac(E, T, mu, Emin=None):
    # see Zwierlein MakingProbingUnderstanding, Equ. 17, lower sign
    return (1.0/(np.exp((E-mu)/T)+1.0))

def BoseEinstein_thermal(E, T, mu, Emin):
    # see Zwierlein MakingProbingUnderstanding, Equ. 17, upper sign
    # note: mu <= 0.0
    # we must limit E >= Emin (ground state energy) otherwise function diverges!
    return np.where(E < Emin, 0.0, 1.0/(np.exp((E-mu)/T)-1.0))

def BoseEinstein_condensed(E, T, mu, Emin):
    # see T. Yamakoshi, S. Watanabe, Ch. Zhang, and Ch. Greene, "Stochastic and equilibrium pictures of the ultracold Fano-Feshbach-resonance molecular conversion rate", Physical Review A 87, 053604 (2013)
    return np.exp(-E/Emin)

def MaxwellBoltzmann(E, T, mu, Emin=None):
    # thermal gas
    return np.exp(-E/T)

def TFermi(N, omega_bar): 
    # Fermi temperature in nK in harmonic trap
    # see Zwierlein MakingProbingUnderstanding, Equ. 33
    return omega_bar*(6*N)**(1.0/3.0)*k_nK; # nK
    
def Tcrit_(N, omega_bar): 
    # BEC critical temperature in nK in harmonic trap
    # see: Ketterle MakingProbingUnderstanding BEC, Equ. 51 (note: Equ. in text on top of page 40 is wrong!)
    #      Zwierlein MakingProbingUnderstanding Fermi Gases, Equ. 29
    #      zeta_3=Li_3(1) must be inside sqrt.
    return omega_bar*((N/zeta)**(1.0/3.0))*k_nK; # nK. k_nK = (hbar*1e9)/kB

dist_func  = {'FermiDirac'             : FermiDirac, 
              'BoseEinstein'           : BoseEinstein_thermal, 
              'BoseEinstein(thermal)'  : BoseEinstein_thermal,
              'BoseEinstein(condensed)': BoseEinstein_condensed, 
              'MaxwellBoltzmann'       : MaxwellBoltzmann       }
Tcrit_func = {'FermiDirac'             : TFermi, 
              'BoseEinstein'           : Tcrit_, 
              'MaxwellBoltzmann'       : lambda N,omega_bar: 0.0}
    
def ChemicalPotential(info, dist, N, T, omega, mass, error_mu=1e-3, error_N=0.1, error_Emax = 1.0, epsilon_N=None, error_Nint=0.01, error_BECfrac=0.01, kdE=0.8, tmax=5.0, show=False):
    # get chemical potential such that integrating distribution gives atom number N
    # info          = species name used for printing
    # dist          = 'MaxwellBoltzmann' or 'MB', 'FermiDirac' or 'FD', 'BoseEinstein' or 'BE', 'BEC'
    # N             = atom number
    # T             = temperature in nK
    # omega         = trap frequency [x,y,z] in rad/s
    # error_mu      = goal error in mu in nK
    # error_N       = goal error in atom number
    # error_Emax    = goal error of Emax. this is not critical.
    # epsilon_N     = resolution in atom number in order to get Emax. if is None uses error_N/N.
    # error_Nint    = criterion for dE loop: np.abs(Nint - Nold)/Nint <= error_Nint.
    # error_BECfrac = gives error when mu=0 and |Nc/N-frac| > error_BECfrac.
    # kdE           = multiplication factor for dE when search direction changes or mu = 0.
    # tmax          = maximum time in seconds before function is aborted.
    # show          = if True prints output during iteration, if False prints only result at end, if None prints nothing
    # returns mu, Tc, Nint, Rth, Rc, Emax
    # note: 
    # - Rth    = for all distributions = sigma Gaussian of thermal (MB) distribution
    # - for MB = Maxwell-Boltzmann statistics returns mu = Tc = Rc = 0.0, Nint = N
    # - for FD = Fermi-Dirac statistics returns Tc = EF, Nint ~ N, Rc = RF
    # - for BE = Bose-Einstein non-condensed statistics returns Tc, mu <= 0.0, Rc = condensate size
    #            if T<=Tc returns mu=0 and Nint = thermal fraction, N-Nint = BEC fraction 
    #            if T> Tc returns mu<0 and Nint ~ N = thermal fraction 
    # - for BEC = Bose-Einstein condensed statistics: returns Tc, mu = 0.0, Rc = condensate size
    # TODO: 
    # - here we assume harmonic trap! 
    # - check the prefactors

    # epsilon_N defines Emax
    # TODO: verify that this is good in all cases with real distribution for different N! 
    #       for high T gives quite high Emax and for T<Tc gives low Emax.
    if epsilon_N is None:
        epsilon_N = error_N/N
    
    # average trap frequency in rad/s
    omega_bar = (omega[0]*omega[1]*omega[2])**(1/3)
    
    # ground-state energy in nK
    Emin = 0.5*omega_bar*k_nK

    # thermal size = Gauss sigma
    Rth = np.sqrt(T/(k_Epot*mass))/omega  # um. [1/2*kB*T = m*omega^2*x^2/2]
    
    if dist == 'MB' or dist == 'MaxwellBoltzmann': # Maxwell-Boltzmann: mu = Tc = 0.0
        Emax = T*np.log(N/epsilon_N)
        if show is not None:
            print(info, '%s: mu = 0.0, Emax = %.3f' % (dist, Emax))
        return 0.0, 0.0, N, Rth, 0.0, Emax
    elif dist == 'FD' or dist == 'FermiDirac': # Fermi-Dirac: Tc = EF
        dist_func = FermiDirac
        Tc = TFermi(N, omega_bar); # nK
        Rc = np.sqrt(2.0*T/(k_Epot*mass))/omega  # Radius in Polylog in um. [kB*T = m*omega^2*x^2/2]
        mu = Tc # start with Fermi energy
        #print('mu start %.3f' % mu)
        dE = Tc # initial energy steps
        limit = Tc # limit mu <= TF. Natoms should be still reached.
        Emax = Tc+T
    elif dist == 'BE' or dist == 'BEC' or dist == 'BoseEinstein': # Bose-Einstein: Tc, mu <= 0
    
        dist_func = BoseEinstein_thermal
        Tc = Tcrit_(N, omega_bar); # nK
        limit = 0.0 # limit mu <= 0.0
        # Radius of non-condensed part in Polylog in um. [kB*T = m*omega^2*x^2/2]
        Rth = np.sqrt(2.0*T/(k_Epot*mass))/omega
        # Radius of condensed part = ground state of non-interacting HO in um. 
        # see Ketterle, Making-Probing-Understanding, Equ. 38. [k_Rc = 1e6*np.sqrt(hbar/mn)] 
        # factor 0.5 to match Gauss-sigma factor 1/2 which is not in Equ. 38.
        Rc = np.sqrt(0.5/(mass*omega))*k_Rc
        if T < Tc: # mu close to 0: see Zwierlein, Making-Probing-Understanding Fermi Gases, Eq. 30 
            Nc = N*(1-(T/Tc)**3)
            mu = -T/Nc # start with Tc/N0, see Zwierlein M-P-U, p.53.
        else: # mu < 0
            mu = -Tc # start with critical temperature
        dE  = Tc # initial energy steps
        Emax = Tc+T

        # if scattering length is given calculate chemical potential [not tested so far]
        # see Ketterle Making Probing Understanding, Eq. 50
        # units: [hbar^2*sqrt(mn)*a0]^2/5*10^9/kB = [hbar^4*mn*a0^2*(10^9/kB)^5]^1/5 = [k_nK^4 * k_Ekin * 10^6*a0^2]^1/5
        #        with k_nK = (hbar*1e9)/kB, k_Ekin = mn/kB*1e3
        #try:
        #    mu = np.sign(a_BEC[info]) * ( (15 * N * omega_bar**3 * a_BEC[info])**2 * mass/32 * 1e6 * k_nK**4 * k_Ekin * a_Bohr**2 )**(1/5)
        #    print(info, 'note: ChemicalPotential calculated from a = %.3f a0 gives mu = %.3f nK' % (a_BEC[info], mu))
        #    return mu, Tc, 0.0, Rth, Rc, Emax
        #except KeyError:
        #    pass
    
    else:
        print(info, "error ChemicalPotential does not support '%s' distribution!" % (dist))
        exit()
    
    if show is not None and show:
        print(info, "%s: get mu with N = %i at T = %.3f nK = %.3f Tc, Tc = %.3f nK, ΔE = %.3f" % (dist, N, T, T/Tc, Tc, dE))

    # in a loop adapt mu until error_N and error_mu reached
    above       = True
    loops_mu    = 0
    loops_dE    = 0
    loops_Emax  = 0
    mu_0_count  = 0
    k = 1.0/(2.0*(omega_bar*k_nK)**3)
    t_start = datetime.now()
    while True:
        
        # 1. find Emax such that Emax^2*dist_func(Emax, T, mu)*dE*k == epsilon_N
        # the step size is dE but not smaller than error_Emax (1nK) in order to have fast calculation.
        # note: this does not affect dE
        Emax_start = Emax
        up = False
        while True:
            p = Emax*Emax*dist_func(Emax, T, mu, Emin)*dE*k
            if p > epsilon_N:
                up = True
                Emax += max(dE, error_Emax)
            elif p < epsilon_N/2:
                if up: break
                else : Emax -= max(dE, error_Emax)
            else: break
            loops_Emax += 1
        if show is not None and show and Emax != Emax_start:
            print(info, '%s: Emax = %.3f nK (start %.3f nK, %i loops, p = %.3e)' % (dist, Emax, Emax_start, loops_Emax, p))
        
        # 2. numerically integrate distribution over energy
        # this runs minimum 2 loops to check if Nint changes when dE is reduced by a factor of 2.
        # note: this might reduce dE. happens rarely, only around T = Tc at beginning of loop.
        loops_dE -= 1
        Nint = 0
        while True:
            x = np.arange(0, Emax, dE)
            Nold = Nint
            if mu == 0.0 and (dist == 'BE' or dist == 'BEC'): 
                # E=0 gives 1/0. excess atoms = BEC fraction
                x=x[1:]
            Nint = np.sum(x*x*dist_func(x, T, mu, Emin))*dE*k
            tact = (datetime.now()-t_start).total_seconds()
            if tact >= tmax:
                # happens sometimes close to critical temperature
                print(info, '%s: chemical potential does not converge after %.1f seconds (%i/%i mu/Emax loops)!' % (dist, tact, loops_mu, loops_Emax))
                print([mu, dE, Emax, T, N, Nint, np.abs(Nint - N), Nold, np.abs(Nint - Nold), error_N, above])
                exit()
            elif np.abs(Nint - Nold)/Nint <= error_Nint: 
                dE *= 2
                break
            dE /= 2
            loops_dE += 1

        if show is not None and show:
            print(info, "%s: %i/%i mu = %+8.1f: ΔN/N = %.3e, δN/N = %.3e, steps = %8i, ΔE = %8.3f" % (dist, loops_mu, loops_dE, mu, np.abs(Nint - N)/N, np.abs(Nint - Nold)/Nint, len(x), dE))

        # 3. abort if error_N and error_mu reached
        # when mu = 0 only error_mu must be fulfilled.
        if (dE <= error_mu) and (mu == 0.0 or np.abs(Nint - N) <= error_N): break
            
        # 4. adapt mu
        # this increments/decrements mu by dE 
        # when direction changes or when mu=0 reduces dE
        if ( Nint > N ):
            if not above:
                above = True
                dE *= kdE
            mu -= dE;
        else:
            if above:
                dE *= kdE
            mu += dE
            if mu > limit: 
                # limit mu and recalculate Emax and Nint=Nthermal and reduce dE. 
                mu = 0.0
                if not above:
                    dE *= kdE                
                if mu_0_count == 0:
                    print(info, '%s: ChemicalPotential warning: mu limited to %.3f' % (dist, limit))
                mu_0_count += 1
            above = False
        loops_mu += 1

    if dist == 'FD':
        Rx = np.sqrt(2*Tc/(k_Epot*mass))/omega  # um. [kB*T = m*omega^2*x^2/2]
        print(info, '%s: T/TF = %.3f nK / %.3f nK = %.3f, N = %.3f' % (dist, T, Tc, T/Tc, N))
    elif dist == 'BE' or dist == 'BEC':
        Rx = np.sqrt(2*Tc/(k_Epot*mass))/omega  # um. [kB*T = m*omega^2*x^2/2]
        if True: #mu == 0.0: # BEC
            # T <= Tc
            Nc = N - Nint # condensed atoms >= 0
            frac = max(1-(T/Tc)**3,0) # expected condensed fraction, Zwierlein M.P.U., Equ. 30
            # check if BEC fraction corresponds to T/Tc within 10% error
            ok = np.abs(Nc/N-frac) < error_BECfrac
            if show is not None or not ok:
                show = None # stop further output
                print(info, '%s: mu = %.3f nK +/- %.3e nK, Nc/N = %.1f/%.1f = %.3f, Emax = %.3f (%i/%i mu/Emax loops, %.1f ms)' % (dist, mu, dE, N - Nint, Nint, (N - Nint)/N, Emax, loops_mu, loops_Emax, (datetime.now()-t_start).total_seconds()*1000))
                print(info, '%s: T/Tc = %.3f, Nc/N = %.3f, expected %.3f, error %.1e (%s)' % (dist, T/Tc, Nc/N, frac, np.abs(Nc/N-frac), 'ok' if ok else 'error'))
            if not ok: exit()
            if mu == 0 and Nint > N:
                print(info, "%s error: for mu = 0 unexpected negative BEC fraction! Nth > Tc" % (dist))
                exit()
            elif mu == 0 and T > Tc:
                print(info, "%s error: for mu = 0 unexpected T > Tc! T/Tc = %.3f nK / %.3f nK = %.3f" % (dist, T, Tc, T/Tc))
                exit()
            elif dist == 'BEC' and abs(mu) > 1e-10:
                print(info, "error: mu = %.3e != 0 for BEC!?" % (mu))
                exit()

    if show is not None:
        # generic result unless BEC
        print(info, '%s: mu = %.3f nK +/- %.3e nK, N error %.3e, Emax = %.3f (%i/%i/%i loops, %.1f ms)' % (dist, mu, dE, np.abs(Nint - N), Emax, loops_mu, loops_dE, loops_Emax, (datetime.now()-t_start).total_seconds()*1000))
            
    return mu, Tc, Nint, Rth, Rc, Emax
    
def Etot(info, dist, N, T, mu, Emax, Emin, bins):
    # returns total energy distribution in harmonic trap
    x = np.linspace(0, Emax, bins)
    y = x*x*dist_func[dist](x, T, mu, Emin)
    dx = x[1]-x[0]
    y_norm = (dx*np.sum(y)*num_bins)/(Emax*N)
    print(info, 'Etot max distribution = %.3f' % (dx*np.sum(y)))
    print(info, 'Etot discretization   = %.3f' % (Emax/num_bins))
    print(info, 'Etot y_norm           = %.3f' % y_norm)
    return [x, y/y_norm]
            
def profiles(info, dist, velocity, N, T, mu, mass, omega, xmax, bins, points):
    # generate [x,y,z] or [vx,vy,vz] profiles velocity == False/True
    # this requires mpmath.polylog and is a bit slow
    print(info, 'calculating', dist, 'profiles ...')

    prof = [None, None, None]
    
    if dist == 'FermiDirac' or dist == 'BoseEinstein(thermal)':
        # FD has - sign in front of exp while BEC has positive sign. rest is the same.
        # note: this is inverse sign of energy distribution functions!
        sign = -1 if dist == 'FermiDirac' else +1
        # cloud size in um. note: this is not Sigma! which is without factor 2!
        R = np.sqrt(2.0*T/(k_Epot*mass))/omega

        if velocity: 
            # map to velocity space
            # note: this collapses all axes on top of each other!
            R *= omega*np.sqrt(k_Epot/k_Ekin)

        for i in range(3):
            x = np.linspace(-xmax[i], xmax[i], points) # um or mm/s
            dx = x[1]-x[0]

            yi = np.array([float(re(sign*polylog(5/2,sign*np.exp(mu/T-(xi/R[i])**2)))) for xi in x])
            area = R[i]*np.sqrt(np.pi)*float(re(sign*polylog(3, sign*np.exp(mu/T))))
            
            prof[i] = [x, yi*N*dx/area*points/bins]        
    
    elif dist == 'BoseEinstein(condensed)':

        if np.abs(mu) <= 1e-3:
            # non-interacting BEC (harmonic oscillator ground state): 2x integration Equ. 38 Ketterle, MPU
            # oscillator length in mu.
            Rc = np.sqrt(1.0/(mass*omega))*k_Rc 

            if velocity: 
                # map to velocity space
                # note: this does not collapse all axes on top of each other
                Rc *= omega*np.sqrt(k_Epot/k_Ekin)

            for i in range(3):
                x = np.linspace(-xmax[i], xmax[i], points) # um or mm/s
                dx = x[1]-x[0]
                        
                # theoretical profile along one spatial coordinate (integrated other 2 coordinates)
                y1d = np.exp(-(x/Rc[i])**2) 
                
                prof[i] = [x, y1d*N/(Rc[i]*np.sqrt(np.pi))*dx*points/bins]
        else:
            # interacting BEC (Thomas-Fermi): 2x integration Equ. 40 Ketterle, MPU
            # interacting size in mu.
            # note: to calculate mu from 'a' see Equ. 50 Ketterle, MPU
            Rc = np.sqrt(2.0*np.abs(mu)/(k_Epot*mass))/omega

            if velocity: 
                # map to velocity space
                # note: this does not collapse all axes on top of each other
                Rc *= omega*np.sqrt(k_Epot/k_Ekin)

            for i in range(3):
                x = np.linspace(-xmax[i], xmax[i], points) # um or mm/s
                dx = x[1]-x[0]
                        
                # theoretical profile along one spatial coordinate (integrated other 2 coordinates)
                y1d = np.max(1.0-(x/Rc[i])**2, 0)**(3.0/2.0) 
                
                prof[i] = [x, y1d*10*N/(3*np.sqrt(np.pi)**(3.0/2.0)*Rc[i])*dx*points/bins]
    
    elif dist == 'MaxwellBoltzmann':
        # cloud size in um
        sigma = np.sqrt(T/(k_Epot*mass))/omega

        if velocity: 
            # map to velocity space
            # note: this collapses all axes on top of each other!
            sigma *= omega*np.sqrt(k_Epot/k_Ekin)

        # 3d peak density in atoms/mu^3
        n3d = N/((2*np.pi)**(3/2)*sigma[0]*sigma[1]*sigma[2])

        for i in range(3):
            x = np.linspace(-xmax[i], xmax[i], points) # um or mm/s
            dx = x[1]-x[0]
                    
            # theoretical profile along one spatial coordinate (integrated other 2 coordinates)
            y1d = np.exp(-0.5*(x/sigma[i])**2)
            
            prof[i] = [x, y1d*n3d*2.0*np.pi*sigma[(i+1)%3]*sigma[(i+2)%3]*dx*points/bins]

    else:
        print(info, 'profile for ', dist, 'not implemented!')
        exit()

    return prof

if __name__ == '__main__':

    # load parameter file
    print('loading parameter file ', param_file, ' ...')
    params = load_params(param_file, skip=['//', '#'], max_lines=None)
    print(len(params), 'parameters loaded\n')
    
    # load specific parameters
    # might throw KeyError for non-existent parameters.
    num_species      = int(params['num_species'     ][0])
    repeat           = int(params['repeat'          ][0])
    num_bins         = int(params['num_bins'        ][0])
    species          = [   params['species_0'       ][0]]
    statistics       = [   params['statistics_0'    ][0]]
    mass             = [   params['mass_0'          ][0]]
    N                = [   params['N_0'             ][0]]
    T                = [   params['T_0'             ][0]]
    fx               = [   params['fx_0'            ][0]]
    fy               = [   params['fy_0'            ][0]]
    fz               = [   params['fz_0'            ][0]]
    if num_species  == 2:
        species     += [   params['species_1'       ][0]] 
        statistics  += [   params['statistics_1'    ][0]]
        mass        += [   params['mass_1'          ][0]]
        N           += [   params['N_1'             ][0]]
        T           += [   params['T_1'             ][0]]
        fx          += [   params['fx_1'            ][0]]
        fy          += [   params['fy_1'            ][0]]
        fz          += [   params['fz_1'            ][0]]
    mol_name         = ''.join(species)
    mol_statistics   = [   params['statistics_mol'  ][0]]
    distance_measure =     params['distance_measure'][0]
    gamma            =     params['gamma'           ][0]
    
    if repetition <= 0 or repetition > repeat:
        print('repetition ', repetition, ' is outside of range 1 ..', repeat)
        exit()
        
    # calculate expected distribution
    omega     = [None for _ in range(num_species)]
    omega_bar = [None for _ in range(num_species)]
    mu        = [None for _ in range(num_species)]
    Emax_theo = [None for _ in range(num_species)]
    Tcrit     = [None for _ in range(num_species)]
    TTcrit    = [None for _ in range(num_species)]
    Etot_theo = [None for _ in range(num_species)]
    Rth       = [None for _ in range(num_species)]
    Rc        = [None for _ in range(num_species)]
    for i in range(num_species):
        omega          [i] = 2*np.pi*np.array([fx[i], fy[i], fz[i]])
        omega_bar      [i] = 2*np.pi*(fx[i]*fy[i]*fz[i])**(1.0/3.0)
        Tcrit          [i] = Tcrit_func[statistics[i]](N[i], omega_bar[i])
        TTcrit         [i] = T[i]/Tcrit[i] if Tcrit[i] != 0.0 else 0.0

        mu[i], Tcrit[i], _, Rth[i], Rc[i], Emax_theo[i] = ChemicalPotential(species[i], statistics[i], N[i], T[i], omega[i], mass[i],show=False)

        print()

    try:
        # load histogram files
        hist_file = params['histogram_file'][0]
    except KeyError:
        print('no histogram file in parameters!')
        hist_file = None
        
    if hist_file is not None:
        split = hist_file.split('.')
        if len(split) < 2:
            print('histogram file without file extension? ', hist_file)
            exit()
        hist_file = '.'.join(split[:-1])
        hist_ext  = split[-1]
        
        species_list = []
        for sp in range(num_species):
            if statistics[sp] == 'BoseEinstein':
                Nc = int(np.round(N[sp]*(1.0-(TTcrit[sp])**3)))
                if Nc <= 0:
                    species_list.append([sp, species[sp]+'(thermal)', N[sp], statistics[sp]+'(thermal)'])
                else:
                    if N[sp] - Nc > 0:
                        species_list.append([sp, species[sp]+'(thermal)', N[sp]-Nc, statistics[sp]+'(thermal)'])
                    if Nc > 0:
                        species_list.append([sp, species[sp]+'(condensed)', Nc, statistics[sp]+'(condensed)'])
            else:
                species_list.append([sp, species[sp], N[sp], statistics[sp]])
        
        hist_data = [{} for _ in range(num_species+1)]
        for sp, species_name, species_N, species_stat in species_list:
            if fig_species is None or species[sp] in fig_species:

                # header contains max values used to scale data
                f = hist_file + '_' + species_name + '_' + str(repetition) + '.' + hist_ext
                print(species_name, 'loading histogram file header', f, ' ...')
                header = load_params(f, skip=[], max_lines=hist_header_len)

                stat  =     header['// stat'    ][0]            
                rep   = int(header['// rep'     ][0])
                N_    = int(header['// N'       ][0])
                T_    =     header['// T'       ][0]
                mu_   =     header['// mu'      ][0]
                xmax  =     header['// max x'   ][0]
                ymax  =     header['// max y'   ][0]
                zmax  =     header['// max z'   ][0]
                vxmax =     header['// max vx'  ][0]
                vymax =     header['// max vy'  ][0]
                vzmax =     header['// max vz'  ][0]
                Emax  =     header['// max Etot'][0]

                if i < num_species and (                                \
                         stat    != species_stat     or \
                         rep     != repetition       or \
                         N_      != species_N        or \
                   round(T_ , 3) != round(T [sp], 3) or \
                   round(mu_, 2) != round(mu[sp], 2) 
                   ):
                    print(species_name, 'histogram file inconsistent with parameter file!')
                    print('statistics:', stat, "vs.", species_stat)
                    print('repetition:', rep , "vs.", repetition)
                    print('N         :', N_  , "vs.", species_N)
                    print('T         :', T_  , "vs.", T[sp])
                    print('mu        :', mu_ , "vs.", mu[sp])
                    exit()
                
                for h in hist_header:
                    key,unit = h[:2]
                    if key is not None and header[key][1] != unit:
                        print(species_name, 'histogram file unexpected unit for entry %s:' % (key), "'%s' != '%s'" % (header[key][1], unit))
                        exit()

                print(species_name, 'histogram file header ok\n')

                # load data from histogram file
                # returned lines correspond to columns in histogram
                print(species_name, 'loading histogram file data', f, ' ...')
                hist = np.loadtxt(f, unpack=True, delimiter=' ', skiprows=hist_header_len, 
                    usecols    = (i for i in range(hist_num_cols)), 
                    converters = {i:get_int for i in range(hist_num_cols)})
                if len(hist[0]) != num_bins:
                    print(species_name, 'number of bins', len(hist[0]), '!=', num_bins)
                    print(hist_header_len)
                    exit()
                elif len(hist) != hist_num_cols:
                    print(species_name, 'number of columns', len(hist), '!=', hist_num_cols)
                    exit()
                print(species_name, "%i bins x %i cols loaded\n" % (len(hist[0]), len(hist)))
                    
                # save each data column separately with scaling as defined in header
                # position and momentum are symmetric, energy is single-sided.
                data         = {label:[] for label in hist_figure if label is not None}
                data_labels  = {label:[] for label in hist_figure if label is not None}
                data_args    = {label:[] for label in hist_figure if label is not None}
                curves       = {label:[] for label in hist_figure if label is not None}
                curve_labels = {label:[] for label in hist_figure if label is not None}
                curve_args   = {label:[] for label in hist_figure if label is not None}
                for j,col in enumerate(hist_colum_index):
                    if col is not None and hist_header_entry[j] is not None:
                        label = hist_figure[j] # figure
                        max_  = header[hist_header_entry[j]][0] # maximum value
                        rng   = [0, max_] if hist_units[j] == 'nK' else [-max_, max_]
                        xj = np.linspace(rng[0], rng[1], num_bins+1)
                        xj = (xj[1:] + xj[0:-1])/2 # center of bins
                        yj = hist[col] # histogram counts
                        hist_data[sp][hist_colum_name[j]]  = [xj, yj] # data for each species and column
                        if label is not None: # save plot data
                            data        [label] += [[xj, yj]]
                            data_labels [label] += [species_name + ' ' + hist_colum_name[j]]
                            data_args   [label] += [{'color':hist_colors[j]}]
                            
                # add theory curves
                if fig_sel is None or fig_E in fig_sel:
                    Emin = 0.5*omega_bar[sp]*k_nK
                    curves      [fig_E] += [Etot(species_name, species_stat, species_N, T[sp], mu[sp], Emax, Emin, num_bins)]
                    curve_labels[fig_E] += [species_name + ' Etot']
                    curve_args  [fig_E] += [{'color': 'Red'}]
                
                if fig_sel is None or fig_x in fig_sel:
                    curves      [fig_x] += profiles(species_name, species_stat, False, species_N, T[sp], mu[sp], mass[sp], omega[sp], [xmax, ymax, zmax], num_bins, points=100)
                    curve_labels[fig_x] += [species_name + ' x', species_name + ' y', species_name + ' z']
                    curve_args  [fig_x] += [{'color': 'Red'}, {'color': 'Blue'}, {'color': 'Green'}]
                    
                    if True: # add thermal gas profile for comparison
                        name = species[sp] + '(MB)'
                        curves      [fig_x] += profiles(name, 'MaxwellBoltzmann', False, species_N, T[sp], mu[sp], mass[sp], omega[sp], [xmax, ymax, zmax], num_bins, points=100)
                        curve_labels[fig_x] += [name + ' x', name + ' y', name + ' z']
                        curve_args  [fig_x] += [{'color': 'Red'  , 'linestyle':'dashed'}, 
                                                {'color': 'Blue' , 'linestyle':'dashed'}, 
                                                {'color': 'Green', 'linestyle':'dashed'}]
                    
                if fig_sel is None or fig_v in fig_sel:
                    curves      [fig_v] += profiles(species_name, species_stat, True, species_N, T[sp], mu[sp], mass[sp], omega[sp], [vxmax, vymax, vzmax], num_bins, points=100)
                    curve_labels[fig_v] += [species_name + ' vx', species_name + ' vy', species_name + ' vz']
                    curve_args  [fig_v] += [{'color': 'Red'}, {'color': 'Blue'}, {'color': 'Green'}]
                    
                    if True: # add thermal gas profile for comparison
                        name = species[sp] + '(MB)'
                        curves      [fig_v] += profiles(name, 'MaxwellBoltzmann', True, species_N, T[sp], mu[sp], mass[sp], omega[sp], [vxmax, vymax, vzmax], num_bins, points=100)
                        curve_labels[fig_v] += [name + ' vx', name + ' vy', name + ' vz']
                        curve_args  [fig_v] += [{'color': 'Red'  , 'linestyle':'dashed'}, 
                                                {'color': 'Blue' , 'linestyle':'dashed'}, 
                                                {'color': 'Green', 'linestyle':'dashed'}]
                
                # plot data for each species and selected label
                for label in data.keys():
                    if label is not None and (fig_sel is None or label in fig_sel):
                        plot(
                            title        = species[sp] + ' histogram ' + label, 
                            data         = data        [label], 
                            data_labels  = data_labels [label],
                            data_args    = data_args   [label], 
                            curves       = curves      [label],   
                            curve_labels = curve_labels[label],
                            curve_args   = curve_args  [label], 
                            xlabel       = label,
                            x_range      = None,
                            ylabel       = 'count',
                            y_range      = None,
                            label_pos    = 'upper right', 
                            label_cols   = 1,
                            fig_size     = (10*4/3,7), # figure (width,height)
                            fig_pos      = [0.07, 0.08, 0.87, 0.87] # sub plot (left,bottom,width,height)
                            ) 
                
        plt.show()
            

