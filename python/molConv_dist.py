#!/usr/local/bin/python3

# plots distribution in histograms and atoms/molecule files
# run molConv_run.py and export histograms and atoms/molecules files
# give here parameter file name and repetition number (>=1) to be analyzed
# plots distribution of all atoms and molecule obtained from histogram files or/and atoms files
# and compares with expected distribution using parameter file.
# if load_atoms == True generates additionally histogram directly from atoms/molecule files.
# last changed 23/9/2026 by Andi

# TODO: at the moment reads only histogram files and does not load atoms csv files

param_file  = './tmp/LiCr_test_delta_x*delta_p_3.8e-01/LiCr_test_params_0.txt'
repetition  = 1
load_atoms  = False # has no effect at the moment!

# statistics
stat_FD = 'FermiDirac'
stat_BE = 'BoseEinstein'
stat_MB = 'MaxwellBoltzmann'

# sub species identifier for Bose gas
sub_thermal = '(thermal)'
sub_BEC     = '(condensed)'

# figure identifier
fig_x = 'position'
fig_v = 'velocity'
fig_E = 'energy'

# figures and species to be plotted. None = all, [] = none
fig_sel     = [] #[fig_x, fig_v, fig_E]
fig_species = ['Li6', 'Cr52'] #['Li6', 'Cr52', 'Li6Cr52(thermal)', 'Li6Cr52(condensed)']

# points for theory curves. keep not too big since polylog is a bit slow
points = 100

# if not None show classic thermal gas (MaxwellBoltzmann) for comparison with this label added
show_thermal = '(th. gas)'

# if True make a nice overview figure for selected species
show_panels = True

# select which panels to show (top to bottom)
panels = [fig_E, fig_x, fig_v]

# if not None scale panels to data range with this relative margin 
# when not None plots panel with y-axis in log scale, otherwise not
scale_to_data = 0.05 # 0.05

# panel customization
def panel_adjust(species_name):
    if species_name == 'Li6':
        title     = species_name + ' (Fermion) histogram'
        xrange    = {fig_E:[0.5, 6e2], fig_x:[-330,330], fig_v:[-30, 30]}
        yrange    = {fig_E:[0.5, 5e5], fig_x:[0.5, 5e5], fig_v:[0.5, 5e5]}
        label_pos = {fig_E:['upper right',(0.97, 0.97), 2],
                     fig_x:['upper right',(0.97, 0.97), 3],
                     fig_v:['upper right',(0.97, 0.97), 3]}
    elif species_base_name == 'Cr52':
        title     = species_name + ' (Boson) histogram'
        xrange    = {fig_E:[0.5, 6e2], fig_x:[-110,110], fig_v:[-10,10]}
        yrange    = {fig_E:[0.5, 5e5], fig_x:[0.5, 5e5], fig_v:[0.5, 5e5]}
        label_pos = {fig_E:['upper right', (0.97, 0.97), 3],
                     fig_x:['upper right', (0.97, 0.97), 5],
                     fig_v:['upper right', (0.97, 0.97), 5]}
    else:
        title     = species_name + ' histogram'
        xrange    = {fig:None for fig in panels}
        yrange    = {fig:None for fig in panels}
        label_pos = {fig:['upper left',(1.01, 1.04)] for fig in panels}
    return [title, xrange, yrange, label_pos]
    
# displayed sub-species names in panel labels
def panel_label(species_name, sub_name):
    if sub_name == species_name + sub_thermal:
        label = 'BEC (th.)'
    elif sub_name == species_name + sub_BEC:
        label = 'BEC (cond.)'
    elif sub_name == species_name + show_thermal:
        label = 'thermal gas'
    elif sub_name == species_name:
        label = species_name
    else:
        print('unrecognized sub-species?', sub_name)
        exit()
    return label


from molConv_run import plot
import numpy as np
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
hist_column_index   = [h[2] for h in hist_header]
hist_column_name    = [h[3] for h in hist_header]
hist_colors         = [h[4] for h in hist_header]
hist_figure         = [h[5] for h in hist_header]
hist_header_len     = len(set([h for h in hist_header_entry if h is not None])) + 3
hist_num_cols       = len([h for h in hist_column_index if h is not None])

# convert text units to Latex labels
unit_to_label = {'nK': r'nK', 'um': r'$\mu$m', 'mm/s': r'mm/s'}

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

dist_func  = {stat_FD             : FermiDirac, 
              stat_BE             : BoseEinstein_thermal, 
              stat_BE+sub_thermal : BoseEinstein_thermal,
              stat_BE+sub_BEC     : BoseEinstein_condensed, 
              stat_MB             : MaxwellBoltzmann       }
Tcrit_func = {stat_FD             : TFermi, 
              stat_BE             : Tcrit_, 
              stat_MB             : lambda N,omega_bar: 0.0}
    
def ChemicalPotential(info, dist, N, T, omega, mass, error_mu=1e-3, error_N=0.1, error_Emax = 1.0, epsilon_N=None, error_Nint=0.01, error_BECfrac=0.01, kdE=0.8, tmax=5.0, show=False):
    # get chemical potential such that integrating distribution gives atom number N
    # info          = species name used for printing
    # dist          = stat_FD, stat_BE or stat_MB
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
    
    if dist == stat_MB: # Maxwell-Boltzmann: mu = Tc = 0.0
        Emax = T*np.log(N/epsilon_N)
        if show is not None:
            print(info, '%s: mu = 0.0, Emax = %.3f' % (dist, Emax))
        return 0.0, 0.0, N, Rth, 0.0, Emax
    elif dist == stat_FD: # Fermi-Dirac: Tc = EF
        dist_func = FermiDirac
        Tc = TFermi(N, omega_bar); # nK
        Rc = np.sqrt(2.0*T/(k_Epot*mass))/omega  # Radius in Polylog in um. [kB*T = m*omega^2*x^2/2]
        mu = Tc # start with Fermi energy
        #print('mu start %.3f' % mu)
        dE = Tc # initial energy steps
        limit = Tc # limit mu <= TF. Natoms should be still reached.
        Emax = Tc+T
    elif dist == stat_BE: # Bose-Einstein: Tc, mu <= 0
    
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
            if mu == 0.0 and (dist == stat_BE): 
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

    if dist == stat_FD:
        Rx = np.sqrt(2*Tc/(k_Epot*mass))/omega  # um. [kB*T = m*omega^2*x^2/2]
        print(info, '%s: T/TF = %.3f nK / %.3f nK = %.3f, N = %.3f' % (dist, T, Tc, T/Tc, N))
    elif dist == stat_BE:
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
                print(info, '%s: T/Tc = %.3f / %.3f = %.3f, Nc/N = %.3f, expected %.3f, error %.1e (%s)' % (dist, T, Tc, T/Tc, Nc/N, frac, np.abs(Nc/N-frac), 'ok' if ok else 'error'))
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
    if isinstance(bins, (list, tuple, np.ndarray)):
        x = np.array(bins)
    else:
        x = np.linspace(0, Emax, bins)
    y = x*x*dist_func[dist](x, T, mu, Emin)
    dx = x[1]-x[0]
    y_norm = (dx*np.sum(y)*num_bins)/(Emax*N)
    print(info, 'Etot max distribution = %.3f' % (dx*np.sum(y)))
    print(info, 'Etot discretization   = %.3f' % (Emax/num_bins))
    print(info, 'Etot y_norm           = %.3f' % y_norm)
    return [x, y/y_norm]
            
def profiles(info, dist, velocity, N, T, mu, mass, omega, xmax, dx, bins, points):
    # generate [x,y,z] or [vx,vy,vz] profiles velocity == False/True
    # this requires mpmath.polylog and is a bit slow
    # xmax = +/-maximum range used for plotting
    # dx   = step size between bins for each coordinate
    # bins = number of bins used to obtain histogram. this is used to scale counts to match the histogram.
    # points = number of points to plot. calculation is a bit slow, so do not use too many points.
    print(info, 'calculating', dist, 'profiles ...')

    prof = [None, None, None]
    
    if dist == stat_FD or dist == stat_BE+sub_thermal:
        # FD has - sign in front of exp while BEC has positive sign. rest is the same.
        # note: this is inverse sign of energy distribution functions!
        sign = -1 if dist == stat_FD else +1
        # cloud size in um. note: this is not Sigma! which is without factor 2!
        R = np.sqrt(2.0*T/(k_Epot*mass))/omega

        if velocity: 
            # map to velocity space
            # note: this collapses all axes on top of each other!
            R *= omega*np.sqrt(k_Epot/k_Ekin)

        for i in range(3):
            x = np.linspace(-xmax[i], xmax[i], points) # um or mm/s
            #dx = x[1]-x[0]
            
            yi = np.array([float(re(sign*polylog(5/2,sign*np.exp(mu/T-(xi/R[i])**2)))) for xi in x])
            area = R[i]*np.sqrt(np.pi)*float(re(sign*polylog(3, sign*np.exp(mu/T))))
            
            prof[i] = [x, yi*N*dx[i]/area*points/bins]        
    
    elif dist == stat_BE+sub_BEC:

        if np.abs(mu) <= 1e-3:
            # non-interacting BEC (harmonic oscillator ground state): 2x integration Equ. 38 Ketterle, MPU
            # oscillator length in mu.
            Rc = np.sqrt(1.0/(mass*omega))*k_Rc 

            if velocity: 
                # map to velocity space
                # note: this does NOT collapse all axes on top of each other!
                Rc *= omega*np.sqrt(k_Epot/k_Ekin)

            for i in range(3):
                x = np.linspace(-xmax[i], xmax[i], points) # um or mm/s
                #dx = x[1]-x[0]
                        
                # theoretical profile along one spatial coordinate (integrated other 2 coordinates)
                y1d = np.exp(-(x/Rc[i])**2) 
                
                prof[i] = [x, y1d*N/(Rc[i]*np.sqrt(np.pi))*dx[i]*points/bins]
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
                #dx = x[1]-x[0]
                        
                # theoretical profile along one spatial coordinate (integrated other 2 coordinates)
                y1d = np.max(1.0-(x/Rc[i])**2, 0)**(3.0/2.0) 
                
                prof[i] = [x, y1d*10*N/(3*np.sqrt(np.pi)**(3.0/2.0)*Rc[i])*dx[i]*points/bins]
    
    elif dist == stat_MB:
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
            #dx = x[1]-x[0]
                    
            # theoretical profile along one spatial coordinate (integrated other 2 coordinates)
            y1d = np.exp(-0.5*(x/sigma[i])**2)
            
            prof[i] = [x, y1d*n3d*2.0*np.pi*sigma[(i+1)%3]*sigma[(i+2)%3]*dx[i]*points/bins]

    else:
        print(info, 'profile for', dist, 'not implemented!')
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
    mol_statistics   =     params['statistics_mol'  ][0] 
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
    
        # get list of species  
        # list contains [species index, species name, atom/molecule number, statistics, linked species name for MB or None]
        atoms_list = []
        MB_list    = []
        for sp in range(num_species):
            if statistics[sp] == stat_BE:
                Nc = int(np.round(N[sp]*(1.0-(TTcrit[sp])**3)))
                if Nc <= 0:
                    name = sub_thermal
                    atoms_list .append([sp, species[sp]+name, N[sp], statistics[sp]+name, None])
                    link = species[sp]+name
                    #Nth  = N[sp]
                else:
                    if N[sp] - Nc > 0:
                        name = sub_thermal
                        atoms_list.append([sp, species[sp]+name, N[sp]-Nc, statistics[sp]+name, None])
                        link = species[sp]+name
                        #Nth  = N[sp]
                    if Nc > 0:
                        name = sub_BEC
                        atoms_list.append([sp, species[sp]+name, Nc, statistics[sp]+name, None])
                        if N[sp] - Nc <= 0:
                            link = species[sp]+name
                            #Nth  = Nc
            else:
                atoms_list.append([sp, species[sp], N[sp], statistics[sp], None])
                link = species[sp]
                #Nth  = N[sp]
            if show_thermal:
                MB_list.append([sp, species[sp]+show_thermal, N[sp], stat_MB, link])
        
        # add list of molecules = combinations of species
        mol_list = []
        for sp0 in atoms_list:
            for sp1 in atoms_list:
                if sp0[0] < sp1[0]: 
                    mol_list.append([[sp0[0], sp1[0]], sp0[1]+sp1[1], 0, mol_statistics, None])

        print('atoms    :', [sp [1] for sp  in atoms_list])
        print('molecules:', [mol[1] for mol in mol_list  ])
        print()
        
        species_list = atoms_list + mol_list + MB_list
        hist_data   = {sp[1]:{} for sp in species_list} # species data for each hist column name
        hist_curves = {sp[1]:{} for sp in species_list} # species theory curves for each hist column name (if available)
        for species_list_entry in species_list:
            sp, species_name, species_N, species_stat, species_link = species_list_entry # entry allows to insert species_N into mol_list
            
            is_molecule = isinstance(sp, (list, tuple))
            species_base_name = species_name if is_molecule else species[sp]
            
            if fig_species is None or species_base_name in fig_species:
            
                # histogram available
                # header contains max values used to scale data
                f = hist_file + '_' + (species_name if species_link is None else species_link) + '_' + str(repetition) + '.' + hist_ext
                print(species_name, 'loading histogram file header', f, ' ...')
                header = load_params(f, skip=[], max_lines=hist_header_len)

                species_stat_ =     header['// stat'    ][0]            
                species_rep   = int(header['// rep'     ][0])
                species_N_    = int(header['// N'       ][0])
                species_T     =     header['// T'       ][0]
                species_mu    =     header['// mu'      ][0]
                species_xmax  =     header['// max x'   ][0]
                species_ymax  =     header['// max y'   ][0]
                species_zmax  =     header['// max z'   ][0]
                species_vxmax =     header['// max vx'  ][0]
                species_vymax =     header['// max vy'  ][0]
                species_vzmax =     header['// max vz'  ][0]
                species_Emax  =     header['// max Etot'][0]

                if species_link is None:
                    if species_stat_ != species_stat     or \
                       species_rep   != repetition       :
                        print(species_name, 'histogram file inconsistent with parameter file!')
                        print('statistics:', species_stat_, "vs.", species_stat)
                        print('repetition:', species_rep  , "vs.", repetition)
                        exit()

                    if is_molecule:
                        # we get from histogram N, T, and mu and mass, omega, omega_bar we can calculate
                        # Epot = m0*w0^2*x^2/2+m1*w1^2*x^2/2 = (m0*w0^2+m1*w1^2)*x^2/2 =(m0+m1)*<w>^2*x^2 
                        # -> <w> = sqrt((m0*w0^2+m1*w1^2)/(m0+m1))
                        if species_N_ == 0:
                            print(species_name, 'zero molecules found (skip)')
                            continue
                        species_N         = species_N_
                        species_mass      = mass[sp[0]] + mass[sp[1]]
                        species_omega     = np.sqrt((mass[sp[0]]*omega[sp[0]]**2+mass[sp[1]]*omega[sp[1]]**2)/species_mass)
                        species_omega_bar = (species_omega[0]*species_omega[1]*species_omega[2])**(1.0/3.0)
                        # insert species_N back into mol_list
                        species_list_entry[2] = species_N
                    else:
                        if species_N_           != species_N        or \
                           round(species_T , 3) != round(T [sp], 3) or \
                           round(species_mu, 2) != round(mu[sp], 2) :
                            print(species_name, 'histogram file (atoms) inconsistent with parameter file!')
                            print('N         :', species_N_ , "vs.", species_N)
                            print('T         :', species_T  , "vs.", T[sp])
                            print('mu        :', species_mu , "vs.", mu[sp])
                            exit()
                        species_mass      = mass     [sp]
                        species_omega     = omega    [sp]
                        species_omega_bar = omega_bar[sp]
                else:
                    # linked species might have different N and statistics
                    if species_rep          != repetition       or \
                       species_N            != N[sp]            or \
                       round(species_T , 3) != round(T [sp], 3) or \
                       round(species_mu, 2) != round(mu[sp], 2) :
                        print(species_name, 'histogram file (atoms) inconsistent with parameter file!')
                        print('N         :', species_N  , "vs.", N[sp])
                        print('rep       :', species_rep, "vs.", repetition)
                        print('T         :', species_T  , "vs.", T[sp])
                        print('mu        :', species_mu , "vs.", mu[sp])
                        exit()
                    species_mass      = mass     [sp]
                    species_omega     = omega    [sp]
                    species_omega_bar = omega_bar[sp]
                                
                # check units
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
                data         = {label:[]   for label in hist_figure if label is not None}
                data_labels  = {label:[]   for label in hist_figure if label is not None}
                data_args    = {label:[]   for label in hist_figure if label is not None}
                data_unit    = {label:None for label in hist_figure if label is not None}
                curves       = {label:[]   for label in hist_figure if label is not None}
                curve_labels = {label:[]   for label in hist_figure if label is not None}
                curve_args   = {label:[]   for label in hist_figure if label is not None}
                for j,col in enumerate(hist_column_index):
                    if col is not None and hist_header_entry[j] is not None:
                        label = hist_figure[j] # figure
                        max_  = header[hist_header_entry[j]][0] # maximum value
                        rng   = [0, max_] if hist_units[j] == 'nK' else [-max_, max_]
                        if data_unit[label] is None:
                            data_unit[label] = hist_units[j]
                        elif data_unit[label] != hist_units[j]:
                            print(species_name, label, "inconsistent unit '%s' != '%s'" % (data_unit[label], hist_units[j]))
                            exit()
                        xj = np.linspace(rng[0], rng[1], num_bins+1)
                        xj = (xj[1:] + xj[0:-1])/2 # center of bins
                        yj = hist[col] # histogram counts
                        if species_link is None:
                            hist_data[species_name][hist_column_name[j]] = [xj, yj] # data for each species and hist column
                        if label is not None: # save plot data
                            data        [label] += [[xj, yj]]
                            data_labels [label] += [species_name + ' ' + hist_column_name[j]]
                            data_args   [label] += [{'color':hist_colors[j]}]

                # energy theory curve [Etot]
                # TODO: load color from hist_color
                species_Emin = 0.5*species_omega_bar*k_nK
                if scale_to_data is not None: 
                    # scale to data range where count > 0
                    x,y = hist_data[species_name if species_link is None else species_link]['Etot']
                    index = np.where(y > 0)
                    rng = [x[index][0], x[index][-1]]
                    bins = np.linspace(max(rng[0]-(rng[1]-rng[0])*scale_to_data, species_Emin), 
                                           rng[1]+(rng[1]-rng[0])*scale_to_data, num_bins)
                else:
                    bins = num_bins
                curve = Etot(species_name, species_stat, species_N, species_T, 
                             species_mu, species_Emax, species_Emin, bins)
                curves      [fig_E] += [curve]
                curve_labels[fig_E] += [species_name + ' Etot']
                curve_args  [fig_E] += [{'color': 'Red'}]
                hist_curves[species_name]['Etot'] = curve
            
                # position theory curves [x,y,z]
                # TODO: load color from hist_color
                if scale_to_data is not None: 
                    # scale to data range where count > 0
                    xmax = [0,0,0]
                    for i,label in enumerate(['x', 'y', 'z']):
                        x,y = hist_data[species_name if species_link is None else species_link][label]
                        index = np.where(y > 0)
                        if xmax[i] < abs(x[index][-1]): xmax[i] = abs(x[index][-1])
                else:
                    xmax = [species_xmax, species_ymax, species_zmax]    
                curve = profiles(species_name, species_stat, False, species_N, species_T, species_mu, 
                                 species_mass, species_omega, xmax,
                                 [species_xmax*2/(points-1), species_ymax*2/(points-1), species_zmax*2/(points-1)],
                                 num_bins, points=points)
                curves      [fig_x] += curve
                curve_labels[fig_x] += [species_name + ' x', species_name + ' y', species_name + ' z']
                curve_args  [fig_x] += [{'color': 'Red'}, {'color': 'Blue'}, {'color': 'Green'}]
                hist_curves[species_name]['x'] = curve[0]
                hist_curves[species_name]['y'] = curve[1]
                hist_curves[species_name]['z'] = curve[2]
                
                # velocity theory curve
                # TODO: load color from hist_color
                if scale_to_data is not None: 
                    # scale to data range where count > 0
                    xmax = [0,0,0]
                    for i,label in enumerate(['vx', 'vy', 'vz']):
                        x,y = hist_data[species_name if species_link is None else species_link][label]
                        index = np.where(y > 0)
                        if xmax[i] < abs(x[index][-1]): xmax[i] = abs(x[index][-1])
                else:
                    xmax = [species_vxmax, species_vymax, species_vzmax]
                curve = profiles(species_name, species_stat, True, species_N, species_T, species_mu, 
                                 species_mass, species_omega, xmax,
                                 [species_vxmax*2/(points-1), species_vymax*2/(points-1), species_vzmax*2/(points-1)],
                                 num_bins, points=points)
                curves      [fig_v] += curve
                curve_labels[fig_v] += [species_name + ' vx', species_name + ' vy', species_name + ' vz']
                curve_args  [fig_v] += [{'color': 'Red'}, {'color': 'Blue'}, {'color': 'Green'}]
                hist_curves[species_name]['vx'] = curve[0]
                hist_curves[species_name]['vy'] = curve[1]
                hist_curves[species_name]['vz'] = curve[2]
                
                # plot data for each species and selected label
                for label in data.keys():
                    if label is not None and (fig_sel is None or label in fig_sel):
                        plot(
                            title        = species_name + ' ' + label, 
                            data         = data        [label], 
                            data_labels  = data_labels [label],
                            data_args    = data_args   [label], 
                            curves       = curves      [label],   
                            curve_labels = curve_labels[label],
                            curve_args   = curve_args  [label], 
                            xlabel       = label + (' (%s)'%unit_to_label[data_unit[label]]),
                            x_range      = None,
                            ylabel       = 'count',
                            y_range      = None,
                            label_pos    = 'upper right', 
                            label_cols   = 1,
                            fig_size     = (10*4/3,7), # figure (width,height)
                            fig_pos      = [0.07, 0.08, 0.87, 0.87] # sub plot (left,bottom,width,height)
                            ) 
                            
        if show_panels:
            # overview figure for each fig_species with panels
            fig_species_names = [sp[1] for sp in species_list] if fig_species is None else fig_species
            for species_base_name in fig_species_names:

                # find all sub-species names with same base name
                sub_list = []
                for sp, species_name, species_N, species_stat, species_link in species_list:
                    is_molecule = isinstance(sp, (list, tuple))
                    species_base_name_ = species_name if is_molecule else species[sp]
                    if species_base_name_ == species_base_name:
                        sub_list.append(species_name)
                        
                print(species_base_name, sub_list)
                
                ax = None
                for i,fig in enumerate(panels):
                
                    # get list of rows in hist_header for colors and units
                    rows = []
                    for row,f in enumerate(hist_figure):
                        if f == fig:
                            rows.append(row)
                                
                    title, xrange, yrange, label_pos = panel_adjust(species_base_name)    

                    left   = 0.075
                    width  = 1-left-0.01
                    bottom = 0.05
                    height = 1-bottom-0.035
                    gap    = 0.060
                    h      = (height-(len(panels)-1)*gap)/len(panels)
                    b      = bottom + (len(panels)-i-1)*(h+gap)
                            
                    ax = plot( 
                        ax           = ax,
                        title        = title if ax is None else None, 
                        data         = [hist_data[sp][hist_column_name[r]] 
                                        for sp in sub_list for r in rows if hist_column_name[r] in hist_data[sp]], 
                        data_labels  = ['%s %s'%(panel_label(species_base_name, sp), hist_column_name[r]) 
                                        for sp in sub_list for r in rows if hist_column_name[r] in hist_data[sp]],
                        data_args    = [{'color'    :hist_colors[r], 
                                         'facecolor':[hist_colors[r],'White','Gray'][j],
                                         'edgecolor':hist_colors[r]
                                        } for j,sp in enumerate(sub_list) for r in rows if hist_column_name[r] in hist_data[sp]], 
                        curves       = [hist_curves[sp][hist_column_name[r]] 
                                        for sp in sub_list for r in rows if hist_column_name[r] in hist_curves[sp]],   
                        curve_labels = ['%s %s'%(panel_label(species_base_name, sp), hist_column_name[r]) 
                                        for sp in sub_list for r in rows if hist_column_name[r] in hist_curves[sp]],
                        curve_args   = [{'color':hist_colors[r], 
                                         'linestyle':['solid','dashed','dotted'][j],
                                        } for j,sp in enumerate(sub_list) for r in rows if hist_column_name[r] in hist_curves[sp]], 
                        xlabel       = fig + (' (%s)'%unit_to_label[data_unit[fig]]),
                        x_range      = xrange[fig],
                        ylabel       = 'count',
                        y_range      = yrange[fig],
                        log_scale    = [True, True] if fig==fig_E else [False, scale_to_data],
                        label_pos    = label_pos[fig],
                        fig_size     = (10,10), # figure (width,height)
                        fig_pos      = [left, b, width, h] # sub plot (left,bottom,width,height)
                        ) 
                            
        plt.show()
            

