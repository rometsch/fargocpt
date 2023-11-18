#!/usr/bin/env python3
import subprocess
import sys
import os
from argparse import ArgumentParser

import numpy as np
import matplotlib.pyplot as plt


from astropy import units as u
from astropy import constants
import scipy.special
import scipy.integrate
import matplotlib.pyplot as plt

def compile_fargo(fargo_path, Nthreads, silent=False):

    if silent:
        out = open("make.log", "w")
        err = out
    else:
        out = sys.stdout
        err = sys.stderr

    subprocess.run(['make', '-j', str(Nthreads), '-C' 'src/'], cwd=fargo_path, stdout=out, stderr=err)

    if silent: 
        out.close()

def run(fargo_path, par_file, Nthreads, Nprocs, silent=False):
    from fargocpt import run as run_fargo

    if silent:
        out = open("sim.log", "w")
        err = open("sim.err", "w")
    else:
        out = sys.stdout
        err = sys.stderr

    run_fargo(["start", par_file], np=Nprocs, nt=Nthreads, stdout=out, stderr=err)

    if silent:
        out.close()
        err.close()
        


sigma_slope = 1
#gr by single ring
def gr_single(r, r_0):    
    
    #smoothing parameters
    h = 0.05
    lamb = 0.4571*h + 0.6737*np.sqrt(h)
    chi = -0.7543*h*h + 0.6472*h
    
    eps = lamb**2 * (r_0-r)**2 + chi**2 * r_0*r

    a = (r_0 + r)**2 + eps**2
    m = float(4*r_0*r / a)
    
    da = 2*(r_0 + r) + 2*lamb**2 * (r_0 - r) + chi**2 * r
    dm = 4*r*(a-r_0*da)/a**2

    Km = scipy.special.ellipk(m)
    Em = scipy.special.ellipe(m)

    sigma = 1/r**sigma_slope
    
    g_r =  sigma*r * ( 1/np.sqrt(a) * ( (m-1)*Km + Em )/(2*m*(1-m))*dm - da/(2*a**(3.0/2)) * Km)

    return g_r

        
def get_a_r_theo(r):
    
    #parameters
    disk_mass = 1
    
    #calculate disk density
    mass_norm,_ = scipy.integrate.quad(lambda x: 1/x**(sigma_slope-1), r[0], r[-1])
    mass_norm *= 2*np.pi
    
    gr = np.empty_like(r)

    for idx, r_0 in enumerate(r):
        gr[idx], err = scipy.integrate.quad(gr_single, r[0],
                                        r[-1], args=(r_0,))
        
    #add units
    gr_qty = gr * 4*constants.G * disk_mass/mass_norm * constants.M_sun/constants.au**2
    cms2 = u.cm / u.s / u.s
    gr_qty = gr_qty.to(cms2)
    return gr_qty.value

def test(out1, dt, interactive=False, log=True):


    file_name1 = out1 + f"snapshots/{dt}/a_sg_rad1D.dat"
    data1 = np.fromfile(file_name1)
    a_r_code = data1[1::4].flatten() * 1.6862658982125107 # code to cgs acceleration

    fig, axs = plt.subplots(2,1, gridspec_kw={'height_ratios': [3, 1]})
    fig.subplots_adjust(hspace=0.0)
    ax = axs[0]
    ax1diff = axs[1]

    r = data1[::4]
    r = r.flatten()

    a_r_theo = get_a_r_theo(r)

    ax.axis('auto')
    ax.set_title('Radial Acceleration', color='black', y = 1.06)
    ax.plot(r, a_r_code, '.r', label='Code', lw=2.5)
    ax.plot(r, a_r_theo, '--k', label='Theory', lw=2.5)

    adiff = np.abs(a_r_code / a_r_theo - 1)
    ax1diff.plot(r, adiff)
    ax1diff.set_ylabel('Relative difference')
    
    ## TODO: remove
    ax1diff.set_xlim(5,10)
    ax1diff.set_ylim(np.min(adiff[np.argmin(np.abs(r-5)):]), np.max(adiff[np.argmin(np.abs(r-5)):]))

    ax.legend(loc='upper right')

    if log:
        ax.set_yscale("log", nonpositive='clip')
        ax.set_xscale("log")

    if log:
        ax1diff.set_yscale("log", nonpositive='clip')
        ax1diff.set_xscale("log")

    fig.savefig("plot.jpg", dpi=150)
    if interactive:
        plt.show()


    # check radial acceleration deviation
    threshold = 0.01
    rmin = 1
    rmax = 12.5
    radial_range = np.logical_and(r > rmin, r < rmax)
    max_diff = np.max(adiff[radial_range])
    pass_test = max_diff < threshold
    with open("test.log", "w") as f:
        print(f"Test name: Self gravity", file=f)
        print(f"Max radial acceleration deviation: {max_diff}", file=f)
        print(f"Threshold: {threshold}", file=f)
        print(f"Radial range: {rmin} - {rmax}", file=f)
        print(f"Pass test: {pass_test}", file=f)

    if pass_test:
        print("SUCCESS: Self gravity test")
    else:
        print("FAIL: Self gravity test")

def main():

    parser = ArgumentParser()
    parser.add_argument('-nt', type=int, default=2, help='Number of threads per process.')
    parser.add_argument('-np', type=int, default=1, help='Number of MPI processes.')
    parser.add_argument('-i', '--interactive', action='store_true', help='Interactive mode. Show the plot in a window.')
    parser.add_argument('-p', '--plot', action='store_true', help='Only plot the results.')
    parser.add_argument('-l', '--linear', action='store_true', help='Linear y scaling in plot.')
    parser.add_argument('-s', '--silent', action='store_true', help="Don't print anything.")

    opts = parser.parse_args()

    Nthreads = opts.nt
    Nprocs = opts.np

    if not opts.plot:
        compile_fargo('../../', Nprocs*Nthreads, silent=opts.silent)
        run('../../', 'setup.yml', Nthreads, Nprocs, silent=opts.silent)   
    test('../../output/tests/self_gravity/out/', 0, interactive=opts.interactive, log=False)

if __name__=='__main__':
    main()