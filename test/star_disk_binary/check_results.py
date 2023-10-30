import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import yaml

def test(_, Nsnap=10):

    path =  "../../output/tests/star_disk_binary/"
    path1 = "euler"
    path2 = "leapfrog_dkd"

    quantity1 = 'Sigma'

    string = ''

    filename1 = path + path1 + f"/snapshots/{Nsnap}/" + quantity1 + string + '.dat'
    filename2 = path + path2 + f"/snapshots/{Nsnap}/" + quantity1 + string + '.dat'

    data1 = np.fromfile(filename1)
    data2 = np.fromfile(filename2)

    Nr = 84
    Nphi = 1182
    Nphi = 400

    Nr, Nphi = np.genfromtxt(path + path1 + '/dimensions.dat', usecols=(4,5), dtype=int)


    data1 = np.log(data1.reshape(Nr, Nphi))
    data2 = np.log(data2.reshape(Nr, Nphi))

    r = np.loadtxt(path + path1 + '/used_rad.dat')
    #r = 0.5*(r[1:] + r[:-1])
    #print(np.argmin(np.abs(r-1)))
    phi = np.linspace(0., 2*np.pi, Nphi+1) - np.pi/Nphi

    PHI, R = np.meshgrid(phi, r)

    X = np.cos(PHI)*R
    Y = np.sin(PHI)*R

    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    ax1.set_title(path1)
    ax2.set_title(path2)

    dens = ax1.pcolormesh(X,Y,data1, cmap=plt.get_cmap('plasma'))
    dens2 = ax2.pcolormesh(X,Y,data2, cmap=plt.get_cmap('plasma'))

    ax1.set_xlim(0.7, 1.3)
    ax2.set_xlim(0.7, 1.3)

    ax1.set_ylim(-0.4, 0.4)
    ax2.set_ylim(-0.4, 0.4)


    fig.colorbar(dens, ax=ax1)
    fig.colorbar(dens2, ax=ax2)

    fig.savefig("plot.jpg", dpi=150)


    with open("testconfig.yml", "r") as f:
        testconfig = yaml.safe_load(f)
    testname = testconfig['testname']
    print(f"NO PASSFAIL: {testname}")

if __name__=="__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('Nsnap', type=int, help='Snapshot number.')
    opts = parser.parse_args()

    main("dummy", Nsnap=opts.Nsnap)