import numpy as np
import matplotlib.pyplot as plt

def test(outputdir, Nsnapshot, interactive=False, log=True):


    file_name1 = outputdir + f"snapshots/{Nsnapshot}/Temperature1D.dat"
    data1 = np.fromfile(file_name1)
    r1 = data1[::4]
    T0 = 106700.1843026118
    Tnum = data1[1::4].flatten() * T0


    fig, axs = plt.subplots(2,2, gridspec_kw={'height_ratios': [1, 1]})
    ax = axs[0,0]
    ax1diff = axs[1,0]
    ax2 = axs[0,1]
    ax2diff = axs[1,1]

    rmax = np.max(r1)
    rmin = np.min(r1)
    vmin = np.min(Tnum)/1.01
    vmax = np.max(Tnum)*1.01


    r = data1[::4]
    r = r.flatten()

    mu = 2.35
    m_H = 1.66054e-24
    k_B = 1.38065e-16
    gamma = 1.4
    l0 = 14959787070000
    rcgs = r * l0
    m0 = 1.98847e+33
    G = 6.6743e-08 # dyne cm^2/g^2

    eta = 2/7
    eps = 0.5
    Rs = 4.6505e-05 * l0

    Ts = 100000

    # note missing factor gamma from eq.(16) fron D'Angelo & Marzari 2012
    # difference comes from them using the adiabatic sound speed for their scaleheight and our code
    # uses the locally isothermal soundspeed for h
    htheo = (eta * (1 - eps) * (k_B * Ts / (mu * m_H))**4 * (Rs / (G * m0))**4 * (rcgs / Rs)**2)**(1/7)
    WG = 0.4 * (Rs / rcgs) + htheo * eta

    Ttheo = Ts * np.sqrt(Rs /rcgs) * ((1-eps)*WG)**(1/4)

    ax.axis('auto')
    # ax.set_title('Temperature', color='black', y = 1.06)

    ax.plot(r, Tnum, '.', color='C0', label='T code', lw=2.5)
    ax.plot(r, Ttheo, '-b', label='T theory', lw=1)
    ax.set_ylabel('log10 T [K]')

    Tdiff = np.abs(Tnum - Ttheo) / Ttheo
    ax1diff.plot(r, Tdiff)
    ax1diff.set_ylabel('Relative difference')

    vmin = min(vmin, np.min(Ttheo))
    vmax = max(vmax, np.max(Ttheo))


    file_name = outputdir + f"snapshots/{Nsnapshot}/aspectratio1D.dat"
    data_dens = np.fromfile(file_name)
    hnum = data_dens[1::4]

    ax2.axis('auto')
    ax2.plot(r, hnum, '.', color='C1', label='$h$', lw=2.5)
    ax2.plot(r, htheo, '-r', label='$h$ theory', lw=1)
    ax2.set_ylabel('h')
    ax.set_xlabel('log10 R [au]')

    vmin2 = np.min(hnum)
    vmax2 = np.max(hnum)

    hdiff = np.abs(hnum - htheo) / htheo
    ax2diff.plot(r, hdiff)
    ax2diff.set_ylabel('Relative difference')




    ax2.legend(loc='best')
    ax.legend(loc='best')

    if log:
        ax.set_yscale("log", nonpositive='clip')
        ax.set_xscale("log")

    if log:
        ax1diff.set_yscale("log", nonpositive='clip')
        ax1diff.set_xscale("log")

    if log:
        ax2.set_yscale("log", nonpositive='clip')
        ax2.set_xscale("log")

    if log:
        ax2diff.set_yscale("log", nonpositive='clip')
        ax2diff.set_xscale("log")


    fig.savefig("plot.jpg", dpi=150)
    if interactive:
        plt.show()

    Tdiff = np.abs(Tnum - Ttheo) / Ttheo

    # check temperature deviation
    testname = "irradiation"
    threshold = 0.03
    rmin = 2
    rmax = 15
    radial_range = np.logical_and(r > rmin, r < rmax)
    max_diff = np.max(Tdiff[radial_range])
    pass_test = max_diff < threshold
    with open("test.log", "w") as f:
        from datetime import datetime
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"{current_time}", file=f)
        print(f"Test name: {testname}", file=f)
        print(f"Max temperature deviation: {max_diff}", file=f)
        print(f"Threshold: {threshold}", file=f)
        print(f"Radial range: {rmin} - {rmax}", file=f)
        print(f"Pass test: {pass_test}", file=f)

    if pass_test:
        print(f"SUCCESS: {testname}")
    else:
        print(f"FAIL: {testname}")

