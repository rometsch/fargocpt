import numpy as np
import scipy
from scipy import optimize
import matplotlib.pyplot as plt
import sys
import os


def load_config(file_path):

    par = {}

    if os.path.isfile(file_path):
        with open(file_path) as f:
            units = f.readlines()

            for line in units:
                l = line.split(' ')
                # l = [x.replace(' ', '') for x in l]
                l = [x.split('\t') for x in l]
                l = [item for sublist in l for item in sublist]
                l = [x.split('\n') for x in l]
                l = [item for sublist in l for item in sublist]
                l = list(filter(None, l))
                # print(l)

                if len(l) > 1:
                    if '#' not in l[0]:
                        if l[1].replace('.','',1).isdigit():
                            if '.' in l[1]:
                                par[l[0]] = float(l[1])
                            else:
                                par[l[0]] = int(l[1])
                        else:
                            par[l[0]] = l[1]

    return par


def get_units(data_folder):
    units_file_name = data_folder + "units.dat"

    unit_dict = {}
    if os.path.isfile(units_file_name):
        with open(units_file_name) as f:
            units = f.readlines()

        for line in units:
            l = line.replace('\t', '')
            l = l.replace('\n', '')
            l = l.split(' ')
            l = list(filter(None, l)) # fastest
            if '=' in l:
                ind = l.index('=')
                l2 = np.array([l[ind-1], l[ind+1]])
            else:
                l2 = ''

            if len(l2) > 1:
                if l2[0] == 'm0':
                    unit_dict['m0'] = np.float(l2[1])
                if l2[0] == 'T0':
                    unit_dict['T0'] = np.float(l2[1])
                if l2[0] == 'Sigma0':
                    unit_dict['Sigma0'] = np.float(l2[1])
                if l2[0] == 'l0':
                    unit_dict['l0'] = np.float(l2[1])
                if l2[0] == 'rho0':
                    unit_dict['rho0'] = np.float(l2[1])
                if l2[0] == 'E0':
                    unit_dict['E0'] = np.float(l2[1])
                if l2[0] == 'v0':
                    unit_dict['v0'] = np.float(l2[1])
                if l2[0] == 'p0':
                    unit_dict['p0'] = np.float(l2[1])
                if l2[0] == 't0':
                    unit_dict['t0'] = np.float(l[ind+1])
    else:
        print("Error in get_units 2D, could not find : " + data_folder + " !\n")
        return 0

    return unit_dict

def cutoff(point, width, x):
	return 1.0/(1.0+np.exp((x-point)/width))


def Ppoly(K, gamma, density):
    return K * (density**gamma)

def Ppoly_optimizer(Kgamma, density, pressure):
    K, gamma = Kgamma
    return np.sum((pressure - Ppoly(K, gamma, density))**2)


def SigmaPoly(Sigma_inner, K, gamma, r, r_inner):

    sigma = ((gamma-1)/(K*gamma) * (1/r - 1/r_inner) + Sigma_inner**(gamma-1))**(1 / (gamma -1))
    return sigma


if __name__ == "__main__":
    args = sys.argv

    if len(args) < 2:
        print("Error: No input file")
        print("Correct usage:\npython get_polytropic_constants.py <path to par file>\n")
        sys.exit()


    file_path = args[1]
    if os.path.isfile(file_path):
        par = load_config(file_path)
    else:
        print("Error: could not find file: ", file_path)
        print("Correct usage:\npython get_polytropic_constants.py <path to par file>")
        sys.exit()

    GlobalNRadial = par['Nrad']
    Radii = np.zeros(GlobalNRadial+1)
    Rmed = np.zeros(GlobalNRadial)
    RMIN = par['Rmin']
    RMAX = par['Rmax']

    for nRadial in range(GlobalNRadial+1):
        Radii[nRadial] = RMIN*np.exp((nRadial-1.0)/(GlobalNRadial-2.0)*np.log(RMAX/RMIN))

    Rsup = Radii[1:]
    Rinf = Radii[:-1]

    for nRadial in range(GlobalNRadial):
        Rmed[nRadial] = 2.0/3.0*(Radii[nRadial+1]*Radii[nRadial+1]*Radii[nRadial+1]-Radii[nRadial]*Radii[nRadial]*Radii[nRadial])
        Rmed[nRadial] = Rmed[nRadial] / (Radii[nRadial+1]*Radii[nRadial+1]-Radii[nRadial]*Radii[nRadial])


    sigma0 = par['Sigma0']
    SIGMASLOPE = par['SigmaSlope']

    density = np.empty_like(Rmed)
    for n_radial in range(GlobalNRadial):
        density[n_radial] = sigma0*pow(Rmed[n_radial],-SIGMASLOPE)


    if par['ProfileDamping'] == 'Yes':
        print('ProfileDamping')
        for n_radial in range(len(density)):
            density[n_radial] *= cutoff(par['ProfileDampingPoint'], par['ProfileDampingWidth'], Rmed[n_radial])

    if par['SetSigma0'] == 'Yes':
        ring_area = np.pi * (Rsup**2 - Rinf**2)
        total_mass = 0
        for dens, surface in zip(density[1:-1], ring_area[1:-1]):
            total_mass += dens*surface

        sigma0 *= par['DiscMass']/total_mass

        for n_radial in range(len(density)):
            density[n_radial] *= par['DiscMass']/total_mass


    soundspeed = np.empty_like(density)
    for n_radial in range(len(density)):
        soundspeed[n_radial] = par['AspectRatio']*np.sqrt(1/Rmed[n_radial])*pow(Rmed[n_radial], par['FlaringIndex'])

    pressure = np.empty_like(density)
    for n_radial in range(len(density)):
        pressure[n_radial] = density[n_radial] * soundspeed[n_radial]**2

    temperature = np.empty_like(density)
    for n_radial in range(len(density)):
        temperature[n_radial] = par['mu']*pressure[n_radial]/density[n_radial]


    xatol = 1e-14
    f = lambda Kgamma: Ppoly_optimizer(Kgamma, density, pressure)
    K, gamma = scipy.optimize.minimize(f, x0=[12, 2], method='Powell', tol=xatol).x

    print('Input file: ', file_path)
    print('PolytropicConstant   ', K, '\nAdiabaticIndex       ', gamma)



    # Plotting for Testing
    DO_PLOT = False
    if DO_PLOT:
        print(par)
        out_folder = 'out'
        units = get_units('../' + out_folder + '/')

        time = [100]

        for dt in time:
            fargo_v_azi = np.fromfile('../' + out_folder + '/gasvtheta1D' + str(dt) + '.dat')
            fargo_v_azi_radius = fargo_v_azi[::4]
            fargo_v_azi = fargo_v_azi[1::4]

            print(dt, 'v_azi = ', np.sum(np.abs(fargo_v_azi)*fargo_v_azi_radius))

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title('density')
        ax.plot(Rmed, density, '-b', label='Iso')

        for dt in time:
            fargo_dens = np.fromfile('../' + out_folder + '/gasdens1D' + str(dt) + '.dat')
            fargo_dens_radius = fargo_dens[::4]
            fargo_dens = fargo_dens[1::4]/units['Sigma0']
            ax.plot(fargo_dens_radius, fargo_dens, '--m')#, label='Fargo')
            total_mass = 0
            for dens, surface in zip(fargo_dens[1:-1], ring_area[1:-1]):
                total_mass += dens*surface
            print('Mass = ', total_mass)

            Sigma0 = fargo_dens[1]
            r_inner = fargo_dens_radius[1]
            theo_dens = SigmaPoly(Sigma0, par['PolytropicConstant'], par['AdiabaticIndex'], fargo_dens_radius, r_inner)
            ax.plot(fargo_dens_radius, theo_dens, '--g')#, label='Fargo')
        ax.legend(loc='upper right')


        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title('Pressure')
        ax.plot(Rmed, pressure, '-b', label='Iso')
        ax.plot(Rmed, Ppoly(K, gamma, density), '--r', label='Poly')

        for dt in time:
            fargo_press = np.fromfile('../' + out_folder + '/gaspressure1D' + str(dt) + '.dat')
            fargo_press_radius = fargo_press[::4]
            fargo_press = fargo_press[1::4]
            ax.plot(fargo_press_radius, fargo_press/units['p0'], '--m')#, label='Fargo')

        ax.legend(loc='upper right')


        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title('Temperature')
        ax.plot(Rmed, temperature, '-b', label='Iso')
        ax.plot(Rmed, (par['mu'] * Ppoly(K, gamma, density)/density), '--r', label='Poly')

        for dt in time:
            fargo_temp = np.fromfile('../' + out_folder + '/gasTemperature1D' + str(dt) +  '.dat')
            fargo_temp_radius = fargo_temp[::4]
            fargo_temp = fargo_temp[1::4]
            ax.plot(fargo_temp_radius, fargo_temp/units['T0'], '--m')#, label='Fargo')

        ax.legend(loc='upper right')
        plt.show()



        # RH2D = np.loadtxt("/home/jordan/Downloads/1dx.00000", skiprows=1, unpack=True)
        # RH2D_radius = RH2D[2,:]

        # RH2D_old = np.loadtxt("/home/jordan/Downloads/1dx(1).00000", skiprows=1, unpack=True)
        # RH2D_radius_old = RH2D_old[2,:]



        # fig=plt.figure()
        # ax = fig.add_subplot(111)
        # ax.plot(Radii/RH2D_radius, '-b', label='Fargo/RH2D NEW')
        # ax.plot(Radii/RH2D_radius_old, '--r', label='Fargo/RH2D OLD')
        # # ax.plot(Radii, '-b', label='Fargo')
        # # ax.plot(RH2D_radius, '--r', label='RH2D')
        # ax.legend(loc='upper left')
        # ax.set_xlabel('Fargo radius [au]')
        # ax.set_ylabel('Radius/Radius')
        # # plt.savefig('compare_radii.png', dpi=120)
        # plt.show()

