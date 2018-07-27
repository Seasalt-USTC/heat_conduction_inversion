from inversion import *
import scipy.integrate
import time

def get_thermal_history(T, u):
    """
    Get thermal history of the samples at surface today.
    T is the temperature profile at all time.
    u is the velocity (1D, Function of time.).
    """
    Td = np.zeros(globalVar.Nt + 1, dtype=np.float64)
    depth = scipy.integrate.cumtrapz(np.flip(-u, 0), dx=globalVar.deltat, initial=0)
    depth = np.flip(depth, 0)
    for n in range(globalVar.Nt, -1, -1):
        i = int(depth[n] / globalVar.deltaz)
        # Linear interpolation.
        Td[n] = (depth[n] - i * globalVar.deltaz) / globalVar.deltaz * (T[n, i + 1] - T[n, i]) + T[n, i]
        # Td[n] = T[n, i]
    return Td

def get_velocity_transverse(T, Td):
    """
    Get velocity history from the thermal history of samples at surface today.
    T is the temperature profile at all time.
    Td is the thermal history

    Here use a transversing search method:
        T(t0) for any t0 is monotonic increasing with z.
        Search all z to find a z0 which Td(t0) == T(t0, z0)
    """

    z = np.zeros(globalVar.Nt + 1, dtype=np.float64)
    for n in range(globalVar.Nt + 1):
        for i in range(globalVar.Nz + 1):
            if T[n, i-1] < Td[n] < T[n, i]:
                # Linear interpolation
                z[n] = (i - (T[n, i] - Td[n]) / (T[n, i] - T[n, i-1])) * globalVar.deltaz
    u = np.gradient(z)
    return u

def main():
    n = 1
    PATH = 'invent_data/' + 'case' + '{:0>3}'.format(n)
    while os.path.exists(PATH):
        n += 1
        PATH = 'invent_data/' + 'case' + '{:0>3}'.format(n)
    os.mkdir(PATH)

    u = globalVar.u_gaussian_time
    Tic_real = globalVar.Tic_continent(p=globalVar.pm, sh0=globalVar.sh0, hr=globalVar.hr)

    # Inversion part start from here (invent data from Tp)
    ###################################################################################
    # Tp = CN_N(Ts=globalVar.Ts, p=globalVar.pm, u=u, Tic=Tic_real, kappa=globalVar.kappa,
    #          sh=globalVar.sh_continental(globalVar.sh0, globalVar.hr, u))[-1, :]
    #
    # Inversion_method = Inversion_N_BFGS_modified
    # start_time = time.perf_counter()
    # Tic = Inversion_method(Ts=globalVar.Ts, p=globalVar.pm, kappa=globalVar.kappa,
    #                                 u=u, Td=Tp, Tic0=globalVar.Tic_guess, epsilon=globalVar.epsilon, MAX=globalVar.MAX,
    #                                 PATH=PATH, sh=globalVar.sh_continental(globalVar.sh0, globalVar.hr, u))
    # # plot inversion
    # plt.plot(Tp, np.linspace(0, globalVar.zTotal, globalVar.Nz + 1), 'r-', label='Tp')
    # plt.plot(globalVar.Tic_real,np.linspace(0, globalVar.zTotal, globalVar.Nz + 1),  'g--', label='Tic_real')
    # plt.plot(Tic,np.linspace(0, globalVar.zTotal, globalVar.Nz + 1),  'b-', label='Tic')
    # plt.plot(globalVar.Tic_guess, np.linspace(0, globalVar.zTotal, globalVar.Nz + 1), 'b--', label='Tic_guess')
    # plt.gca().invert_yaxis()
    # plt.xlabel('T')
    # plt.ylabel('z')
    # plt.text(0.05, 0.05, 'deltaz = {:<5.2e}\n'
    #                     'deltat = {:<5.2e}\n'
    #                     'tTotal = {:<7.2f}\n'
    #                     'epsilon = {:<5.2e}\n'
    #                     'kappa = {:<7}\n'
    #                     'qm = {:<7}\n'
    #                     'rho*H0 = {:<7}\n'
    #                     'hr = {:<7}\n'
    #                     'u = {:<7}\n'
    #                     'Pe = {:<5.3}'.
    #          format(globalVar.deltaz,
    #                 globalVar.deltat,
    #                 globalVar.tTotal,
    #                 globalVar.epsilon,
    #                 globalVar.kappa,
    #                 globalVar.qm,
    #                 globalVar.rho_H0,
    #                 globalVar.hr,
    #                 globalVar.u_mag,
    #                 globalVar.Pe),
    #          transform=plt.gca().transAxes)
    # plt.legend(loc='upper right')
    # plt.title('{} Result'.format(Inversion_method.__name__))
    #
    # plt.savefig(PATH + '/case' + '{:0>3}'.format(n) + 'inversion.png')
    #
    # end_time = time.perf_counter()
    # total_time = end_time - start_time
    # print(total_time)
    # with open(PATH + '/log.txt', 'a') as file:
    #     file.write('Total time used: {:<7.2f}\n'.format(total_time))
    #
    # plot_log(PATH + '/log.txt', PATH)
    ###################################################################################
    # Inversion part end here


    # Invent data from Tic_real
    Tic = Tic_real


    '''Invent data.'''
    T = CN_N(Ts=globalVar.Ts, p=globalVar.pm, u=u, Tic=Tic, kappa=globalVar.kappa,
             sh=globalVar.sh_continental(globalVar.sh0, globalVar.hr, u))
    Td = get_thermal_history(T, u[:, 0])
    depth = scipy.integrate.cumtrapz(u[:, 0], dx=globalVar.deltat, initial=0)

    plot_velocity(u, 'time', PATH, tTotal=globalVar.tTotal)

    plt.plot(T[0, :], np.linspace(0, globalVar.zTotal, globalVar.Nz + 1),
             'r-', label='    0 Myr')
    plt.plot(T[int(globalVar.Nt / 2), :], np.linspace(0, globalVar.zTotal, globalVar.Nz + 1),
             'g-', label='{:>5}'.format(int(globalVar.tTotal / 2)) + ' Myr')
    plt.plot(T[-1, :], np.linspace(0, globalVar.zTotal, globalVar.Nz + 1),
             'b-', label='{:>5}'.format(globalVar.tTotal) + ' Myr')
    plt.gca().invert_yaxis()
    plt.legend(loc='lower right')
    plt.title('T')
    plt.savefig(PATH + '/T.png')
    plt.cla()

    plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nt + 1), depth)
    plt.savefig(PATH + '/depth.png')
    plt.cla()

    plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nt + 1), Td)
    plt.savefig(PATH + '/Td.png')
    plt.cla()

    np.savetxt(PATH + '/T.txt', T, fmt='%10.5f')
    np.savetxt(PATH + '/Td.txt', Td, fmt='%10.5f')
    pass

if __name__ == '__main__':
    main()