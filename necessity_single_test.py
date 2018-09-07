from invent_data import *
import globalVar
import matplotlib.pyplot as plt

# Inversion_method = Inversion_N_BFGS_modified

def main():
    n = 5
    u  = np.ones((globalVar.Nt + 1, globalVar.Nz + 1), dtype=np.float32)
    PATH = 'necessity/' + 'case' + '{:0>3}'.format(n)
    n += 1

    sigma = 20
    u0 = 5

    # Build a synthetic velocity data.
    for k in range(globalVar.Nt + 1):
        tk = k * globalVar.deltat
        u[k, :] = - globalVar.gauss(sigma ** 2, 50, tk) * u0

    # Construct the temperature field.
    T_real = CN_N(Ts=globalVar.Ts, p=globalVar.pm, kappa=globalVar.kappa, u=u, Tic=globalVar.Tic_real,
              sh=globalVar.sh_continental(sh0=globalVar.sh0, hr=globalVar.hr, u=u))
    plt.plot(np.linspace(0, globalVar.zTotal, globalVar.Nz + 1), T_real[0, :],
             'r-', label='    0 Ma')
    plt.plot(np.linspace(0, globalVar.zTotal, globalVar.Nz + 1), T_real[int(globalVar.Nt / 2), :],
             'g-', label='{:>5}'.format(int(globalVar.tTotal / 2)) + ' Ma')
    plt.plot(np.linspace(0, globalVar.zTotal, globalVar.Nz + 1), T_real[-1, :],
             'b-', label='{:>5}'.format(globalVar.tTotal) + ' Ma')
    plt.legend(loc='lower right')
    plt.title('T_real')
    plt.savefig(PATH + 'T_real.png')
    plt.cla()

    # Obtain the thermal history.
    Td = get_thermal_history(T_real, u[:, 0])
    plt.plot(Td)
    plt.title('Td')
    plt.savefig(PATH + 'Td.png')
    plt.cla()

    # Construct a geotherm from surface heat flow.
    # Use this geotherm with the assumption that (partial T/partial t) == 0 to recover the velocity.
    Tp = T_real[-1, :]
    ps = (Tp[1] - Tp[0]) / globalVar.deltaz
    sh0 = (ps - globalVar.pm) * globalVar.kappa / globalVar.hr
    T_reconstructed_by_qs =  globalVar.Tic_continent(globalVar.pm, sh0, globalVar.hr)
    T_reconstructed_by_qs = np.tile(T_reconstructed_by_qs, (globalVar.Nt + 1, 1))

    # plt.plot(Tp, color='k', label='Tp')
    # plt.plot(T_reconstructed_by_qs[0, :], color='r', label='T_qs')
    # plt.legend()
    # plt.savefig(PATH + 'T_qs.png')
    # plt.cla()

    # u_r_present = get_velocity_transverse(np.tile(Tp, (globalVar.Nt + 1, 1)), Td)
    # plt.plot(u_r_present, color='g', label='u_r_present')


    # Tic = Inversion_method(Ts=globalVar.Ts, p=globalVar.pm, kappa=globalVar.kappa, u=u,
    #                        Tp=Tp, Tic0=globalVar.Tic_guess, epsilon=globalVar.epsilon, MAX=globalVar.MAX, PATH=PATH,
    #                        sh=globalVar.sh_continental(sh0=globalVar.sh0, hr=globalVar.hr, u=u))
    # T = CN_N(Ts=globalVar.Ts, p=globalVar.pm, kappa=globalVar.kappa, u=u, Tic=Tic,
    #           sh=globalVar.sh_continental(sh0=globalVar.sh0, hr=globalVar.hr, u=u))
    # np.savetxt(PATH + '/T.txt', T)
    #
    # T = np.loadtxt(PATH + '/T.txt')
    # u_r_inversion = get_velocity_transverse(T, Td)
    # plt.plot(u_r_inversion, color='r', label='u_r_inversion')

    u_r_qs = get_velocity_transverse(T_reconstructed_by_qs, Td)
    plt.plot(u_r_qs, color='b', label='u_r_qs')
    plt.plot(u[:, 0], color='k', label='u_real')
    plt.legend()
    plt.title('velocity reconstructed')
    plt.savefig(PATH + 'u.png')
    plt.cla()

    J_value = np.sqrt(scipy.integrate.simps(np.array(u_r_qs - u[:, 0]) ** 2, dx=globalVar.deltat) / globalVar.tTotal)
    print(sigma, u0)
    plt.figure(figsize=(10, 7))
    # plt.imshow(J, extent=[0, 5, 1, 20], origin='lower')
    grid_sigma, grid_u0 = np.mgrid[1:20:5j, 0:5:6j]
    plt.contourf(grid_sigma, grid_u0, J)
    plt.xlabel('sigma')
    plt.xticks(np.append([1], np.linspace(0, 20, 6)))
    plt.xlim(1, 20)
    plt.ylabel('u0')
    plt.yticks(np.linspace(0, 5, 6))
    plt.ylim(0, 5)
    plt.title('J')
    plt.colorbar()
    plt.savefig('necessity/J.png')
    plt.show()

if __name__ == '__main__':
    main()