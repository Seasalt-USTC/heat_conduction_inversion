from invent_data import *
import globalVar
import matplotlib.pyplot as plt
from scipy import integrate

# Inversion_method = Inversion_N_BFGS_modified

def main():
    n = 1
    u  = np.ones((globalVar.Nt + 1, globalVar.Nz + 1), dtype=np.float32)
    J = []
    sigmas = np.linspace(1, 20, 39)
    u0s = np.linspace(0, 1, 21)
    for sigma in sigmas:
        J_u0 = []
        for u0 in u0s:
            # PATH = 'necessity/' + 'case' + '{:0>3}'.format(n)
            n += 1

            # Build a synthetic velocity data.
            for k in range(globalVar.Nt + 1):
                tk = k * globalVar.deltat
                u[k, :] = - globalVar.gauss(sigma ** 2, 70, tk) * u0

            # Construct the temperature field.
            T_real = CN_N(Ts=globalVar.Ts, p=globalVar.pm, kappa=globalVar.kappa, u=u, Tic=globalVar.Tic_real,
                      sh=globalVar.sh_continental(sh0=globalVar.sh0, hr=globalVar.hr, u=u))
            # plt.plot(np.linspace(0, globalVar.zTotal, globalVar.Nz + 1), T_real[0, :],
            #          'r-', label='    0 Ma')
            # plt.plot(np.linspace(0, globalVar.zTotal, globalVar.Nz + 1), T_real[int(globalVar.Nt / 2), :],
            #          'g-', label='{:>5}'.format(int(globalVar.tTotal / 2)) + ' Ma')
            # plt.plot(np.linspace(0, globalVar.zTotal, globalVar.Nz + 1), T_real[-1, :],
            #          'b-', label='{:>5}'.format(globalVar.tTotal) + ' Ma')
            # plt.legend(loc='lower right')
            # plt.title('T_real')
            # plt.savefig(PATH + 'T_real.png')
            # plt.cla()

            # Obtain the thermal history.
            Td = get_thermal_history(T_real, u[:, 0])
            # plt.plot(Td)
            # plt.title('Td')
            # plt.savefig(PATH + 'Td.png')
            # plt.cla()

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
            # plt.plot(u_r_qs, color='b', label='u_r_qs')
            # plt.plot(u[:, 0], color='k', label='u_real')
            # plt.legend()
            # plt.title('velocity reconstructed')
            # plt.savefig(PATH + 'u.png')
            # plt.cla()

            J_value = np.sqrt(scipy.integrate.simps(np.array(u_r_qs - u[:, 0])**2, dx=globalVar.deltat) / globalVar.tTotal)
            J_u0.append(J_value)
            print(sigma, u0)
        J.append(J_u0)
    J = np.array(J)
    np.savetxt('J.txt', J)
    J = np.loadtxt('J.txt')

    # plt.figure(figsize=(10, 7))
    # plt.imshow(J, extent=[0, 1, 1, 20], origin='lower', aspect='auto')
    # # plt.imshow(J_u0, extent=[0, 1, 1, 20], origin='lower', aspect='auto')
    # plt.colorbar()
    # grid_u0, grid_sigma = np.meshgrid(u0s, sigmas)
    # # plt.contourf(grid_sigma, grid_u0, J)
    # CS = plt.contour(grid_u0, grid_sigma, J, colors='k', levels=np.linspace(0, 1, 51))
    # # CS = plt.contour(grid_u0, grid_sigma, J_u0, colors='k', levels=np.linspace(0, 1, 101))
    # plt.clabel(CS)
    # plt.xlabel('u0')
    # plt.xticks(np.linspace(0, 1, 6))
    # plt.xlim(0, 1)
    # plt.ylabel('sigma')
    # plt.yticks(np.append([1], np.linspace(0, 20, 6)))
    # plt.ylim(1, 20)
    # # plt.title('J')
    # # plt.savefig('necessity/J.png')
    # plt.title('J')
    # plt.savefig('necessity/J.png')
    # plt.show()
    #
    # J_u0 = J / np.tile(u0s, (39, 1))
    # J_u0[:, 0] = 0
    # plt.figure(figsize=(10, 7))
    # # plt.imshow(J, extent=[0, 5, 1, 20], origin='lower', aspect='auto')
    # plt.imshow(J_u0, extent=[0, 1, 1, 20], origin='lower', aspect='auto')
    # plt.colorbar()
    # grid_u0, grid_sigma = np.meshgrid(u0s, sigmas)
    # # plt.contourf(grid_sigma, grid_u0, J)
    # # CS = plt.contour(grid_u0, grid_sigma, J, colors='k', levels=np.linspace(0, 5, 11))
    # CS = plt.contour(grid_u0, grid_sigma, J_u0, colors='k', levels=np.linspace(0, 1, 51))
    # plt.clabel(CS)
    # plt.xlabel('u0')
    # plt.xticks(np.linspace(0, 1, 6))
    # plt.xlim(0, 1)
    # plt.ylabel('sigma')
    # plt.yticks(np.append([1], np.linspace(0, 20, 6)))
    # plt.ylim(1, 20)
    # # plt.title('J')
    # # plt.savefig('necessity/J.png')
    # plt.title('J/u0')
    # plt.savefig('necessity/J_u0.png')
    # plt.show()

    FWHMs = sigmas * 2.355
    grid_u0, grid_FWHMs = np.meshgrid(u0s, FWHMs)
    plt.figure(figsize=(10, 7))
    plt.imshow(J, extent=[0, 1, 1, 20*2.355], origin='lower', aspect='auto')
    # plt.imshow(J_u0, extent=[0, 1, 1, 20], origin='lower', aspect='auto')
    plt.colorbar()
    CS = plt.contour(grid_u0, grid_FWHMs, J, colors='k', levels=np.linspace(0, 1, 51))
    plt.clabel(CS)
    plt.xlabel('u0')
    plt.xticks(np.linspace(0, 1, 6))
    plt.xlim(0, 1)
    plt.ylabel('FWHM')
    plt.yticks(np.arange(1, 20*2.355, 4))
    # plt.ylim(1, 20)
    plt.title('J')
    plt.savefig('necessity/J.png')
    plt.show()

    J_u0 = J / np.tile(u0s, (39, 1))
    J_u0[:, 0] = 0
    plt.figure(figsize=(10, 7))
    # plt.imshow(J, extent=[0, 5, 1, 20], origin='lower', aspect='auto')
    plt.imshow(J_u0, extent=[0, 1, 1, 20*2.355], origin='lower', aspect='auto')
    plt.colorbar()
    CS = plt.contour(grid_u0, grid_FWHMs, J_u0, colors='k', levels=np.linspace(0, 1, 51))
    plt.clabel(CS)
    plt.xlabel('u0')
    plt.xticks(np.linspace(0, 1, 6))
    plt.xlim(0, 1)
    plt.ylabel('FWHM')
    plt.yticks(np.arange(1, 20*2.355, 4))
    # plt.ylim(1, 20)
    # plt.title('J')
    # plt.savefig('necessity/J.png')
    plt.title('J/u0')
    plt.savefig('necessity/J_u0.png')
    plt.show()

if __name__ == '__main__':
    main()