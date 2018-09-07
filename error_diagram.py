from invent_data import *
import globalVar
import matplotlib.pyplot as plt
# from scipy import integrate

# velocity_model = 'Gaussian'
velocity_model = 'Linear'

def main():
    n = 1
    u  = np.ones((globalVar.Nt + 1, globalVar.Nz + 1), dtype=np.float32)
    J = []
    # sigmas = np.linspace(0, 100, 41)
    # u0s = np.linspace(0, 1, 41)
    sigmas = np.linspace(0, 100, 6)
    u0s = np.linspace(0, 1, 6)
    for sigma in sigmas:
        J_u0 = []
        for u0 in u0s:
            n += 1

            # Build a synthetic velocity data.
            if velocity_model == 'Gaussian':
                for k in range(globalVar.Nt + 1):
                    tk = k * globalVar.deltat
                    u[k, :] = - globalVar.gauss(sigma ** 2, 50, tk) * u0
            elif velocity_model == 'Linear':
                for k in range(globalVar.Nt + 1):
                    tk = k * globalVar.deltat
                    if tk <= globalVar.tTotal - sigma:
                        u[k, :] = 0
                    else:
                        u[k, :] = - (tk - (globalVar.tTotal - sigma)) / sigma * u0

            # Construct the temperature field.
            T_real = CN_N(Ts=globalVar.Ts, p=globalVar.pm, kappa=globalVar.kappa, u=u, Tic=globalVar.Tic_real,
                      sh=globalVar.sh_continental(sh0=globalVar.sh0, hr=globalVar.hr, u=u))

            # Obtain the thermal history.
            Td = get_thermal_history(T_real, u[:, 0])

            # Construct a geotherm from surface heat flow.
            # Use this geotherm with the assumption that (partial T/partial t) == 0 to recover the velocity.
            Tp = T_real[-1, :]
            ps = (Tp[1] - Tp[0]) / globalVar.deltaz
            sh0 = (ps - globalVar.pm) * globalVar.kappa / globalVar.hr
            T_reconstructed_by_qs =  globalVar.Tic_continent(globalVar.pm, sh0, globalVar.hr)
            T_reconstructed_by_qs = np.tile(T_reconstructed_by_qs, (globalVar.Nt + 1, 1))

            u_r_qs = get_velocity_transverse(T_reconstructed_by_qs, Td)

            plt.plot(np.linspace(0, 120, 61), u_r_qs, color='b', label='u_r_qs')
            plt.plot(np.linspace(0, 120, 61), u[:, 0], color='k', label='u_real')
            plt.legend()
            plt.title('velocity reconstructed')
            plt.savefig(f'necessity/cases/u{n}.png')
            plt.cla()

            J_value = np.sqrt(scipy.integrate.simps(np.array(u_r_qs - u[:, 0])**2, dx=globalVar.deltat) / globalVar.tTotal)
            J_u0.append(J_value)
            print(sigma, u0)
        J.append(J_u0)
    J = np.array(J)
    np.savetxt('J.txt', J)
    J = np.loadtxt('J.txt')
    plt.figure(figsize=(10, 7))
    plt.imshow(J, extent=[0, 1, 0, 100], origin='lower', aspect='auto')
    # plt.imshow(J_u0, extent=[0, 1, 0, 100], origin='lower', aspect='auto')
    plt.colorbar()
    grid_u0, grid_sigma = np.meshgrid(u0s, sigmas)
    CS = plt.contour(grid_u0, grid_sigma, J, colors='k', levels=np.linspace(0, 0.2, 21))
    # CS = plt.contour(grid_u0, grid_sigma, J_u0, colors='k', levels=np.linspace(0, 1, 11))
    plt.clabel(CS)
    plt.xlabel('u0')
    plt.xticks(np.linspace(0, 1, 6))
    plt.xlim(0, 1)
    plt.ylabel('sigma')
    plt.yticks(np.linspace(0, 100, 6))
    plt.ylim(0, 100)
    plt.title('J')
    plt.savefig('necessity/J.png')

    J_u0 = J / np.tile(u0s, (41, 1))
    J_u0[:, 0] = 0
    plt.figure(figsize=(10, 7))
    plt.imshow(J_u0, extent=[0, 1, 0, 100], origin='lower', aspect='auto')
    plt.colorbar()
    grid_u0, grid_sigma = np.meshgrid(u0s, sigmas)
    CS = plt.contour(grid_u0, grid_sigma, J_u0, colors='k', levels=np.linspace(0, 1, 51))
    plt.clabel(CS)
    plt.xlabel('u0')
    plt.xticks(np.linspace(0, 1, 6))
    plt.xlim(0, 1)
    plt.ylabel('sigma')
    plt.yticks(np.linspace(0, 100, 6))
    plt.ylim(0, 100)
    plt.title('J/u0')
    plt.savefig('necessity/J_.png')
    # plt.show()

if __name__ == '__main__':
    main()