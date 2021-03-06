from inversion import *
from utils import *
import time
import globalVar
import matplotlib.pyplot as plt

Inversion_method = Inversion_N_BFGS_modified

def main():
    n = 1
    PATH = 'inverte_T_field/' + 'case' + '{:0>3}'.format(n)
    while os.path.exists(PATH):
        n += 1
        PATH = 'inverte_T_field/' + 'case' + '{:0>3}'.format(n)
    os.mkdir(PATH)


    start_time = time.perf_counter()

    # Td = CN_D(Ts=globalVar.Ts, Tb=globalVar.Tb, kappa=globalVar.kappa, u=globalVar.u, Tic=globalVar.Tic_real)[-1, :]
    # Tic = Inversion_method(Ts=globalVar.Ts, Tb=globalVar.Tb, kappa=globalVar.kappa, u=globalVar.u,
    #                        Td = Td, Tic0=globalVar.Tic_guess, epsilon=globalVar.epsilon, MAX=globalVar.MAX, PATH=PATH)

    # Td = CN_N(Ts=globalVar.Ts, p=globalVar.p, kappa=globalVar.kappa, u=globalVar.u, Tic=globalVar.Tic_real)[-1, :]
    # Tic = Inversion_method(Ts=globalVar.Ts, p=globalVar.p, kappa=globalVar.kappa, u=globalVar.u,
    #                        Td = Td, Tic0=globalVar.Tic_guess, epsilon=globalVar.epsilon, MAX=globalVar.MAX, PATH=PATH)

    # Td = CN_N(Ts=globalVar.Ts, p=globalVar.p, kappa=globalVar.kappa, u=globalVar.u, Tic=globalVar.Tic_real,
    #           sh=globalVar.sh_uniform * 0.3)[-1, :]
    # Td = np.loadtxt('Td.txt')
    # Tic = Inversion_method(Ts=globalVar.Ts, p=globalVar.p, kappa=globalVar.kappa, u=globalVar.u,
    #                        Td = Td, Tic0=globalVar.Tic_guess, epsilon=globalVar.epsilon, MAX=globalVar.MAX, PATH=PATH,
    #                        sh=globalVar.sh_uniform * 0.3)

    Tp = CN_N(Ts=globalVar.Ts, p=globalVar.pm, kappa=globalVar.kappa, u=globalVar.u, Tic=globalVar.Tic_real,
              sh=globalVar.sh_continental(sh0=globalVar.sh0, hr=globalVar.hr, u=globalVar.u))[-1, :]
    Tic = Inversion_method(Ts=globalVar.Ts, p=globalVar.pm, kappa=globalVar.kappa, u=globalVar.u,
                           Tp = Tp, Tic0=globalVar.Tic_guess, epsilon=globalVar.epsilon, MAX=globalVar.MAX, PATH=PATH,
                           sh=globalVar.sh_continental(sh0=globalVar.sh0, hr=globalVar.hr, u=globalVar.u))

    np.savetxt(PATH + '/T.txt', Tic, fmt='%10.5f')

    # plot
    plt.plot(Tp, np.linspace(0, globalVar.zTotal, globalVar.Nz + 1), 'r-', label='Tp')
    plt.plot(globalVar.Tic_real,np.linspace(0, globalVar.zTotal, globalVar.Nz + 1),  'g--', label='Tic_real')
    plt.plot(Tic,np.linspace(0, globalVar.zTotal, globalVar.Nz + 1),  'b-', label='Tic')
    plt.plot(globalVar.Tic_guess, np.linspace(0, globalVar.zTotal, globalVar.Nz + 1), 'b--', label='Tic_guess')
    plt.gca().invert_yaxis()
    plt.xlabel('T')
    plt.ylabel('z')
    plt.text(0.05, 0.05, 'deltaz = {:<5.2e}\n'
                        'deltat = {:<5.2e}\n'
                        'tTotal = {:<7.2f}\n'
                        'epsilon = {:<5.2e}\n'
                        'kappa = {:<7}\n'
                        'qm = {:<7}\n'
                        'rho*H0 = {:<7}\n'
                        'hr = {:<7}\n'
                        'u = {:<7}\n'
                        'Pe = {:<5.3}'.
             format(globalVar.deltaz,
                    globalVar.deltat,
                    globalVar.tTotal,
                    globalVar.epsilon,
                    globalVar.kappa,
                    globalVar.qm,
                    globalVar.rho_H0,
                    globalVar.hr,
                    globalVar.u_mag,
                    globalVar.Pe),
             transform=plt.gca().transAxes)
    plt.legend(loc='upper right')
    plt.title('{} Result'.format(Inversion_method.__name__))

    plt.savefig(PATH + '/case' + '{:0>3}'.format(n) + 'inversion.png')

    end_time = time.perf_counter()
    total_time = end_time - start_time
    print(total_time)
    with open(PATH + '/log.txt', 'a') as file:
        file.write('Total time used: {:<7.2f}\n'.format(total_time))

    plot_log(PATH + '/log.txt', PATH)

    plot_velocity(globalVar.u, mode='time', PATH=PATH, tTotal=globalVar.tTotal)
    pass

if __name__ == '__main__':
    main()