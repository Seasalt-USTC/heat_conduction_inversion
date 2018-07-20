from inversion import *
from utils import *
import time
import globalVar
import matplotlib.pyplot as plt

Inversion_method = Inversion_N_Steepest

def main():
    n = 1
    PATH = 'case/' + 'case' + '{:0>3}'.format(n)
    while os.path.exists(PATH):
        n += 1
        PATH = 'case/' + 'case' + '{:0>3}'.format(n)
    os.mkdir(PATH)


    start_time = time.clock()

    # Td = CN_D(Ts=globalVar.Ts, Tb=globalVar.Tb, kappa=globalVar.kappa, u=globalVar.u, Tic=globalVar.Tic_real)[-1, :]
    # Tic = Inversion_method(Ts=globalVar.Ts, Tb=globalVar.Tb, kappa=globalVar.kappa, u=globalVar.u,
    #                   Td = Td, Tic0=globalVar.Tic_guess, epsilon=globalVar.epsilon, MAX=globalVar.MAX, PATH=PATH)

    # Td = CN_N(Ts=globalVar.Ts, p=globalVar.p, kappa=globalVar.kappa, u=globalVar.u, Tic=globalVar.Tic_real)[-1, :]
    # Tic = Inversion_method(Ts=globalVar.Ts, p=globalVar.p, kappa=globalVar.kappa, u=globalVar.u,
    #                   Td = Td, Tic0=globalVar.Tic_guess, epsilon=globalVar.epsilon, MAX=globalVar.MAX, PATH=PATH)

    Td = CN_N(Ts=globalVar.Ts, p=globalVar.pm, kappa=globalVar.kappa, u=globalVar.u, Tic=globalVar.Tic_real,
              sh=globalVar.sh_continental(sh0=globalVar.sh0, hr=globalVar.hr, u=globalVar.u))[-1, :]
    Tic = Inversion_method(Ts=globalVar.Ts, p=globalVar.pm, kappa=globalVar.kappa, u=globalVar.u,
                           Td = Td, Tic0=globalVar.Tic_guess, epsilon=globalVar.epsilon, MAX=globalVar.MAX, PATH=PATH,
                           sh=globalVar.sh_continental(sh0=globalVar.sh0, hr=globalVar.hr, u=globalVar.u))

    np.savetxt(PATH + '/T.txt', Tic, fmt='%10.5f')

    # plot
    plt.plot(np.linspace(0, globalVar.zTotal, globalVar.Nz + 1), Td, 'r-', label='Td')
    plt.plot(np.linspace(0, globalVar.zTotal, globalVar.Nz + 1), globalVar.Tic_real, 'g--', label='Tic_real')
    plt.plot(np.linspace(0, globalVar.zTotal, globalVar.Nz + 1), Tic, 'b-', label='Tic')
    plt.plot(np.linspace(0, globalVar.zTotal, globalVar.Nz + 1), globalVar.Tic_guess, 'b--', label='Tic_guess')
    plt.xlabel('z')
    plt.ylabel('T')
    plt.text(0.05, 0.6, 'zTotal = {:<7}\n'
                        'qm = {:<7}\nk = {:<7}\nkappa = {:<7}\n'
                        'u = {:<7}\nrho*H0 = {:<7}\nhr = {:<7}\n'
                        'Pe = {:<5.3}'.
             format(globalVar.zTotal,
                    globalVar.qm, globalVar.k, globalVar.kappa,
                    globalVar.u_mag, globalVar.rho_H0, globalVar.hr,
                    globalVar.Pe),
             transform=plt.gca().transAxes)
    plt.legend(loc='upper right')
    plt.title('{} Result'.format(Inversion_method.__name__))

    plt.savefig(PATH + '/case' + '{:0>3}'.format(n) + 'inversion.png')

    end_time = time.clock()
    total_time = end_time - start_time
    print(total_time)
    with open(PATH + '/log.txt', 'a') as file:
        file.write('Total time used: {}\n'.format(total_time))

    plot_log(PATH + '/log.txt', PATH)
    pass

if __name__ == '__main__':
    main()