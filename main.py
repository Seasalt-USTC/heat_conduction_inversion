from inversion import *
from utils import *
import time
import globalVar
import matplotlib.pyplot as plt

Inversion_method = Inversion_N_BFGS_root

def main():
    n = 1
    PATH = 'case/' + 'case' + '{:0>3}'.format(n)
    while os.path.exists(PATH):
        n += 1
        PATH = 'case/' + 'case' + '{:0>3}'.format(n)
    os.mkdir(PATH)

    start_time = time.clock()

    # Td = CN_D(Ts=globalVar.Ts, Tb=globalVar.Tb, kappa=globalVar.kappa, u=u(PATH, True), Tic=globalVar.Tic_real())[-1, :]
    # Tic = Inversion_method(Ts=globalVar.Ts, Tb=globalVar.Tb, kappa=globalVar.kappa, u=u(),
    #                   Td = Td, Tic0=globalVar.Tic_guess(), epsilon=globalVar.epsilon, MAX=globalVar.MAX, PATH=PATH)

    Td = CN_N(Ts=0, p=1, kappa=0.05, u=globalVar.u0, Tic=globalVar.Tic_real())[-1, :]
    Tic = Inversion_method(Ts=0, p=1, kappa=0.05, u=globalVar.u0,
                      Td = Td, Tic0=globalVar.Tic_guess(), epsilon=globalVar.epsilon, MAX=globalVar.MAX, PATH=PATH)

    # plot
    plt.plot(np.linspace(0, 1, globalVar.Nz + 1), Td, 'r-', label='Td')
    plt.plot(np.linspace(0, 1, globalVar.Nz + 1), globalVar.Tic_real(), 'g--', label='globalVar.Tic_real')
    plt.plot(np.linspace(0, 1, globalVar.Nz + 1), Tic, 'b-', label='Tic')
    plt.plot(np.linspace(0, 1, globalVar.Nz + 1), globalVar.Tic_guess(), 'b--', label='globalVar.Tic_guess')
    plt.xlabel('z')
    plt.ylabel('T')
    plt.text(0, 0.7, 't = {:<7}\nkappa = {:<7.5}\n'
                     'u = {:<}\nglobalVar.Tic_real = {:<}\nglobalVar.Tic_guess = {:<}\n'
                     'Set Ts always positive: {}\n'
                     'Epsilon = {:<10.0e}'.
             format(globalVar.tTotal, globalVar.kappa,
                    globalVar.u.__name__, globalVar.Tic_real.__name__, globalVar.Tic_guess.__name__,
                    globalVar.set_T0_always_positive,
                    globalVar.epsilon))
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