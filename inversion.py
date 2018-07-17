import os
import time
from forward import *

# TODO: modify
def Inversion_D_Steepest(T0, T1, kappa, u, Td, Tic0, epsilon, MAX, PATH):
    """
    Inverse heat conduction equation with the iteration of adjoint equation method
    from Td back to 0.
    Td is the data obtained today.
    Tic0 is a initial guess of the iteration.
    """

    logFile = PATH + '/log.txt'
    T0k = Tic0
    for k in range(MAX):
        T1k = CN_D(T0=T0, T1=T1, kappa=kappa, u=u, Tic=T0k)[-1, :]

        if not os.path.exists(PATH + '/T0'):
            os. mkdir(PATH + '/T0')
        if not os.path.exists(PATH + '/T1'):
            os. mkdir(PATH + '/T1')
        if k % 50 == 0:
            plt.ylim(-0.1, 1.1)
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), T0k, 'r-', label='Tick')
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), globalVar.Tic_real(), 'b-', label='Tic_real')
            no = '{:0>4}'.format(str(k))
            plt.xlabel('z')
            plt.ylabel('T0')
            plt.legend(loc='upper right')
            plt.savefig(PATH + '/T0/'+ no + 'T0.png')
            plt.cla()

            plt.ylim(-0.1, 1.1)
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), T1k, 'r-', label='T1n')
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), Td, 'b-', label='Td')
            no = '{:0>4}'.format(str(k))
            plt.xlabel('z')
            plt.ylabel('T1')
            plt.legend(loc='upper right')
            plt.savefig(PATH + '/T1/'+ no + 'T1.png')
            plt.cla()

        Jk = norm_2(T1k - Td)
        with open(logFile, 'a') as file:
            file.write('J{:<4} = {:<9.7}\n'.format(k, Jk))
        if Jk <= epsilon:
            Tic = T0k
            break

        lambda0k = CN_D_B(T0=0, T1=0, kappa=kappa, u=u, Tec=2 * (T1k - Td))[0, :]
        alphak = 1
        T0k = T0k - alphak * lambda0k

        if globalVar.set_T0_always_positive:
            for i in range(globalVar.Nz + 1):
                if T0k[i] < 0:
                    T0k[i] = 0
    else:
        Tic = T0k
    return Tic
def Inversion_D_DFP(T0, T1, kappa, u, Td, Tic0, epsilon, MAX, PATH):
    """
    Inverse heat conduction equation with the iteration of adjoint equation method.
    Minimize J(objective function wich represents the mismach between the prediction and the data)
        with BFGS Quasi-Newton Method.
    From Td back to 0.
        Td is the data obtained today.
        Tic0 is a initial guess of the iteration.
    """

    logFile = PATH + '/log.txt'
    for k in range(MAX):
        if k == 0:
            T0k = Tic0
            Vk = np.eye(globalVar.Nz + 1, dtype=np.float64)
            T1k = CN_D(T0=T0, T1=T1, kappa=kappa, u=u, Tic=T0k)[-1, :]
            lambda0k = CN_D_B(T0=0, T1=0, kappa=kappa, u=u, Tec=2 * (T1k - Td))[0, :]
        else:
            T0k = T0k1
            Vk = Vk1
            T1k = T1k1
            lambda0k = lambda0k1
        if not os.path.exists(PATH + '/T0'):
            os. mkdir(PATH + '/T0')
        if not os.path.exists(PATH + '/T1'):
            os. mkdir(PATH + '/T1')
        if k % 50 == 0:  # plot
            plt.ylim(-0.1, 1.1)
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), T0k, 'r-', label='Tick')
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), globalVar.Tic_real(), 'b-', label='Tic_real')
            no = '{:0>4}'.format(str(k))
            plt.xlabel('z')
            plt.ylabel('T0')
            plt.legend(loc='upper right')
            plt.savefig(PATH + '/T0/'+ no + 'T0.png')
            plt.cla()

            plt.ylim(-0.1, 1.1)
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), T1k, 'r-', label='T1n')
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), Td, 'b-', label='Td')
            no = '{:0>4}'.format(str(k))
            plt.xlabel('z')
            plt.ylabel('T1')
            plt.legend(loc='upper right')
            plt.savefig(PATH + '/T1/'+ no + 'T1.png')
            plt.cla()

        Jk = norm_2(T1k - Td)
        with open(logFile, 'a') as file:
            file.write('J{:<4} = {:<9.7}\n'.format(k, Jk))
        if Jk < epsilon:
            Tic = T0k
            break

        pk = - np.dot(Vk, lambda0k)
        alphak = 1  # TODO: line search to get a step length alphak.
        sk = alphak * pk
        T0k1 = T0k + sk
        if globalVar.set_T0_always_positive:
            for i in range(globalVar.Nz + 1):
                if T0k1[i] < 0:
                    T0k1[i] = 0
        T1k1 = CN_D(T0=T0, T1=T1, kappa=kappa, u=u, Tic=T0k1)[-1, :]
        # Jk1 = norm_2(T1k1 - Td)
        # while Jk1 > Jk:
        #     alphak /= 2
        #     sk = alphak * pk
        #     T0k1 = T0k + sk
        #     if globalVar.set_T0_always_positive:
        #         for i in range(globalVar.Nz + 1):
        #             if T0k1[i] < 0:
        #                 T0k1[i] = 0
        #     T1k1 = CN_D(T0=T0, T1=T1, kappa=kappa, u=u, Tic=T0k1)[-1, :]
        #     Jk1 = norm_2(T1k1 - Td)
        # print(alphak)
        lambda0k1 = CN_D_B(T0=0, T1=0, kappa=kappa, u=u, Tec=2 * (T1k1 - Td))[0, :]
        yk = lambda0k1 - lambda0k
        Vk1 = Vk \
              - np.dot( np.outer( np.dot(Vk, yk), yk), Vk) / np.inner(np.dot(yk, Vk), yk) \
              + np.outer(sk, sk) / np.inner(yk, sk)
    else:
        Tic = T0k
    return Tic
def Inversion_D_BFGS(T0, T1, kappa, u, Td, Tic0, epsilon, MAX, PATH):
    """
    Inverse heat conduction equation with the iteration of adjoint equation method.
    Minimize J(objective function wich represents the mismach between the prediction and the data)
        with BFGS Quasi-Newton Method.
    From Td back to 0.
        Td is the data obtained today.
        Tic0 is a initial guess of the iteration.
    """

    logFile = PATH + '/log.txt'
    for k in range(MAX):
        if k == 0:
            T0k = Tic0
            Vk = np.eye(globalVar.Nz + 1, dtype=np.float64)
            T1k = CN_D(T0=T0, T1=T1, kappa=kappa, u=u, Tic=T0k)[-1, :]
            lambda0k = CN_D_B(T0=0, T1=0, kappa=kappa, u=u, Tec=2 * (T1k - Td))[0, :]
            Jk = norm_2(T1k - Td)
        else:
            T0k = T0k1
            Vk = Vk1
            T1k = T1k1
            lambda0k = lambda0k1
            Jk = Jk1
        if not os.path.exists(PATH + '/T0'):
            os. mkdir(PATH + '/T0')
        if not os.path.exists(PATH + '/T1'):
            os. mkdir(PATH + '/T1')
        if k % 50 == 0:  # plot
            plt.ylim(-0.1, 1.1)
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), T0k, 'r-', label='Tick')
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), globalVar.Tic_real(), 'b-', label='Tic_real')
            no = '{:0>4}'.format(str(k))
            plt.xlabel('z')
            plt.ylabel('T0')
            plt.legend(loc='upper right')
            plt.savefig(PATH + '/T0/'+ no + 'T0.png')
            plt.cla()

            plt.ylim(-0.1, 1.1)
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), T1k, 'r-', label='T1n')
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), Td, 'b-', label='Td')
            no = '{:0>4}'.format(str(k))
            plt.xlabel('z')
            plt.ylabel('T1')
            plt.legend(loc='upper right')
            plt.savefig(PATH + '/T1/'+ no + 'T1.png')
            plt.cla()

        with open(logFile, 'a') as file:
            file.write('J{:<4} = {:<9.7}\n'.format(k, Jk))
        if Jk < epsilon:
            Tic = T0k
            break

        pk = - np.dot(Vk, lambda0k)
        alphak = 1  # TODO: line search to get a step length alphak.
        sk = alphak * pk
        T0k1 = T0k + sk
        if globalVar.set_T0_always_positive:
            for i in range(globalVar.Nz + 1):
                if T0k1[i] < 0:
                    T0k1[i] = 0
        T1k1 = CN_D(T0=T0, T1=T1, kappa=kappa, u=u, Tic=T0k1)[-1, :]
        Jk1 = norm_2(T1k1 - Td)
        while Jk1 > Jk:
            if alphak < 1e-3:
                sk = - lambda0k
                T0k1 = T0k + sk
                if globalVar.set_T0_always_positive:
                    for i in range(globalVar.Nz + 1):
                        if T0k1[i] < 0:
                            T0k1[i] = 0
                T1k1 = CN_D(T0=T0, T1=T1, kappa=kappa, u=u, Tic=T0k1)[-1, :]
                Jk1 = norm_2(T1k1 - Td)
                break
            else:
                alphak /= 2
                sk = alphak * pk
                T0k1 = T0k + sk
                if globalVar.set_T0_always_positive:
                    for i in range(globalVar.Nz + 1):
                        if T0k1[i] < 0:
                            T0k1[i] = 0
                T1k1 = CN_D(T0=T0, T1=T1, kappa=kappa, u=u, Tic=T0k1)[-1, :]
                Jk1 = norm_2(T1k1 - Td)
        if alphak != 1:
            print(alphak)
        lambda0k1 = CN_D_B(T0=0, T1=0, kappa=kappa, u=u, Tec=2 * (T1k1 - Td))[0, :]
        yk = lambda0k1 - lambda0k
        Vk1 = np.dot(np.dot( (np.eye(globalVar.Nz + 1) - np.outer(sk, yk)/np.inner(sk, yk)),
                             Vk ),
                     (np.eye(globalVar.Nz + 1) - np.outer(yk, sk) / np.inner(sk, yk))) \
              + np.outer(sk, sk) / np.inner(sk, yk)
    else:
        Tic = T0k
    return Tic
def Inversion_D_BFGS_root(T0, T1, kappa, u, Td, Tic0, epsilon, MAX, PATH):
    """
    Inverse heat conduction equation with the iteration of adjoint equation method.
    Minimize J(objective function wich represents the mismach between the prediction and the data)
        with BFGS Quasi-Newton Method.
    From Td back to 0.
        Td is the data obtained today.
        Tic0 is a initial guess of the iteration.
    """

    logFile = PATH + '/log.txt'
    for k in range(MAX):
        if k == 0:
            T0k = Tic0
            Vk = np.eye(globalVar.Nz + 1, dtype=np.float64)
            T1k = CN_D(T0=T0, T1=T1, kappa=kappa, u=u, Tic=T0k)[-1, :]
            Jk = norm_2(T1k - Td)
            lambda0k = CN_D_B(T0=0, T1=0, kappa=kappa, u=u, Tec=2 * (T1k - Td))[0, :]
            gk = lambda0k / 2 / Jk
        else:
            T0k = T0k1
            Vk = Vk1
            T1k = T1k1
            gk = gk1
            Jk = Jk1
        if not os.path.exists(PATH + '/T0'):
            os. mkdir(PATH + '/T0')
        if not os.path.exists(PATH + '/T1'):
            os. mkdir(PATH + '/T1')
        if k % 50 == 0:  # plot
            plt.ylim(-0.1, 1.1)
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), T0k, 'r-', label='Tick')
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), globalVar.Tic_real(), 'b-', label='Tic_real')
            no = '{:0>4}'.format(str(k))
            plt.xlabel('z')
            plt.ylabel('T0')
            plt.legend(loc='upper right')
            plt.savefig(PATH + '/T0/'+ no + 'T0.png')
            plt.cla()

            plt.ylim(-0.1, 1.1)
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), T1k, 'r-', label='T1n')
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), Td, 'b-', label='Td')
            no = '{:0>4}'.format(str(k))
            plt.xlabel('z')
            plt.ylabel('T1')
            plt.legend(loc='upper right')
            plt.savefig(PATH + '/T1/'+ no + 'T1.png')
            plt.cla()

        with open(logFile, 'a') as file:
            file.write('J{:<4} = {:<9.7}\n'.format(k, Jk))
        if Jk < epsilon:
            return T0k

        pk = - np.dot(Vk, gk)
        alphak = 0.5  # TODO: line search to get a step length alphak.
        sk = alphak * pk
        T0k1 = T0k + sk
        if globalVar.set_T0_always_positive:
            for i in range(globalVar.Nz + 1):
                if T0k1[i] < 0:
                    T0k1[i] = 0
        T1k1 = CN_D(T0=T0, T1=T1, kappa=kappa, u=u, Tic=T0k1)[-1, :]
        Jk1 = norm_2(T1k1 - Td)
        while Jk1 > Jk:
            if alphak < 1e-15:
                return T0k
            else:
                alphak /= 2
                sk = alphak * pk
                T0k1 = T0k + sk
                if globalVar.set_T0_always_positive:
                    for i in range(globalVar.Nz + 1):
                        if T0k1[i] < 0:
                            T0k1[i] = 0
                T1k1 = CN_D(T0=T0, T1=T1, kappa=kappa, u=u, Tic=T0k1)[-1, :]
                Jk1 = norm_2(T1k1 - Td)
        print(alphak)
        lambda0k1 = CN_D_B(T0=0, T1=0, kappa=kappa, u=u, Tec=2 * (T1k1 - Td))[0, :]
        gk1 = lambda0k1 / 2 / Jk1
        yk = gk1 - gk
        Vk1 = np.dot(np.dot( (np.eye(globalVar.Nz + 1) - np.outer(sk, yk)/np.inner(sk, yk)),
                             Vk ),
                     (np.eye(globalVar.Nz + 1) - np.outer(yk, sk) / np.inner(sk, yk))) \
              + np.outer(sk, sk) / np.inner(sk, yk)
    else:
        return T0k

# TODO: complete the N part.
def Inversion_N_BFGS_root(T0, p, kappa, u, Td, Tic0, epsilon, MAX, PATH):
    """
    Inverse heat conduction equation with the iteration of adjoint equation method
    from Td back to 0.
    Td is the data obtained today.
    Tic0 is a initial guess of the iteration.
    """

    logFile = PATH + '/log.txt'
    for k in range(MAX):
        if k == 0:
            T0k = Tic0
            Vk = np.eye(globalVar.Nz + 1, dtype=np.float64)
            T1k = CN_N(T0=T0, p=p, kappa=kappa, u=u, Tic=T0k)[-1, :]
            Jk = norm_2(T1k - Td)
            lambda0k = CN_N_B(T0=0, p=0, kappa=kappa, u=u, Tec=2 * (T1k - Td))[0, :]
            gk = lambda0k / 2 / Jk
        else:
            T0k = T0k1
            Vk = Vk1
            T1k = T1k1
            gk = gk1
            Jk = Jk1
        if not os.path.exists(PATH + '/T0'):
            os. mkdir(PATH + '/T0')
        if not os.path.exists(PATH + '/T1'):
            os. mkdir(PATH + '/T1')
        if k % 50 == 0:  # plot
            plt.ylim(-0.1, 1.1)
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), T0k, 'r-', label='Tick')
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), globalVar.Tic_real(), 'b-', label='Tic_real')
            no = '{:0>4}'.format(str(k))
            plt.xlabel('z')
            plt.ylabel('T0')
            plt.legend(loc='upper right')
            plt.savefig(PATH + '/T0/'+ no + 'T0.png')
            plt.cla()

            plt.ylim(-0.1, 1.1)
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), T1k, 'r-', label='T1n')
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), Td, 'b-', label='Td')
            no = '{:0>4}'.format(str(k))
            plt.xlabel('z')
            plt.ylabel('T1')
            plt.legend(loc='upper right')
            plt.savefig(PATH + '/T1/'+ no + 'T1.png')
            plt.cla()

        with open(logFile, 'a') as file:
            file.write('J{:<4} = {:<9.7}\n'.format(k, Jk))
        if Jk < epsilon:
            return T0k

        pk = - np.dot(Vk, gk)
        alphak = 0.5  # TODO: line search to get a step length alphak.
        sk = alphak * pk
        T0k1 = T0k + sk
        if globalVar.set_T0_always_positive:
            for i in range(globalVar.Nz + 1):
                if T0k1[i] < 0:
                    T0k1[i] = 0
        T1k1 = CN_N(T0=T0, p=p, kappa=kappa, u=u, Tic=T0k1)[-1, :]
        Jk1 = norm_2(T1k1 - Td)
        while Jk1 > Jk:
            if alphak < 1e-15:
                return T0k
            else:
                alphak /= 2
                sk = alphak * pk
                T0k1 = T0k + sk
                if globalVar.set_T0_always_positive:
                    for i in range(globalVar.Nz + 1):
                        if T0k1[i] < 0:
                            T0k1[i] = 0
                T1k1 = CN_N(T0=T0, p=p, kappa=kappa, u=u, Tic=T0k1)[-1, :]
                Jk1 = norm_2(T1k1 - Td)
        print(alphak)
        lambda0k1 = CN_N_B(T0=0, p=0, kappa=kappa, u=u, Tec=2 * (T1k1 - Td))[0, :]
        gk1 = lambda0k1 / 2 / Jk1
        yk = gk1 - gk
        Vk1 = np.dot(np.dot( (np.eye(globalVar.Nz + 1) - np.outer(sk, yk)/np.inner(sk, yk)),
                             Vk ),
                     (np.eye(globalVar.Nz + 1) - np.outer(yk, sk) / np.inner(sk, yk))) \
              + np.outer(sk, sk) / np.inner(sk, yk)
    else:
        return T0k

def main():
    n = 1
    PATH = 'case/' + 'case' + '{:0>3}'.format(n)
    while os.path.exists(PATH):
        n += 1
        PATH = 'case/' + 'case' + '{:0>3}'.format(n)
    os.mkdir(PATH)

    Inversion_method = Inversion_N_BFGS_root

    start_time = time.clock()

    # Td = CN_D(T0=globalVar.T0, T1=globalVar.T1, kappa=globalVar.kappa, u=globalVar.u(PATH, True), Tic=globalVar.Tic_real())[-1, :]
    # Tic = Inversion_method(T0=globalVar.T0, T1=globalVar.T1, kappa=globalVar.kappa, u=globalVar.u(),
    #                   Td = Td, Tic0=globalVar.Tic_guess(), epsilon=globalVar.epsilon, MAX=globalVar.MAX, PATH=PATH)

    Td = CN_N(T0=globalVar.T0, p=globalVar.p, kappa=globalVar.kappa, u=globalVar.u(PATH, True), Tic=globalVar.Tic_real())[-1, :]
    Tic = Inversion_method(T0=globalVar.T0, p=globalVar.p, kappa=globalVar.kappa, u=globalVar.u(),
                      Td = Td, Tic0=globalVar.Tic_guess(), epsilon=globalVar.epsilon, MAX=globalVar.MAX, PATH=PATH)

    # plot
    plt.plot(np.linspace(0, 1, globalVar.Nz + 1), Td, 'r-', label='Td')
    plt.plot(np.linspace(0, 1, globalVar.Nz + 1), globalVar.Tic_real(), 'g--', label='Tic_real')
    plt.plot(np.linspace(0, 1, globalVar.Nz + 1), Tic, 'b-', label='Tic')
    plt.plot(np.linspace(0, 1, globalVar.Nz + 1), globalVar.Tic_guess(), 'b--', label='Tic_guess')
    plt.xlabel('z')
    plt.ylabel('T')
    plt.text(0, 0.7, 't = {:<7}\nkappa = {:<7.5}\n'
                     'u = {:<}\nTic_real = {:<}\nTic_guess = {:<}\n'
                     'Set T0 always positive: {}\n'
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