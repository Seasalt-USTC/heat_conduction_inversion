from forward import *
from utils import *
import matplotlib.pyplot as plt
import os

def Inversion_D_Steepest(Ts, Tb, kappa, u, Td, Tic0, epsilon, MAX, PATH):
    """
    Inverse heat conduction equation with the iteration of adjoint equation method
    from Td back to 0.
    Td is the data obtained today.
    Tic0 is a initial guess of the iteration.
    """

    logFile = PATH + '/log.txt'
    T0k = Tic0
    with open(logFile, 'a') as log:
        for k in range(MAX):
            if k == 0:
                T1k = CN_D(Ts=Ts, Tb=Tb, kappa=kappa, u=u, Tic=T0k)[-1, :]
            else:
                T0k = T0k1
                T1k = T1k1
                Jk = Jk1
            if not os.path.exists(PATH + '/T0'):
                os.mkdir(PATH + '/T0')
            if not os.path.exists(PATH + '/T1'):
                os.mkdir(PATH + '/T1')
            if k % 50 == 0:  # plot
                plt.ylim(-0.1, 1.1)
                plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), T0k, 'r-', label='Tick')
                plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), globalVar.Tic_real(), 'b-',
                         label='Tic_real')
                no = '{:0>4}'.format(str(k))
                plt.xlabel('z')
                plt.ylabel('Ts')
                plt.legend(loc='upper right')
                plt.savefig(PATH + '/T0/' + no + 'T0.png')
                plt.cla()

                plt.ylim(-0.1, 1.1)
                plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), T1k, 'r-', label='T1n')
                plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), Td, 'b-', label='Td')
                no = '{:0>4}'.format(str(k))
                plt.xlabel('z')
                plt.ylabel('Tb')
                plt.legend(loc='upper right')
                plt.savefig(PATH + '/T1/' + no + 'T1.png')
                plt.cla()

            log.write('J{:<4} = {:<9.7}\n'.format(k, Jk))  # write log
            if Jk <= epsilon:
                log.write('Return: J is lower than epsilon.')
                return T0k

            lambda0k = CN_D_B(Ts=0, Tb=0, kappa=kappa, u=u, Tec=2 * (T1k - Td))[0, :]
            alphak = 1
            T0k1 = T0k - alphak * lambda0k
            T1k1 = CN_D(Ts=Ts, Tb=Tb, kappa=kappa, u=u, Tic=T0k1)[-1, :]
            Jk1 = norm_2(T1k1 - Td)
            while Jk1 > Jk:
                if not globalVar.line_search:
                    break
                if alphak < 1e-15:
                    log.write('Return: alphak is smaller than 1e-15.')
                    return T0k1
                alphak /= 2
                T0k1 = T0k - alphak * lambda0k
                T1k1 = CN_D(Ts=Ts, Tb=Tb, kappa=kappa, u=u, Tic=T0k1)[-1, :]
                Jk1 = norm_2(T1k1 - Td)

        else:
            log.write('Return: max iterations')
            return T0k
def Inversion_D_DFP(Ts, Tb, kappa, u, Td, Tic0, epsilon, MAX, PATH):
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
            T1k = CN_D(Ts=Ts, Tb=Tb, kappa=kappa, u=u, Tic=T0k)[-1, :]
            lambda0k = CN_D_B(Ts=0, Tb=0, kappa=kappa, u=u, Tec=2 * (T1k - Td))[0, :]
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
            plt.ylabel('Ts')
            plt.legend(loc='upper right')
            plt.savefig(PATH + '/T0/'+ no + 'T0.png')
            plt.cla()

            plt.ylim(-0.1, 1.1)
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), T1k, 'r-', label='T1n')
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), Td, 'b-', label='Td')
            no = '{:0>4}'.format(str(k))
            plt.xlabel('z')
            plt.ylabel('Tb')
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
        alphak = 1
        sk = alphak * pk
        T0k1 = T0k + sk
        T1k1 = CN_D(Ts=Ts, Tb=Tb, kappa=kappa, u=u, Tic=T0k1)[-1, :]
        # Jk1 = norm_2(T1k1 - Td)
        # while Jk1 > Jk:
        #     alphak /= 2
        #     sk = alphak * pk
        #     T0k1 = T0k + sk
        #     T1k1 = CN_D(Ts=Ts, Tb=Tb, kappa=kappa, u=u, Tic=T0k1)[-1, :]
        #     Jk1 = norm_2(T1k1 - Td)
        # print(alphak)
        lambda0k1 = CN_D_B(Ts=0, Tb=0, kappa=kappa, u=u, Tec=2 * (T1k1 - Td))[0, :]
        yk = lambda0k1 - lambda0k
        Vk1 = Vk \
              - np.dot( np.outer( np.dot(Vk, yk), yk), Vk) / np.inner(np.dot(yk, Vk), yk) \
              + np.outer(sk, sk) / np.inner(yk, sk)
    else:
        Tic = T0k
    return Tic
def Inversion_D_BFGS(Ts, Tb, kappa, u, Td, Tic0, epsilon, MAX, PATH):
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
            T1k = CN_D(Ts=Ts, Tb=Tb, kappa=kappa, u=u, Tic=T0k)[-1, :]
            lambda0k = CN_D_B(Ts=0, Tb=0, kappa=kappa, u=u, Tec=2 * (T1k - Td))[0, :]
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
            plt.ylabel('Ts')
            plt.legend(loc='upper right')
            plt.savefig(PATH + '/T0/'+ no + 'T0.png')
            plt.cla()

            plt.ylim(-0.1, 1.1)
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), T1k, 'r-', label='T1n')
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), Td, 'b-', label='Td')
            no = '{:0>4}'.format(str(k))
            plt.xlabel('z')
            plt.ylabel('Tb')
            plt.legend(loc='upper right')
            plt.savefig(PATH + '/T1/'+ no + 'T1.png')
            plt.cla()

        with open(logFile, 'a') as file:
            file.write('J{:<4} = {:<9.7}\n'.format(k, Jk))
        if Jk < epsilon:
            Tic = T0k
            break

        pk = - np.dot(Vk, lambda0k)
        alphak = 1
        sk = alphak * pk
        T0k1 = T0k + sk
        T1k1 = CN_D(Ts=Ts, Tb=Tb, kappa=kappa, u=u, Tic=T0k1)[-1, :]
        Jk1 = norm_2(T1k1 - Td)
        while Jk1 > Jk:
            if alphak < 1e-3:
                sk = - lambda0k
                T0k1 = T0k + sk
                T1k1 = CN_D(Ts=Ts, Tb=Tb, kappa=kappa, u=u, Tic=T0k1)[-1, :]
                Jk1 = norm_2(T1k1 - Td)
                break
            else:
                alphak /= 2
                sk = alphak * pk
                T0k1 = T0k + sk
                T1k1 = CN_D(Ts=Ts, Tb=Tb, kappa=kappa, u=u, Tic=T0k1)[-1, :]
                Jk1 = norm_2(T1k1 - Td)
        if alphak != 1:
            print(alphak)
        lambda0k1 = CN_D_B(Ts=0, Tb=0, kappa=kappa, u=u, Tec=2 * (T1k1 - Td))[0, :]
        yk = lambda0k1 - lambda0k
        Vk1 = np.dot(np.dot( (np.eye(globalVar.Nz + 1) - np.outer(sk, yk)/np.inner(sk, yk)),
                             Vk ),
                     (np.eye(globalVar.Nz + 1) - np.outer(yk, sk) / np.inner(sk, yk))) \
              + np.outer(sk, sk) / np.inner(sk, yk)
    else:
        Tic = T0k
    return Tic
def Inversion_D_BFGS_root(Ts, Tb, kappa, u, Td, Tic0, epsilon, MAX, PATH):
    """
    Inverse heat conduction equation with the iteration of adjoint equation method.
    Minimize J(objective function wich represents the mismach between the prediction and the data)
        with BFGS Quasi-Newton Method.
    From Td back to 0.
        Td is the data obtained today.
        Tic0 is a initial guess of the iteration.
    """

    logFile = PATH + '/log.txt'
    with open(logFile, 'a') as log:
        for k in range(MAX):
            if k == 0:
                T0k = Tic0
                Vk = np.eye(globalVar.Nz + 1, dtype=np.float64)
                T1k = CN_D(Ts=Ts, Tb=Tb, kappa=kappa, u=u, Tic=T0k)[-1, :]
                Jk = norm_2(T1k - Td)
                lambda0k = CN_D_B(Ts=0, Tb=0, kappa=kappa, u=u, Tec=2 * (T1k - Td))[0, :]
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
                plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), globalVar.Tic_real, 'b-', label='Tic_real')
                no = '{:0>4}'.format(str(k))
                plt.xlabel('z')
                plt.ylabel('Ts')
                plt.legend(loc='upper right')
                plt.savefig(PATH + '/T0/'+ no + 'T0.png')
                plt.cla()

                plt.ylim(-0.1, 1.1)
                plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), T1k, 'r-', label='T1n')
                plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), Td, 'b-', label='Td')
                no = '{:0>4}'.format(str(k))
                plt.xlabel('z')
                plt.ylabel('Tb')
                plt.legend(loc='upper right')
                plt.savefig(PATH + '/T1/'+ no + 'T1.png')
                plt.cla()

            log.write('J{:<4} = {:<9.7}\n'.format(k, Jk))  # write log
            if Jk < epsilon:
                log.write('Return: J is lower than epsilon.')
                return T0k

            pk = - np.dot(Vk, gk)
            alphak = 0.5  # start step length of th eline search.
            sk = alphak * pk
            T0k1 = T0k + sk
            T1k1 = CN_D(Ts=Ts, Tb=Tb, kappa=kappa, u=u, Tic=T0k1)[-1, :]
            Jk1 = norm_2(T1k1 - Td)

            # line search
            while Jk1 > Jk:
                if not globalVar.line_search:
                    break
                if alphak < 1e-15:
                    print('Return: alphak is smaller than 1e-15.')
                    return T0k
                else:
                    alphak /= 2
                    sk = alphak * pk
                    T0k1 = T0k + sk
                    T1k1 = CN_D(Ts=Ts, Tb=Tb, kappa=kappa, u=u, Tic=T0k1)[-1, :]
                    Jk1 = norm_2(T1k1 - Td)
            print(alphak)
            lambda0k1 = CN_D_B(Ts=0, Tb=0, kappa=kappa, u=u, Tec=2 * (T1k1 - Td))[0, :]
            gk1 = lambda0k1 / 2 / Jk1
            yk = gk1 - gk
            Vk1 = np.dot(np.dot( (np.eye(globalVar.Nz + 1) - np.outer(sk, yk)/np.inner(sk, yk)),
                                 Vk ),
                         (np.eye(globalVar.Nz + 1) - np.outer(yk, sk) / np.inner(sk, yk))) \
                  + np.outer(sk, sk) / np.inner(sk, yk)
        else:
            log.write('Return: max iterations')
            return T0k

def Inversion_N_Steepest(Ts, p, kappa, u, Td, Tic0, epsilon, MAX, PATH, sh=np.zeros((globalVar.Nt+1, globalVar.Nz+1))):
    """
    Inverse heat conduction equation with the iteration of adjoint equation method
    from Td back to 0.
    Td is the data obtained today.
    Tic0 is a initial guess of the iteration.
    """

    logFile = PATH + '/log.txt'
    T0k = Tic0
    with open(logFile, 'a') as log:
        for k in range(MAX):
            if k == 0:
                T1k = CN_N(Ts=Ts, p=p, kappa=kappa, u=u, Tic=T0k, sh=sh)[-1, :]
                Jk = norm_2(T1k - Td)
            else:
                T0k = T0k1
                T1k = T1k1
                Jk = Jk1
            if not os.path.exists(PATH + '/T0'):
                os. mkdir(PATH + '/T0')
            if not os.path.exists(PATH + '/T1'):
                os. mkdir(PATH + '/T1')
            if k % 50 == 0:  # plot
                # plt.ylim(-0.1, 1.1)
                plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), T0k, 'r-', label='Tick')
                plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), globalVar.Tic_real, 'b-', label='Tic_real')
                no = '{:0>4}'.format(str(k))
                plt.xlabel('z')
                plt.ylabel('Ts')
                plt.legend(loc='upper right')
                plt.savefig(PATH + '/T0/'+ no + 'T0.png')
                plt.cla()

                # plt.ylim(-0.1, 1.1)
                plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), T1k, 'r-', label='T1n')
                plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), Td, 'b-', label='Td')
                no = '{:0>4}'.format(str(k))
                plt.xlabel('z')
                plt.ylabel('Tb')
                plt.legend(loc='upper right')
                plt.savefig(PATH + '/T1/'+ no + 'T1.png')
                plt.cla()

            log.write('J{:<4} = {:<9.7}\n'.format(k, Jk))  # write log
            if Jk <= epsilon:
                log.write('Return: J is lower than epsilon.\n')
                return T0k

            lambda0k = CN_N_B(Ts=0, p=0, kappa=kappa, u=u, Tec=2 * (T1k - Td))[0, :]
            alphak = 1
            T0k1 = T0k - alphak * lambda0k
            T1k1 = CN_N(Ts=Ts, p=p, kappa=kappa, u=u, Tic=T0k1, sh=sh)[-1, :]
            Jk1 = norm_2(T1k1 - Td)
            while Jk1 > Jk:
                if not globalVar.line_search:
                    break
                alphak /= 2
                T0k1 = T0k - alphak * lambda0k
                T1k1 = CN_N(Ts=Ts, p=p, kappa=kappa, u=u, Tic=T0k1, sh=sh)[-1, :]
                Jk1 = norm_2(T1k1 - Td)
            print(alphak)
        else:
            log.write('Return: max iterations.\n')
            return T0k
def Inversion_N_BFGS_root(Ts, p, kappa, u, Td, Tic0, epsilon, MAX, PATH, sh=np.zeros((globalVar.Nt+1, globalVar.Nz+1))):
    """
    Inverse heat conduction equation with the iteration of adjoint equation method
    from Td back to 0.
    Td is the data obtained today.
    Tic0 is a initial guess of the iteration.
    """

    logFile = PATH + '/log.txt'

    if not os.path.exists(PATH + '/T0'):
        os.mkdir(PATH + '/T0')
    if not os.path.exists(PATH + '/T1'):
        os.mkdir(PATH + '/T1')


    np.seterr(all='raise')

    with open(logFile, 'a') as log:
        for k in range(MAX):
            if k == 0:
                T0k = Tic0
                Vk = np.eye(globalVar.Nz + 1, dtype=np.float64)
                T1k = CN_N(Ts=Ts, p=p, kappa=kappa, u=u, Tic=T0k)[-1, :]
                Jk = norm_2(T1k - Td)
                lambda0k = CN_N_B(Ts=0, p=0, kappa=kappa, u=u, Tec=2 * (T1k - Td))[0, :]
                gk = lambda0k / 2 / Jk
            else:
                T0k = T0k1
                Vk = Vk1
                T1k = T1k1
                gk = gk1
                Jk = Jk1
            if k % 50 == 0:  # plot
                plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), T0k, 'r-', label='Tick')
                plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), globalVar.Tic_real, 'b-', label='Tic_real')
                no = '{:0>4}'.format(str(k))
                plt.xlabel('z')
                plt.ylabel('Ts')
                plt.legend(loc='upper right')
                plt.savefig(PATH + '/T0/'+ no + 'T0.png')
                plt.cla()

                plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), T1k, 'r-', label='T1n')
                plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), Td, 'b-', label='Td')
                no = '{:0>4}'.format(str(k))
                plt.xlabel('z')
                plt.ylabel('Tb')
                plt.legend(loc='upper right')
                plt.savefig(PATH + '/T1/'+ no + 'T1.png')
                plt.cla()

            log.write('J{:<4} = {:<9.7}\n'.format(k, Jk))  # write log

            if Jk < epsilon:
                log.write('Return: J is lower than epsilon.\n')
                return T0k

            pk = - np.dot(Vk, gk)
            alphak = 0.5  # start step length of th eline search
            sk = alphak * pk
            T0k1 = T0k + sk
            T1k1 = CN_N(Ts=Ts, p=p, kappa=kappa, u=u, Tic=T0k1, sh=sh)[-1, :]
            Jk1 = norm_2(T1k1 - Td)

            while Jk1 > Jk:  # line search
                if not globalVar.line_search:
                    break
                alphak /= 2
                sk = alphak * pk
                T0k1 = T0k + sk
                T1k1 = CN_N(Ts=Ts, p=p, kappa=kappa, u=u, Tic=T0k1, sh=sh)[-1, :]
                Jk1 = norm_2(T1k1 - Td)
            print(alphak)
            lambda0k1 = CN_N_B(Ts=0, p=0, kappa=kappa, u=u, Tec=2 * (T1k1 - Td))[0, :]
            gk1 = lambda0k1 / 2 / Jk1
            yk = gk1 - gk
            try:
                Vk1 = np.dot(np.dot( (np.eye(globalVar.Nz + 1) - np.outer(sk, yk)/np.inner(sk, yk)), Vk ), (np.eye(globalVar.Nz + 1) - np.outer(yk, sk) / np.inner(sk, yk))) + np.outer(sk, sk) / np.inner(sk, yk)
            except Exception:
                log.write('Return: alphak is too small.\n')
                return T0k
        else:
            log.write('Return: max iterations.\n')
            return T0k
