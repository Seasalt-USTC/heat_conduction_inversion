from forward import *
import os

# TODO: modify
def Inversion_D_Steepest(Ts, Tb, kappa, u, Td, Tic0, epsilon, MAX, PATH):
    """
    Inverse heat conduction equation with the iteration of adjoint equation method
    from Td back to 0.
    Td is the data obtained today.
    Tic0 is a initial guess of the iteration.
    """

    logFile = PATH + '/log.txt'
    T0k = Tic0
    for k in range(MAX):
        T1k = CN_D(Ts=Ts, Tb=Tb, kappa=kappa, u=u, Tic=T0k)[-1, :]

        if not os.path.exists(PATH + '/Ts'):
            os. mkdir(PATH + '/Ts')
        if not os.path.exists(PATH + '/Tb'):
            os. mkdir(PATH + '/Tb')
        if k % 50 == 0:
            plt.ylim(-0.1, 1.1)
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), T0k, 'r-', label='Tick')
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), globalVar.Tic_real(), 'b-', label='Tic_real')
            no = '{:0>4}'.format(str(k))
            plt.xlabel('z')
            plt.ylabel('Ts')
            plt.legend(loc='upper right')
            plt.savefig(PATH + '/Ts/'+ no + 'Ts.png')
            plt.cla()

            plt.ylim(-0.1, 1.1)
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), T1k, 'r-', label='T1n')
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), Td, 'b-', label='Td')
            no = '{:0>4}'.format(str(k))
            plt.xlabel('z')
            plt.ylabel('Tb')
            plt.legend(loc='upper right')
            plt.savefig(PATH + '/Tb/'+ no + 'Tb.png')
            plt.cla()

        Jk = norm_2(T1k - Td)
        with open(logFile, 'a') as file:
            file.write('J{:<4} = {:<9.7}\n'.format(k, Jk))
        if Jk <= epsilon:
            Tic = T0k
            break

        lambda0k = CN_D_B(Ts=0, Tb=0, kappa=kappa, u=u, Tec=2 * (T1k - Td))[0, :]
        alphak = 1
        T0k = T0k - alphak * lambda0k

        if globalVar.set_T0_always_positive:
            for i in range(globalVar.Nz + 1):
                if T0k[i] < 0:
                    T0k[i] = 0
    else:
        Tic = T0k
    return Tic
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
        if not os.path.exists(PATH + '/Ts'):
            os. mkdir(PATH + '/Ts')
        if not os.path.exists(PATH + '/Tb'):
            os. mkdir(PATH + '/Tb')
        if k % 50 == 0:  # plot
            plt.ylim(-0.1, 1.1)
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), T0k, 'r-', label='Tick')
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), globalVar.Tic_real(), 'b-', label='Tic_real')
            no = '{:0>4}'.format(str(k))
            plt.xlabel('z')
            plt.ylabel('Ts')
            plt.legend(loc='upper right')
            plt.savefig(PATH + '/Ts/'+ no + 'Ts.png')
            plt.cla()

            plt.ylim(-0.1, 1.1)
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), T1k, 'r-', label='T1n')
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), Td, 'b-', label='Td')
            no = '{:0>4}'.format(str(k))
            plt.xlabel('z')
            plt.ylabel('Tb')
            plt.legend(loc='upper right')
            plt.savefig(PATH + '/Tb/'+ no + 'Tb.png')
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
        T1k1 = CN_D(Ts=Ts, Tb=Tb, kappa=kappa, u=u, Tic=T0k1)[-1, :]
        # Jk1 = norm_2(T1k1 - Td)
        # while Jk1 > Jk:
        #     alphak /= 2
        #     sk = alphak * pk
        #     T0k1 = T0k + sk
        #     if globalVar.set_T0_always_positive:
        #         for i in range(globalVar.Nz + 1):
        #             if T0k1[i] < 0:
        #                 T0k1[i] = 0
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
        if not os.path.exists(PATH + '/Ts'):
            os. mkdir(PATH + '/Ts')
        if not os.path.exists(PATH + '/Tb'):
            os. mkdir(PATH + '/Tb')
        if k % 50 == 0:  # plot
            plt.ylim(-0.1, 1.1)
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), T0k, 'r-', label='Tick')
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), globalVar.Tic_real(), 'b-', label='Tic_real')
            no = '{:0>4}'.format(str(k))
            plt.xlabel('z')
            plt.ylabel('Ts')
            plt.legend(loc='upper right')
            plt.savefig(PATH + '/Ts/'+ no + 'Ts.png')
            plt.cla()

            plt.ylim(-0.1, 1.1)
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), T1k, 'r-', label='T1n')
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), Td, 'b-', label='Td')
            no = '{:0>4}'.format(str(k))
            plt.xlabel('z')
            plt.ylabel('Tb')
            plt.legend(loc='upper right')
            plt.savefig(PATH + '/Tb/'+ no + 'Tb.png')
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
        T1k1 = CN_D(Ts=Ts, Tb=Tb, kappa=kappa, u=u, Tic=T0k1)[-1, :]
        Jk1 = norm_2(T1k1 - Td)
        while Jk1 > Jk:
            if alphak < 1e-3:
                sk = - lambda0k
                T0k1 = T0k + sk
                if globalVar.set_T0_always_positive:
                    for i in range(globalVar.Nz + 1):
                        if T0k1[i] < 0:
                            T0k1[i] = 0
                T1k1 = CN_D(Ts=Ts, Tb=Tb, kappa=kappa, u=u, Tic=T0k1)[-1, :]
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
        if not os.path.exists(PATH + '/Ts'):
            os. mkdir(PATH + '/Ts')
        if not os.path.exists(PATH + '/Tb'):
            os. mkdir(PATH + '/Tb')
        if k % 50 == 0:  # plot
            plt.ylim(-0.1, 1.1)
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), T0k, 'r-', label='Tick')
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), globalVar.Tic_real(), 'b-', label='Tic_real')
            no = '{:0>4}'.format(str(k))
            plt.xlabel('z')
            plt.ylabel('Ts')
            plt.legend(loc='upper right')
            plt.savefig(PATH + '/Ts/'+ no + 'Ts.png')
            plt.cla()

            plt.ylim(-0.1, 1.1)
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), T1k, 'r-', label='T1n')
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), Td, 'b-', label='Td')
            no = '{:0>4}'.format(str(k))
            plt.xlabel('z')
            plt.ylabel('Tb')
            plt.legend(loc='upper right')
            plt.savefig(PATH + '/Tb/'+ no + 'Tb.png')
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
        T1k1 = CN_D(Ts=Ts, Tb=Tb, kappa=kappa, u=u, Tic=T0k1)[-1, :]
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
        return T0k

# TODO: complete the N part.
def Inversion_N_BFGS_root(Ts, p, kappa, u, Td, Tic0, epsilon, MAX, PATH):
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
        if not os.path.exists(PATH + '/Ts'):
            os. mkdir(PATH + '/Ts')
        if not os.path.exists(PATH + '/Tb'):
            os. mkdir(PATH + '/Tb')
        if k % 50 == 0:  # plot
            plt.ylim(-0.1, 1.1)
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), T0k, 'r-', label='Tick')
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), globalVar.Tic_real(), 'b-', label='Tic_real')
            no = '{:0>4}'.format(str(k))
            plt.xlabel('z')
            plt.ylabel('Ts')
            plt.legend(loc='upper right')
            plt.savefig(PATH + '/Ts/'+ no + 'Ts.png')
            plt.cla()

            plt.ylim(-0.1, 1.1)
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), T1k, 'r-', label='T1n')
            plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nz + 1), Td, 'b-', label='Td')
            no = '{:0>4}'.format(str(k))
            plt.xlabel('z')
            plt.ylabel('Tb')
            plt.legend(loc='upper right')
            plt.savefig(PATH + '/Tb/'+ no + 'Tb.png')
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
        T1k1 = CN_N(Ts=Ts, p=p, kappa=kappa, u=u, Tic=T0k1)[-1, :]
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
                T1k1 = CN_N(Ts=Ts, p=p, kappa=kappa, u=u, Tic=T0k1)[-1, :]
                Jk1 = norm_2(T1k1 - Td)
        print(alphak)
        lambda0k1 = CN_N_B(Ts=0, p=0, kappa=kappa, u=u, Tec=2 * (T1k1 - Td))[0, :]
        gk1 = lambda0k1 / 2 / Jk1
        yk = gk1 - gk
        Vk1 = np.dot(np.dot( (np.eye(globalVar.Nz + 1) - np.outer(sk, yk)/np.inner(sk, yk)),
                             Vk ),
                     (np.eye(globalVar.Nz + 1) - np.outer(yk, sk) / np.inner(sk, yk))) \
              + np.outer(sk, sk) / np.inner(sk, yk)
    else:
        return T0k
