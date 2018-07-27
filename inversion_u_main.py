from invent_data import *

def main():
    u = globalVar.u_gaussian_time
    Tic_real = globalVar.Tic_continent(p=globalVar.pm, sh0=globalVar.sh0, hr=globalVar.hr)

    T = CN_N(Ts=globalVar.Ts, p=globalVar.pm, u=u, Tic=Tic_real, kappa=globalVar.kappa,
             sh=globalVar.sh_continental(globalVar.sh0, globalVar.hr, u))
    Tp = T[-1, :]
    Td = get_thermal_history(T, u[:, 0])

    for countor in range(60, 80 + 1):
        if countor == 60: # initialize
            # Tp_exp = np.transpose(np.tile(np.transpose(np.expand_dims(Tp, 0)), globalVar.Nt + 1))
            # u_guess = get_velocity_transverse(T=Tp_exp, Td=Td)
            u_guess = np.loadtxt('temp/u060.txt')
        else:
            u_guess = get_velocity_transverse(T=T, Td=Td)

        Tic = Inversion_N_BFGS_modified(Ts=globalVar.Ts, p=globalVar.pm, kappa=globalVar.kappa,
                                        u=np.tile(np.transpose(np.expand_dims(u_guess, 0)), globalVar.Nz + 1),
                                        Td = Tp, Tic0=globalVar.Tic_guess, epsilon=globalVar.epsilon, MAX=globalVar.MAX, PATH='temp',
                                        sh=globalVar.sh_continental(sh0=globalVar.sh0, hr=globalVar.hr,
                                                                    u=np.tile(np.transpose(np.expand_dims(u_guess, 0)), globalVar.Nz + 1)))
        T = CN_N(Ts=globalVar.Ts, p=globalVar.pm,
                 u=np.tile(np.transpose(np.expand_dims(u_guess, 0)), globalVar.Nz + 1),
                 Tic=Tic, kappa=globalVar.kappa,
                 sh=globalVar.sh_continental(sh0=globalVar.sh0, hr=globalVar.hr,
                                             u=np.tile(np.transpose(np.expand_dims(u_guess, 0)), globalVar.Nz + 1)))
        Td_guess = get_thermal_history(T=T, u=u_guess)
        plt.plot(u_guess, label='u_guess')
        plt.plot(u[:, 0], label='u')
        plt.legend()
        plt.savefig('temp' + '/u{:0>3}.png'.format(countor))
        plt.cla()
        plt.plot(Td_guess, label='Td_guess')
        plt.plot(Td, label='Td')
        plt.legend()
        plt.savefig('temp' + '/Td{:0>3}.png'.format(countor))
        plt.cla()

        np.savetxt('temp' + '/u{:0>3}.txt'.format(countor), u_guess)
        pass

if __name__ == '__main__':
    main()