from invent_data import *
import globalVar
import matplotlib.pyplot as plt

Inversion_method = Inversion_N_BFGS_modified

PATH = 'fixed-point/from real velocity'

T_real = CN_N(Ts=globalVar.Ts, p=globalVar.pm, kappa=globalVar.kappa, u=globalVar.u, Tic=globalVar.Tic_real,
              sh=globalVar.sh_continental(sh0=globalVar.sh0, hr=globalVar.hr, u=globalVar.u))
Tp = T_real[-1, :]
np.savetxt(PATH + '/Tp.txt', Tp)

Td = get_thermal_history(T_real, globalVar.u[:, 0])
np.savetxt(PATH + '/Td.txt', Td)

Tp = np.loadtxt(PATH + '/Tp.txt')
Td = np.loadtxt(PATH + '/Td.txt')
u_guess = globalVar.u
u = u_guess[:, 0]

# Tp_exp = np.transpose(np.tile(np.transpose(np.expand_dims(Tp, 0)), globalVar.Nt + 1))
# u_guess = get_velocity_transverse(T=Tp_exp, Td=Td)
# u = u_guess

counter = 0

while True:
    plt.plot(u, label='recovered velocity')
    plt.plot(globalVar.u[:, 0], label='real velocity')
    plt.legend()
    plt.savefig(PATH + '/u{:0>3}.png'.format(counter))
    plt.cla()

    u_2D = np.tile(np.transpose(np.expand_dims(u, 0)), globalVar.Nz + 1)
    Tic = Inversion_method(Ts=globalVar.Ts, p=globalVar.pm, kappa=globalVar.kappa, u=u_2D,
                           Tp=Tp, Tic0=globalVar.Tic_guess, epsilon=globalVar.epsilon, MAX=globalVar.MAX, PATH=PATH,
                           sh=globalVar.sh_continental(sh0=globalVar.sh0, hr=globalVar.hr, u=u_2D))
    T = CN_N(Ts=globalVar.Ts, p=globalVar.pm, kappa=globalVar.kappa, u=u_2D, Tic=Tic,
             sh=globalVar.sh_continental(sh0=globalVar.sh0, hr=globalVar.hr, u=u_2D))

    plt.plot(get_thermal_history(T, u), label='recovered thermalhistory')
    plt.plot(Td, label='real thermal history')
    plt.legend()
    plt.savefig(PATH + '/Td{:0>3}.png'.format(counter))
    plt.cla()

    u = get_velocity_transverse(T, Td)
    counter += 1

    if counter >= 30:
        break

# plt.plot(Td)