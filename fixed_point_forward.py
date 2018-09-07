from invent_data import *
import globalVar
import matplotlib.pyplot as plt

MAX = 50

PATH = 'fixed-point_forward/from qs recovered velocity'

Tic = globalVar.Tic_real
T_real = CN_N(Ts=globalVar.Ts, p=globalVar.pm, kappa=globalVar.kappa, u=globalVar.u, Tic=Tic,
         sh=globalVar.sh_continental(sh0=globalVar.sh0, hr=globalVar.hr, u=globalVar.u))
Td = get_thermal_history(T_real, globalVar.u[:, 0])
np.savetxt(PATH + '/Td.txt', Td)
Tp = T_real[-1, :]
ps = (Tp[1] - Tp[0]) / globalVar.deltaz
np.savetxt(PATH + '/ps.txt', [ps])

Td = np.loadtxt(PATH + '/Td.txt')
ps = np.loadtxt(PATH + '/ps.txt')

# Construct a geotherm from surface heat flow.
# Use this geotherm with the assumption that (partial T/partial t) == 0 to recover the velocity.
sh0 = (ps - globalVar.pm) * globalVar.kappa / globalVar.hr
T_reconstructed_by_qs =  globalVar.Tic_continent(globalVar.pm, sh0, globalVar.hr)
T_reconstructed_by_qs = np.tile(T_reconstructed_by_qs, (globalVar.Nt + 1, 1))

# u_guess = globalVar.u
# u = u_guess[:, 0]

u_guess = get_velocity_transverse(T=T_reconstructed_by_qs, Td=Td)
u = u_guess

counter = 0

while True:
    plt.plot(u, label='recovered velocity')
    plt.plot(globalVar.u[:, 0], label='real velocity')
    plt.legend()
    plt.savefig(PATH + '/u{:0>3}.png'.format(counter))
    plt.cla()

    u_2D = np.tile(np.transpose(np.expand_dims(u, 0)), globalVar.Nz + 1)
    T = CN_N(Ts=globalVar.Ts, p=globalVar.pm, kappa=globalVar.kappa, u=u_2D, Tic=Tic,
             sh=globalVar.sh_continental(sh0=globalVar.sh0, hr=globalVar.hr, u=u_2D))

    plt.plot(get_thermal_history(T, u), label='recovered thermal history')
    plt.plot(Td, label='real thermal history')
    plt.legend()
    plt.savefig(PATH + '/Td{:0>3}.png'.format(counter))
    plt.cla()

    u = get_velocity_transverse(T, Td)
    counter += 1

    if counter >= MAX:
        break