from invent_data import *

# for p in [10, 20, 30]:
#     sh0 = (p - globalVar.pm) * globalVar.kappa / globalVar.hr
#     T = globalVar.Tic_continent(globalVar.pm, sh0, globalVar.hr)
#     plt.plot(T, np.linspace(0, globalVar.zTotal, globalVar.Nz + 1), label=f'{p}K/km')
# plt.gca().invert_yaxis()
# plt.xlabel('T / $^{\circ}$C')
# plt.ylabel('z / km')
# legend = plt.legend(loc='upper right')
# plt.title('Continental Geotherm Model')
# legend.get_frame().set_facecolor('none')
# # plt.show()
# plt.savefig('C:/Users/Administrator/Desktop/Centinental geotherm model.png', transparent=True)



# T = CN_N(Ts=0, p=globalVar.pm, kappa=globalVar.kappa, u=globalVar.u,
#          Tic=globalVar.Tic_continent(p=globalVar.pm, sh0=globalVar.sh0, hr=globalVar.hr),
#          sh=globalVar.sh_continental(u=globalVar.u, sh0=globalVar.sh0, hr=globalVar.hr))
# for n in np.arange(0, 51):
#     t = n * globalVar.deltat
#     plt.plot(T[n], np.linspace(0, globalVar.zTotal, globalVar.Nz + 1))
#     # plt.gca().invert_yaxis()
#     plt.title(f'{100-t}Ma')
#     plt.xlabel('T / $^{\circ}$C')
#     plt.ylabel('z / km')
#     plt.xlim(0, 2300)
#     plt.ylim(200, 0)
#     plt.savefig(f'C:/Users/Administrator/Desktop/forward pngs/{n}.png', transparent=True)
#     plt.cla()



T = CN_N(Ts=0, p=globalVar.pm, kappa=globalVar.kappa, u=globalVar.u,
         Tic=globalVar.Tic_continent(p=globalVar.pm, sh0=globalVar.sh0, hr=globalVar.hr),
         sh=globalVar.sh_continental(u=globalVar.u, sh0=globalVar.sh0, hr=globalVar.hr))
plt.plot(T[0, :], np.linspace(0, globalVar.zTotal, globalVar.Nz + 1),
         'r-', label='    0 Myr')
plt.plot(T[int(globalVar.Nt/2), :], np.linspace(0, globalVar.zTotal, globalVar.Nz + 1),
         'g-', label='{:>5}'.format(int(globalVar.tTotal/2)) +' Myr')
plt.plot(T[-1, :], np.linspace(0, globalVar.zTotal, globalVar.Nz + 1),
         'b-', label='{:>5}'.format(globalVar.tTotal) + ' Myr')
legend = plt.legend(loc='upper right')
plt.xlabel('T / $^{\circ}$C')
plt.ylabel('z / km')
plt.xlim(0, 2300)
plt.ylim(200, 0)
plt.title('Forward Modelling with Constant Velocity u = 0.3km/Myr')
legend.get_frame().set_facecolor('none')
plt.savefig('C:/Users/Administrator/Desktop/forward.png', transparent=True)




# PATH = 'C:/Users/Administrator/Desktop/test'
# Tp = CN_N(Ts=globalVar.Ts, p=globalVar.pm, kappa=globalVar.kappa, u=globalVar.u, Tic=globalVar.Tic_real,
#           sh=globalVar.sh_continental(sh0=globalVar.sh0, hr=globalVar.hr, u=globalVar.u))[-1, :]
# Tic = Inversion_N_BFGS_modified(Ts=globalVar.Ts, p=globalVar.pm, kappa=globalVar.kappa, u=globalVar.u,
#                        Tp=Tp, Tic0=globalVar.Tic_guess, epsilon=globalVar.epsilon, MAX=globalVar.MAX, PATH=PATH,
#                        sh=globalVar.sh_continental(sh0=globalVar.sh0, hr=globalVar.hr, u=globalVar.u))
# np.savetxt(PATH + '/T.txt', Tic, fmt='%10.5f')
#
# # plot
# plt.plot(Tp, np.linspace(0, globalVar.zTotal, globalVar.Nz + 1), 'r-', label='Tp')
# plt.plot(globalVar.Tic_real, np.linspace(0, globalVar.zTotal, globalVar.Nz + 1), 'g--', label='Tic_real')
# plt.plot(Tic, np.linspace(0, globalVar.zTotal, globalVar.Nz + 1), 'b-', label='Tic')
# plt.plot(globalVar.Tic_guess, np.linspace(0, globalVar.zTotal, globalVar.Nz + 1), 'b--', label='Tic_guess')
# plt.xlabel('T / $^{\circ}$C')
# plt.ylabel('z / km')
# plt.xlim(0, 2300)
# plt.ylim(200, 0)
# plt.legend(loc='upper right')
# plt.title('Inversion Result')
#
# plt.savefig(PATH + '/inversion.png')
# plot_log(PATH + '/log.txt', PATH)
#
# plot_velocity(globalVar.u, mode='time', PATH=PATH, tTotal=globalVar.tTotal)



# u = globalVar.u_gaussian_time[:, 0] * globalVar.u_mag
# plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nt + 1), u)
# plt.hlines(-globalVar.u_mag / 2, 0, globalVar.tTotal, colors='grey')
# plt.gca().yaxis.grid(ls='--')
# plt.xlim(0, globalVar.tTotal)
# plt.title('Gaussian like Velocity')
# plt.xlabel('t / Myr')
# plt.ylabel('u / km/Myr')
# plt.arrow(38, -0.15, 62-38, 0, length_includes_head=True, head_length=3, head_width=0.01, shape='full', color='k')
# plt.arrow(62, -0.15, 38-62, 0, length_includes_head=True, head_length=3, head_width=0.01, shape='full', color='k')
# plt.text(50, -0.14, r'$FWHM=2.355\sigma$', horizontalalignment='center')
# plt.savefig('C:/Users/Administrator/Desktop/gaussian_velocity.png', transparent=True)



# u = globalVar.u_gaussian_time
# T_real = CN_N(Ts=globalVar.Ts, p=globalVar.pm, kappa=globalVar.kappa, u=u, Tic=globalVar.Tic_real,
#               sh=globalVar.sh_continental(sh0=globalVar.sh0, hr=globalVar.hr, u=u))
# Td = get_thermal_history(T_real, u[:, 0])
# Tp = T_real[-1, :]
# ps = (Tp[1] - Tp[0]) / globalVar.deltaz
# sh0 = (ps - globalVar.pm) * globalVar.kappa / globalVar.hr
# T_reconstructed_by_qs = globalVar.Tic_continent(globalVar.pm, sh0, globalVar.hr)
# T_reconstructed_by_qs = np.tile(T_reconstructed_by_qs, (globalVar.Nt + 1, 1))
# u_r_qs = get_velocity_transverse(T_reconstructed_by_qs, Td)
# J_value = np.sqrt(scipy.integrate.simps(np.array(u_r_qs - u[:, 0]) ** 2, dx=globalVar.deltat) / globalVar.tTotal)
# plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nt + 1), u[:, 0], label='Real velocity')
# plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nt + 1), u_r_qs, label='Recovered velocity')
# plt.xlabel('t / Myr')
# plt.ylabel('u / km/Myr')
# legend = plt.legend()
# legend.get_frame().set_facecolor('none')
# plt.title('Recovered Veloctiy')
# plt.savefig('C:/Users/Administrator/Desktop/recover_gaussian_velocity.png', transparent=True)



# u = globalVar.u_gaussian_time[:, 0] * globalVar.u_mag
# for k in range(globalVar.Nt + 1):
#     tk = k * globalVar.deltat
#     if tk <= globalVar.tTotal - sigma:
#         u[k, :] = 0
#     else:
#         u[k, :] = - (tk - (globalVar.tTotal - sigma)) / sigma * u0
# plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nt + 1), u)
# plt.hlines(-globalVar.u_mag / 2, 0, globalVar.tTotal, colors='grey')
# plt.gca().yaxis.grid(ls='--')
# plt.xlim(0, globalVar.tTotal)
# plt.title('Gaussian like Velocity')
# plt.xlabel('t / Myr')
# plt.ylabel('u / km/Myr')
# plt.arrow(38, -0.15, 62-38, 0, length_includes_head=True, head_length=3, head_width=0.01, shape='full', color='k')
# plt.arrow(62, -0.15, 38-62, 0, length_includes_head=True, head_length=3, head_width=0.01, shape='full', color='k')
# plt.text(50, -0.14, r'$FWHM=2.355\sigma$', horizontalalignment='center')
# plt.savefig('C:/Users/Administrator/Desktop/gaussian_velocity.png', transparent=True)
