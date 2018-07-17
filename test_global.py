from inversion import *
# PATH = 'case/case001'
# plot_log(PATH + '/log.txt', PATH)

# plt.cla()
# logfiles = ['case/case001/log.txt', 'case/case005/log.txt', 'case/case006/log.txt', 'case/case007/log.txt']
# PATH = 'case'
# for logfile in logfiles:
#     J = []
#     with open(logfile) as f:
#         for line in f.readlines()[0: -1]:
#             J.append(float(line.split()[2]))
#     plt.loglog(J, label='J')
#     plt.xlabel('k')
#     plt.ylabel('J')
#     plt.savefig(PATH + '/J_loglog.png')

# plt.cla()
# PATH = 'case/case002'
# plot_log(PATH + '/log.txt', PATH)

# T = CN_D(Ts=0, Tb=0, kappa=globalVar.kappa, u=globalVar.u(), Tic=globalVar.Tic_real())
# TODO
T = CN_N(Ts=0, p=globalVar.p, kappa=globalVar.kappa, u=globalVar.u_uniform * (-0), Tic=globalVar.Tic_continent(hr=10),
         sh=globalVar.sh_continental(u=globalVar.u_uniform * (-0), hr=10))
# T = CN_D_B(Ts=0, Tb=0, kappa=globalVar.kappa, u=globalVar.u, Tec=Tic2())
# T = CN_N_B(Ts=0, p=0, kappa=globalVar.kappa, u=globalVar.u, Tec=Tic2())

# outputnStep = int(globalVar.outputTimeStep / globalVar.deltat)
#
# outputiStep = int(globalVar.outputSpaceStep / globalVar.deltaz)
#
# with open('T.xyz', 'w') as f:
#     for n in range(0, globalVar.Nt + 1, outputnStep):
#         for i in range(0, globalVar.Nz + 1, outputiStep):
#             f.write('{:10f} {:10f} {:10f}\n'.format(n * globalVar.deltat, i * globalVar.deltaz, T[n, i]))

qs = (T[:, 1] - T[:, 0]) / globalVar.deltaz * 3.35

plt.cla()
# plt.xlim(0, 40)
# plt.ylim(0, 600)
# plt.plot(np.linspace(0, globalVar.zTotal, globalVar.Nz + 1), T[0, :], 'r-',
#          np.linspace(0, globalVar.zTotal, globalVar.Nz + 1), T[int(globalVar.Nt/2), :], 'g-',
#          np.linspace(0, globalVar.zTotal, globalVar.Nz + 1), T[-1, :], 'b-',
#          )
plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nt + 1), qs, 'r-', label='qs')
plt.legend()
plt.show()
