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

# T = CN_D(T0=0, T1=0, kappa=globalVar.kappa, u=globalVar.u(), Tic=globalVar.Tic_real())
T = CN_N(T0=0, p=1, kappa=globalVar.kappa, u=globalVar.u0(), Tic=globalVar.Tic11(),
         q=np.ones((globalVar.Nt+1, globalVar.Nz+1), dtype=np.float64) * 0.005)
# T = CN_D_B(T0=0, T1=0, kappa=globalVar.kappa, u=globalVar.u, Tec=Tic2())
# T = CN_N_B(T0=0, p=0, kappa=globalVar.kappa, u=globalVar.u, Tec=Tic2())

# outputnStep = int(globalVar.outputTimeStep / globalVar.deltat)
#
# outputiStep = int(globalVar.outputSpaceStep / globalVar.deltaz)
#
# with open('T.xyz', 'w') as f:
#     for n in range(0, globalVar.Nt + 1, outputnStep):
#         for i in range(0, globalVar.Nz + 1, outputiStep):
#             f.write('{:10f} {:10f} {:10f}\n'.format(n * globalVar.deltat, i * globalVar.deltaz, T[n, i]))


plt.plot(np.linspace(0, 1, globalVar.Nz + 1), T[0, :], 'r-',
         np.linspace(0, 1, globalVar.Nz + 1), T[int(globalVar.Nt/2), :], 'g-',
         np.linspace(0, 1, globalVar.Nz + 1), T[-1, :], 'b-',
         )
plt.show()