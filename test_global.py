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

# T = CN_N(Ts=0, p=globalVar.pm, kappa=globalVar.kappa, u=globalVar.u,
#          Tic=globalVar.Tic_continent(p=globalVar.pm, sh0=globalVar.sh0, hr=globalVar.hr),
#          sh=globalVar.sh_continental(u=globalVar.u, sh0=globalVar.sh0, hr=globalVar.hr))
T = CN_N(Ts=0, p=globalVar.pm, kappa=globalVar.kappa, u=globalVar.u,
         Tic=globalVar.Tic_continent(p=globalVar.pm, sh0=globalVar.sh0, hr=globalVar.hr),
         sh=globalVar.sh_continental(u=globalVar.u, sh0=globalVar.sh0, hr=globalVar.hr))

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

n = 1
PATH = 'Continental_Geotherms/' + '{:0>3}'.format(n)
while os.path.exists(PATH):
    n += 1
    PATH = 'Continental_Geotherms/' + '{:0>3}'.format(n)
os.mkdir(PATH)
qs = (T[:, 1] - T[:, 0]) / globalVar.deltaz * 3.35

plt.cla()
plt.xlim(0, 100)
plt.ylim(0, 1000)
plt.text(0.05, 0.6, 'zTotal = {:<7}\n'
                 'qm = {:<7}\nk = {:<7}\nkappa = {:<7}\n'
                 'u = {:<7}\nrho*H0 = {:<7}\nhr = {:<7}'.
         format(globalVar.zTotal,
                globalVar.qm, globalVar.k, globalVar.kappa,
                globalVar.u_mag, globalVar.rho_H0, globalVar.hr),
         transform=plt.gca().transAxes)
plt.plot(np.linspace(0, globalVar.zTotal, globalVar.Nz + 1), T[0, :], 'r-',
         np.linspace(0, globalVar.zTotal, globalVar.Nz + 1), T[int(globalVar.Nt/2), :], 'g-',
         np.linspace(0, globalVar.zTotal, globalVar.Nz + 1), T[-1, :], 'b-',
         )
plt.savefig(PATH + '/a.png')
plt.cla()

plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nt + 1), qs, 'r-', label='qs')
plt.legend(loc='lower right')
plt.text(0.05, 0.6, 'zTotal = {:<7}\n'
                 'qm = {:<7}\nk = {:<7}\nkappa = {:<7}\n'
                 'u = {:<7}\nrho*H0 = {:<7}\nhr = {:<7}\n'
                 'Pe = {:<5.3}'.
         format(globalVar.zTotal,
                globalVar.qm, globalVar.k, globalVar.kappa,
                globalVar.u_mag, globalVar.rho_H0, globalVar.hr,
                globalVar.Pe),
         transform=plt.gca().transAxes)
plt.savefig(PATH + '/b.png')
plt.cla()

plt.xlim(0, 10)
plt.ylim(0, 100)
plt.plot(np.linspace(0, globalVar.zTotal, globalVar.Nz + 1), T[0, :], 'r-',
         np.linspace(0, globalVar.zTotal, globalVar.Nz + 1), T[int(globalVar.Nt/2), :], 'g-',
         np.linspace(0, globalVar.zTotal, globalVar.Nz + 1), T[-1, :], 'b-',
         )
plt.text(0.05, 0.6, 'zTotal = {:<7}\n'
                 'qm = {:<7}\nk = {:<7}\nkappa = {:<7}\n'
                 'u = {:<7}\nrho*H0 = {:<7}\nhr = {:<7}'.
         format(globalVar.zTotal,
                globalVar.qm, globalVar.k, globalVar.kappa,
                globalVar.u_mag, globalVar.rho_H0, globalVar.hr),
         transform=plt.gca().transAxes)
plt.savefig(PATH + '/c.png')
plt.cla()

