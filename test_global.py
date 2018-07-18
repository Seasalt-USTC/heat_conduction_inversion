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

qm = 30  # mW * m^-2 == kW * km^-2
k = 3.35  # W * m^-1 * K^-1 == kW * km^-1 * k^-1
p = qm / k
u_mag = 0.3  # km * Ma^-1
rho_H0 = 2.5  # rho*H0 muW * m^-3 == kW * km^-3
sh0 = rho_H0 * globalVar.kappa / k  # K * Ma^-1
hr = 10  # km

T = CN_N(Ts=0, p=p, kappa=globalVar.kappa, u=globalVar.u_uniform * (-u_mag),
         Tic=globalVar.Tic_continent(p=p, sh0=sh0, hr=hr),
         sh=globalVar.sh_continental(u=globalVar.u_uniform * (-u_mag), sh0=sh0, hr=hr))
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
plt.text(5, 700, 'zTotal = {:<7}\n'
                 'qm = {:<7}\nk = {:<7}\nkappa = {:<7}\n'
                 'u = {:<7}\nrho*H0 = {:<7}\nhr = {:<7}'.
         format(globalVar.zTotal,
                qm, k, globalVar.kappa,
                u_mag, rho_H0, hr))
plt.plot(np.linspace(0, globalVar.zTotal, globalVar.Nz + 1), T[0, :], 'r-',
         np.linspace(0, globalVar.zTotal, globalVar.Nz + 1), T[int(globalVar.Nt/2), :], 'g-',
         np.linspace(0, globalVar.zTotal, globalVar.Nz + 1), T[-1, :], 'b-',
         )
plt.savefig(PATH + '/a.png')
plt.cla()

plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nt + 1), qs, 'r-', label='qs')
plt.legend()
plt.text(5, qs[0], 'zTotal = {:<7}\n'
                 'qm = {:<7}\nk = {:<7}\nkappa = {:<7}\n'
                 'u = {:<7}\nrho*H0 = {:<7}\nhr = {:<7}'.
         format(globalVar.zTotal,
                qm, k, globalVar.kappa,
                u_mag, rho_H0, hr))
plt.savefig(PATH + '/b.png')
plt.cla()

plt.xlim(0, 10)
plt.ylim(0, 100)
plt.plot(np.linspace(0, globalVar.zTotal, globalVar.Nz + 1), T[0, :], 'r-',
         np.linspace(0, globalVar.zTotal, globalVar.Nz + 1), T[int(globalVar.Nt/2), :], 'g-',
         np.linspace(0, globalVar.zTotal, globalVar.Nz + 1), T[-1, :], 'b-',
         )
plt.text(0.5, 70, 'zTotal = {:<7}\n'
                 'qm = {:<7}\nk = {:<7}\nkappa = {:<7}\n'
                 'u = {:<7}\nrho*H0 = {:<7}\nhr = {:<7}'.
         format(globalVar.zTotal,
                qm, k, globalVar.kappa,
                u_mag, rho_H0, hr))
plt.savefig(PATH + '/c.png')
plt.cla()

