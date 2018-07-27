from invent_data import *

Td = np.loadtxt('invent_data/from_Tic_real_and_u/Td.txt')
plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nt + 1), Td, label='from Tic_real and u')
Td = np.loadtxt('invent_data/from_Tp_and_u/Td.txt')
plt.plot(np.linspace(0, globalVar.tTotal, globalVar.Nt + 1), Td, label='from Tp and u')
plt.legend()
plt.title('Td')
plt.savefig('invent_data/contrst.png')
plt.cla()

# T = np.loadtxt('invent_data/from_Tic_real_and_u/T.txt')
# Td = np.loadtxt('invent_data/from_Tic_real_and_u/Td.txt')
# u = get_velocity_transverse(T=T, Td=Td)
# plt.plot(u)
# plt.show()
# plt.cla()
#
# T = np.loadtxt('invent_data/from_Tp_and_u/T.txt')
# Td = np.loadtxt('invent_data/from_Tp_and_u/Td.txt')
# u = get_velocity_transverse(T=T, Td=Td)
# plt.plot(u)
# plt.show()
# plt.cla()