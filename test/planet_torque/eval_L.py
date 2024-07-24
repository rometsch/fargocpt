# %%
import numpy as np
import matplotlib.pyplot as plt
%matplotlib widget

path = "../../output/tests/planet_torque/out/monitor/nbody1.dat"

data = np.loadtxt(path)

time = data[:, 1]
mask = np.unique(time, return_index=True)[1]

data_new = np.zeros((len(mask), len(data[0,:])))
for i in range(len(data[0,:])):
    data_new[:, i] = data[:, i][mask]

data = data_new

time = np.copy(data[:, 7])
dt = data[:, 7]
dt = dt[1]-dt[0]

L = data[:, 11]

torque = data[:, 18]
acc_torque = data[:, 19]
ind_torque = data[:, 20]


plt.figure()
plt.plot(time, L, label='Angular Momentum L')
#plt.plot(time, torque)
#plt.plot(time, acc_torque)
#plt.plot(time, ind_torque)
plt.plot(time, L - np.cumsum(torque*dt) - np.cumsum(acc_torque*dt) - np.cumsum(ind_torque*dt), label='L + torques integrated')
plt.legend()
plt.show()