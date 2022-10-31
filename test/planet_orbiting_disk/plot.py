import numpy as np
import matplotlib.pyplot as plt


p_LF = np.loadtxt("LF_disk_at_center/monitor/bigplanet2.dat", skiprows=22)
p_e = np.loadtxt("euler_disk_at_center/monitor/bigplanet2.dat", skiprows=22)
#p_e = np.loadtxt("LF_disk_at_center/monitor/bigplanet2.dat", skiprows=22)
p_ref = np.loadtxt("star_at_center/monitor/bigplanet2.dat", skiprows=22)

t_ref = p_ref[:,6]
x_ref = p_ref[:,1]
y_ref = p_ref[:,2]

vx_ref = p_ref[:,3]
vy_ref = p_ref[:,4]



t_e = p_e[:,6]
x_e = p_e[:,1]
y_e = p_e[:,2]
vx_e = p_e[:,3]
vy_e = p_e[:,4]

t_LF = p_LF[:,6]
x_LF = p_LF[:,1]
y_LF = p_LF[:,2]
vx_LF = p_LF[:,3]
vy_LF = p_LF[:,4]


fig, axs = plt.subplots(2)

axs[0].set_title("Position")
axs[0].plot(t_ref,x_ref, '-k', lw=4, label='ref')
axs[0].plot(t_ref,y_ref, '--k', lw=4)
axs[0].plot(t_e,x_e, '-r', lw=3, label='euler')
axs[0].plot(t_e,y_e, '--r', lw=3)
axs[0].plot(t_LF,x_LF, '-b', label='LF')
axs[0].plot(t_LF,y_LF, '--b')
axs[0].legend(loc='upper right')

axs[1].set_title("Velocity")
axs[1].plot(t_ref,vx_ref, '-k', lw=4, label='ref')
axs[1].plot(t_ref,vy_ref, '--k', lw=4)
axs[1].plot(t_e,vx_e, '-r', lw=3, label='euler')
axs[1].plot(t_e,vy_e, '--r', lw=3)
axs[1].plot(t_LF,vx_LF, '-b', label='LF')
axs[1].plot(t_LF,vy_LF, '--b')
axs[1].legend(loc='upper right')


plt.show()
