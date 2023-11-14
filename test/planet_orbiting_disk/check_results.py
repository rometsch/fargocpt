import numpy as np
import matplotlib.pyplot as plt

import yaml

def test(_):

    with open("testconfig.yml", 'r') as ymlfile:
        testconfig = yaml.safe_load(ymlfile)
    testname = testconfig['testname']

    p_LF = np.loadtxt(f"../../output/tests/{testname}/LF_disk_at_center/monitor/nbody1.dat", skiprows=22)
    p_e = np.loadtxt(f"../../output/tests/{testname}/euler_disk_at_center/monitor/nbody1.dat", skiprows=22)
    #p_e = np.loadtxt("LF_disk_at_center/monitor/nbody1.dat", skiprows=22)
    p_ref = np.loadtxt(f"../../output/tests/{testname}/star_at_center/monitor/nbody1.dat", skiprows=22)

    t_ref = p_ref[:,7]
    x_ref = p_ref[:,2]
    y_ref = p_ref[:,3]
    vx_ref = p_ref[:,4]
    vy_ref = p_ref[:,5]

    t_e = p_e[:,7]
    x_e = p_e[:,2]
    y_e = p_e[:,3]
    vx_e = p_e[:,4]
    vy_e = p_e[:,5]

    t_LF = p_LF[:,7]
    x_LF = p_LF[:,2]
    y_LF = p_LF[:,3]
    vx_LF = p_LF[:,4]
    vy_LF = p_LF[:,5]


    fig, axs = plt.subplots(2,2)
    axs = axs.flatten()

    lw_back = 4
    lw_mid = 2

    axs[0].plot(t_ref,x_ref, '-k', label='ref', lw=lw_back)
    axs[0].plot(t_e,x_e, '-r', label='euler', lw=lw_mid)
    axs[0].plot(t_LF,x_LF, '-b', label='LF')
    axs[0].legend(loc='upper right')
    axs[0].set_ylabel("x")
    
    axs[1].plot(t_ref,y_ref, '--k', lw=lw_back)
    axs[1].plot(t_e,y_e, '--r', lw=lw_mid)
    axs[1].plot(t_LF,y_LF, '--b')
    axs[1].set_ylabel("y")

    axs[2].plot(t_ref,vx_ref, '-k', label='ref', lw=lw_back)
    axs[2].plot(t_e,vx_e, '-r', label='euler', lw=lw_mid)
    axs[2].plot(t_LF,vx_LF, '-b', label='LF')
    axs[2].legend(loc='upper right')
    axs[2].set_ylabel("vx")

    axs[3].plot(t_ref,vy_ref, '--k', lw=lw_back)
    axs[3].plot(t_e,vy_e, '--r', lw=lw_mid)
    axs[3].plot(t_LF,vy_LF, '--b')
    axs[3].set_ylabel("vy")

    fig.savefig("plot.jpg", dpi=150)

    diff_e_x = np.max(np.abs(x_e - x_ref))
    diff_e_y = np.max(np.abs(y_e - y_ref))

    diff_LF_x = np.max(np.abs(x_LF - x_ref))
    diff_LF_y = np.max(np.abs(y_LF - y_ref))

    threshold_LF = float(testconfig['threshold_LF'])
    threshold_euler = float(testconfig['threshold_Euler'])
    
    with open("differences.txt", "w") as f:
        print(f"Euler position difference: |x - x_ref| = {diff_e_x}, |y - y_ref| = {diff_e_y}| (threshold = {threshold_euler})", file=f)
        print(f"LF position difference: |x - x_ref| = {diff_LF_x}, |y - y_ref| = {diff_LF_y}| (threshold = {threshold_LF})", file=f)

    euler_pass = diff_e_x < threshold_euler and diff_e_y < threshold_euler
    LF_pass = diff_LF_x < threshold_LF and diff_LF_y < threshold_LF

    with open("test.log", "w") as f:
        from datetime import datetime
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"{current_time}", file=f)
        print(f"Euler position difference: |x - x_ref| = {diff_e_x}, |y - y_ref| = {diff_e_y}| (threshold = {threshold_euler})", file=f)
        print(f"LF position difference: |x - x_ref| = {diff_LF_x}, |y - y_ref| = {diff_LF_y}| (threshold = {threshold_LF})", file=f)



    if euler_pass and LF_pass:
        print(f"SUCCESS: {testname}")
    else:
        print(f"FAIL: {testname}")


if __name__ == "__main__":
    test("foo")
