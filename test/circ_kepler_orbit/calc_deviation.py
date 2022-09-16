#!/usr/bin/env python3
import numpy as np
import os

def main():
    t, x, y = np.genfromtxt("output/monitor/bigplanet2.dat", usecols=(6,1,2), unpack=True)
    t = t

    OmegaP = np.sqrt(1 + 1e-3)
    x_theo = np.cos(t*OmegaP)
    y_theo = np.sin(t*OmegaP)


    delta_x = x - x_theo
    diff_x = np.max(np.abs(delta_x))

    delta_y = y - y_theo
    diff_y = np.max(np.abs(delta_y))

    test_name = "Circular Kepler Problem"

    threshold = 1e-11
    if max(diff_x, diff_y) > threshold:
        print(f"FAIL: {test_name} at {os.getcwd()}")
        print(f"Diff x = {diff_x}, Diff y = {diff_y}, threshold = {threshold}")
    else:
        print(f"SUCCESS: {test_name}")



if __name__=="__main__":
    main()