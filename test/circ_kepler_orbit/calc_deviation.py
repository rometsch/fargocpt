#!/usr/bin/env python3
import numpy as np

def main():
    t, x, y = np.genfromtxt("output/monitor/bigplanet2.dat", usecols=(6,1,2), unpack=True)
    t = t

    OmegaP = np.sqrt(1 + 1e-3)
    x_theo = np.cos(t*OmegaP)
    y_theo = np.sin(t*OmegaP)


    delta_x = x - x_theo
    diff_x = np.max(np.abs(delta_x))
    print(f"Diff x = {diff_x}")

    delta_y = y - y_theo
    diff_y = np.max(np.abs(delta_y))
    print(f"Diff y = {diff_y}")

if __name__=="__main__":
    main()