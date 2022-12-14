#!/usr/bin/env python3
import numpy as np

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('file1')
    parser.add_argument('file2')
    args = parser.parse_args()

    file1 = args.file1
    file2 = args.file2

    compare(file1, file2)

def compare(file1, file2):
    data1 = np.fromfile(file1)
    data2 = np.fromfile(file2)

    delta = np.abs(data1-data2)
    reldelta = delta/np.abs(data1)
    Ndiff = np.sum(delta!=0)
    Ntot = len(data1)

    delta = delta[ delta!=0 ]

    print("Number of doubles = ", Ntot)
    print("Number of different values = ", Ndiff)
    print("min value (file1) = ", np.min(data1))
    print("max value (file1) = ", np.max(data1))
    print("min value (file2) = ", np.min(data2))
    print("max value (file2) = ", np.max(data2))
    if (Ndiff != 0):
        print("min diff = ", np.min(delta))
        print("max diff = ", np.max(delta))
        print("avg diff = ", np.mean(delta))
        print("median diff = ", np.median(delta))
        print("min relative diff = ", np.min(reldelta))
        print("max relative diff = ", np.max(reldelta))
        print("avg relative diff = ", np.mean(reldelta))
        print("median relative diff = ", np.median(reldelta))
        print("max relative diff id = ", np.argmax(reldelta))
    else:
        print("Files are identical")

if __name__=="__main__":
    main()
