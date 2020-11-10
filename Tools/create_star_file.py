#!/usr/bin/env python3
import numpy as np
import pandas as pd
import sys
import os
from shutil import copy2

src = sys.argv[1]
dst = list(src)
dst[-5] = '2'
dst = "".join(dst)

print(src, dst)

copy2(src,dst)

planet = np.loadtxt(src)
ids = planet[:,0]
time = planet[:,6]
mass = np.ones_like(ids)

star = np.zeros_like(planet)
star[:,0] = ids.astype(int)
star[:,5] = mass
star[:,6] = time

df = pd.DataFrame(data=star)
df = df.set_index(0)
df.index = df.index.values.astype(int)

with open(src) as file:
    lines = [next(file) for x in range(16)]

header = "".join(lines)

df.to_csv(src, sep='\t', index=True, header=False)
with open(src, "r+") as f:
    old = f.read() # read everything in the file
    f.seek(0) # rewind
    f.write(header) # write the new line before
    f.write(old) # write the old content

print("Success\n")
