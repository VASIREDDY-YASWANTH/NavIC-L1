import numpy as np
import math as m
from binary_fractions import Binary

c_t = [m.cos((i/255)*2*m.pi) for i in range(256)]
c_tbin = [Binary(i) for i in c_t]

s_t = [m.sin((i/255)*2*m.pi) for i in range(256)]
s_tbin=[Binary(i) for i in s_t]

with open('cos_lookup.dat','w') as f:
    for i in c_tbin :
        f.write(str(i)+'\n')

with open('sine_lookup.dat','w') as f:
    for i in s_tbin:
        f.write(str(i)+'\n')

print(m.sin(m.pi))

