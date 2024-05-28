import numpy as np
#a = np.loadtxt("dummy.dat",float,"\n")
a = np.loadtxt("pilotXsubc.dat",float,"\n")
b=a[:16]
#print(b)
#c=[-1, -1, -1, -1,   -1, 1, 1, 1,   1, 1, 1, 1,   1, 1, 1, -1] #2^4 16 point fft
fftr=np.fft.fft(b)
#print(fftr[:20])
for i in range(16):
    bi[i]=bin(b[i])
print(bi)
