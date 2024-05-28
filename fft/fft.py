import numpy as np
a = np.loadtxt("dummy.dat",float,"\n")


#a=[1, 7, 2, 4,   8, 6, 4, 9,   3, 6, 8, 5,   9, 3, 4, 9,     3, 2, 8, 2,   6, 4, 5, 8,  4, 5, 2, 8,  8, 2, 6, 9]  # 2^5 35 point fft
#a=[1, 7, 2, 4,   8, 6, 4, 9,   3, 6, 8, 5,   9, 3, 4, 9]  # 2^4 16 point fft
#a=[1, 7, 2, 4,   8, 6, 4, 9]  # 2^3  8 point fft
#a=[1, 7, 2, 4]  # 2^2  4 point fft

#a=[-1, -1, -1, -1,   -1, 1, 1, 1,   1, 1, 1, 1,   1, 1, 1, -1] #2^4 16 point fft
#a=[-1, -1, -1, -1,   -1, -1, 1, 1,   1, 1, 1, 1,   1, 1, 1, 1] #2^4 16 point fft
#a=[-1, -1, -1, -1,   -1, 1, 1, 1] #2^3  8 point fft
#a=[-1, -1, -1, -1,   -1, -1, 1, 1] #2^3  8 point fft
#a=[-1, -1, -1, -1] #2^2  4 point fft
b=a[:2**5]

fftr=np.fft.fft(b)
print(fftr)
with open('pyfftout.dat','w') as f:
        for i in fftr:
            f.write(str(i)+'\n')
b.tofile('data2.csv', sep = ',')
