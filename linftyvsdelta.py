import matplotlib.pyplot as plt
import numpy as np






x = np.array([0.01,0.005,0.0025, 0.001, 0.0005])
y = np.array([511.8600618733877,1044.1823505462148,2108.7395648442466,5254.817360315813,10293.840421881721])
plt.ion()
# plt.loglog(x,x,label='x')
# plt.plot(x,y, label='diadata')
plt.loglog(x,y, label='Log-Plot', marker= '.')
m,c = np.polyfit(np.log(x), np.log(y), 1)
yfit = np.exp(m*np.log(x) + c+10**(-1))
plt.plot(x,yfit)
print("c*x^(m) with m={},c={}".format(m,c))
plt.legend()
plt.show()
