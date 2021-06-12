import numpy as np
import matplotlib.pyplot as plt


file_name = '../fve1.26'

data = np.fromfile(file_name,dtype = np.float32)
data = data.reshape(800,int(204800//800))
print(data.shape)
print(data[:,3:9])
plt.plot(data[0,:])
#print(data[:,-10])
plt.ylim(0,50)
plt.savefig('testv.jpg')