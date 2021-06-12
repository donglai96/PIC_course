import numpy as np 
import matplotlib.pyplot as plt
import math


import sys
import math
import numpy


# sys.path.append('../ReadFiles')
# from fgraf import *
# from libcmfield1 import *


run_id = 26
t_num = 4000
#mode_num = 41

INDX = 9
NTB = 5
MODESXB = 256
FBNAME = 'bk1'
nbrec = t_num // NTB
norm = 0
# Not very sure this dim , right now just think it's
mdim = 2
complex_type = numpy.complex64


#allocate complex scalar array
nx = int(math.pow(2,INDX))
nxh = int(nx/2)


vfieldc = numpy.empty((mdim,MODESXB),complex_type,'F')


file_name = '../mbeps1.source/bk1.' + str(run_id)
bk1 = np.fromfile(file_name,dtype = numpy.complex64)
print(bk1.shape)
bk1_reshape = bk1.reshape(nbrec ,MODESXB,mdim)
print(type(bk1_reshape[0,0,0]))
print(bk1_reshape.shape)

# Change the BK to real place


print(bk1_reshape[0,5,1])



# plt.plot(range(MODESXB),np.absolute(bk1_reshape[200,:,0]),label = 'y')
# plt.plot(range(MODESXB),np.absolute(bk1_reshape[200,:,1]),label = 'z')
# plt.legend()
# plt.title('BFIELD MODES')
# plt.savefig('BK_test2.jpg',format = 'jpg')
# plt.show()
# # change this to real space



def decompressVCfield(vfc,nx):
    """This function decompress the field from VC field

    Args:
        vf ([type]): vfield after decompressed
        vfc ([type]):  input vfield
                        unlike fortran process vfc here is (time,k,2)
                        the vf is (time,nx/2,2)
        nx ([type]): length of box
        modex ([type]): modenumber of your k
    return: vf
    """
    print("The input shape is ",vfc.shape)
    if nx < vfc.shape[1]:
        raise ValueError('the box number is less than k number')
    
    nxh = nx//2 
    vf = np.empty((vfc.shape[0],nxh,mdim),complex_type,'F')
    modex = vfc.shape[1]
    print(vf.shape)
    jmax = min(nxh,modex)
    j1 = nxh + 1
    # Mode numbers 0 < kx < nx/2
    print(vfc.shape)
    print(vfc.real)
    for j in np.arange(1,jmax,1):
        for i in range(mdim):
            vf[:,j,i] = vfc[:,j,i]
    
    for j in np.arange(jmax,nxh,1):
        for i in range(mdim):
            vf[:,j,i] = np.complex(0.0)
    #kx = 0
    for i in range(mdim):
        vf[:,0,i] = vfc.real[:,0,i] + 0j
        #print(np.complex(vfc.real[:,0,i]))
        if modex > nxh:
            vf[:,0,i] = vf.real[:,0,i]+ 1j*vfc.real[:,nxh,i]

    return vf

bf = decompressVCfield(bk1_reshape,nx)
print(bf.shape)

#Change the bfc back to real space 
by= np.fft.irfft(bk1_reshape[0,:,0].squeeze(),n = 512)
bz= np.fft.irfft(bk1_reshape[0,:,1].squeeze(),n = 512)

plt.plot(range(512),512*(bz),label = 'z')
plt.plot(range(512),512*(by),label = 'y')
plt.legend()
plt.title('BFIELD MODES')
plt.savefig('Bx_test26.jpg',format = 'jpg')


#Now calculate the bywk and bzwk

# by_all = np.zeros((nbrec,int(2**INDX)))
# bz_all = np.zeros((nbrec,int(2**INDX)))
# for i in range(nbrec):
#     by_all[i,:] = np.fft.irfft(bk1_reshape[i,:,0].squeeze(),n = 512)
#     bz_all[i,:] = np.fft.irfft(bk1_reshape[i,:,0].squeeze(),n = 512)

# plt.imshow(bz_all,cmap = 'rainbow')
# plt.colorbar()
# plt.savefig('test_spec')











# def bwk_plot(figname, bk,tstart,tend):
#     bk = bk[tstart:tend,:]
#     N = bk.shape[0]
#     print(N)
#     freq = np.fft.fftfreq(bk.shape[0] )
#     bkw = bk.shape[0] * np.fft.ifft(bk,axis = 0)
#     bkw_shift = np.roll(bkw,-N//2,axis = 0)
    
#     Z = np.log10(np.absolute(bkw_shift)**2)[N//4:N*3//4,:]
#     x = np.arange(-np.pi/2,np.pi/2,np.pi/(N//2)) 
#     y = np.arange(0, bk.shape[1], 1)  # len = 7
#     X, Y = np.meshgrid(x, y)
#     fig, ax = plt.subplots(figsize=(16,8))
#     dis = ax.pcolormesh(X, Y, Z.T,vmin = -7,vmax = -1,cmap = 'rainbow')
#     ax.set_xlabel(r'$\omega / \omega_{pe}$')
#     ax.set_ylabel('mode number')
#     cbar = fig.colorbar(dis, ax=ax)
#     cbar.ax.set_ylabel(r'$log_{10}(|\delta E_L(\omega,k_x)|)$')
#     #wplot = np.arange(0,200/(7/0.0194),0.01)
#     #disper = ax.plot(wplot,wplot*7/0.0194,'k')
#     plt.savefig(figname + '.jpg')
#     plt.title(figname)
#     plt.show()





# by_k = bk1_reshape[:,:,0]
# bz_k = bk1_reshape[:,:,1]
# #bwk_plot('B_z',bz_k,tstart = 0,tend = 800)