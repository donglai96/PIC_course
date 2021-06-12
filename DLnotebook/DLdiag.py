"""
    Diagnostic file for Vicktor's UPIC Code
    right now only support 1d 
    Donglai Ma, UCLA 2021
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker 
import matplotlib.cm as cm
import matplotlib.colors as colors
import os 
import math
import linecache
from mpl_toolkits.axes_grid1 import make_axes_locatable

class DL_diag(object):
    def __init__(self,IDRUN,mode,folder): 
        self.IDRUN = IDRUN  
        self.mode = mode # ES,EM DW
        self.folder =folder
        print("doing diagnostic...")
        print("The ID run is :",IDRUN) 
        print("The code mode is :", mode)



        # read the diag file, for simple just use hand for now
        ############################### 
        # When it's capital means the parameter could be found in the file
        self.NTB = 5 
        self.TEND = 100 
        self.DT = 0.025
        self.tstep = 4000 
        self.VTX = 1 
        self.VTY = 2 
        self.VTZ = 1 
        self.MODESXB = 256
        self.INDX = 9
        self.mdim =2
        self.nx = int(math.pow(2,self.INDX))
        self.energyplot = [0,5] # 0 total field, 5 is magnetic field



        
        print("support list:" )
        if self.mode == 'EM':
            self.supportlist = ['B_k_t','B_x_t','B_w_k','Energy_t','V_dis']
        # ES DW left here
        print(self.supportlist)
        self.diagfile = self.folder + '/diag1.' + str(IDRUN)
        print('The diagfile is at:',self.diagfile)
    
    def TimeRangeStudy(self,timerange = None,target_list = ['B_k_t','B_x_t','Energy_t','V_dis'],plot = True):
        """
        This function give me the full time analysis, store the total data into self, 
        should warning about the size
        """
        ################################################################
        # right now the time step lets say is the total time steps
        # The multiple plot logic is poor

        self.timerange = timerange # An array with start time and end time
        # Check the number of target_list right now default is 4 figures :
        fig,axs = plt.subplots(3,2,figsize=(16,20),sharex = True)
        #for target in target_list[1]:
        ###############################
        for target,i in zip(target_list[0:4],range(len(target_list[0:4]))):
            print('The %d plot is : '%i, target)
            if target =='B_k_t':
                
                self.bkt = self.getbk()
                print('The bkt shape is:',self.bkt.shape)
                print('The first is time ,the second is mode, the third is dimension')
                print('The bkt type is:',self.bkt.dtype)

                # Make the plot
                # test plot
        #         print(i)
        #         axs[i].plot(range(self.MODESXB),np.absolute(self.bkt[0,:,0]),label = 'y')
        #         axs[i].plot(range(self.MODESXB),np.absolute(self.bkt[0,:,1]),label = 'z')
        #         axs[i].legend()
        #         axs[i].set_xlabel('mode')
        #         axs[i].set_title('Bk field')
        # plt.savefig("testdiag.jpg")
                bykt_abs = np.absolute(self.bkt[:,:,0])
                bzkt_abs = np.absolute(self.bkt[:,:,1])
                Z1 = np.flipud(bykt_abs.T)
                Z2 = np.flipud(bzkt_abs.T)
                im_bykt = axs[0,0].imshow((Z1), cmap=cm.rainbow, extent=[0, int(self.nbrec*self.NTB*self.DT), 1,self.MODESXB],aspect='auto')
                im_bzkt = axs[0,1].imshow((Z2), cmap=cm.rainbow, extent=[0, int(self.nbrec*self.NTB*self.DT), 1,self.MODESXB],aspect='auto')
                #divider00 = make_axes_locatable(axs[0,0])
                #cax00 = divider00.append_axes("right", size="5%", pad=0.05)
                cbar00 = fig.colorbar(im_bykt, shrink=0.9,ax = axs[0,0],orientation = 'horizontal')

                #divider01 = make_axes_locatable(axs[0,1])
                #cax01 = divider01.append_axes("right", size="5%", pad=0.05)
                cbar01 = fig.colorbar(im_bzkt, shrink=0.9,ax = axs[0,1],orientation = 'horizontal')
                axs[0,0].set_yscale('log')
                axs[0,0].set_title('By k-t ')
                axs[0,1].set_yscale('log')
                axs[0,1].set_title('Bz k-t ')
                axs[0,0].set_xlim([0,self.tstep*self.DT])
            
            if target == 'B_x_t':
                self.bkt = self.getbk()

                #change bkt to bxt #c means complex
                self.b_xt_c = self.decompressVCfield(self.bkt,self.nx)
                self.by_xt= self.nx * np.fft.irfft(self.b_xt_c[:,:,0],axis = 1,n = self.nx)
                #test 
                print(self.by_xt.shape)
                self.bz_xt= self.nx * np.fft.irfft(self.b_xt_c[:,:,1],axis = 1,n = self.nx)
                # axs[1,1].plot(range(self.nx),self.by_xt[0,:],label = 'y')
                # axs[1,1].plot(range(self.nx),self.bz_xt[0,:],label = 'z')
                # axs[1,1].legend()
                ZZ1 = np.flipud(self.by_xt.T)
                ZZ2 = np.flipud(self.bz_xt.T)
                im_bykt = axs[1,0].imshow((ZZ1), cmap=cm.bwr, extent=[0, int(self.nbrec*self.NTB*self.DT), 1,self.nx],aspect='auto')
                im_bzkt = axs[1,1].imshow((ZZ2), cmap=cm.bwr, extent=[0, int(self.nbrec*self.NTB*self.DT), 1,self.nx],aspect='auto')
                cbar0 = fig.colorbar(im_bykt, shrink=0.9,ax = axs[1,0],orientation = 'horizontal')
                cbar1 = fig.colorbar(im_bzkt, shrink=0.9,ax = axs[1,1],orientation = 'horizontal')
                #axs[0,0].set_yscale('log')
                axs[1,0].set_title('By x-t ')
                #axs[0,1].set_yscale('log')
                axs[1,1].set_title('Bz x-t ')


            
            if target == 'Energy_t':
                self.energymatrix = self.getenergy()
                print(self.energymatrix.shape)
                axs[2,1].plot(np.arange(0,int(self.tstep*self.DT),self.DT),self.energymatrix[:,0],label = 'Total Energy')
                axs[2,1].plot(np.arange(0,int(self.tstep*self.DT),self.DT),self.energymatrix[:,5],label = 'Magnetic')
                axs[2,1].set_ylabel('Energy')
                axs[2,1].set_xlabel('Time')
                axs[2,1].set_title('ENERGY')
                axs[2,1].legend()
                # How to get magnetic field energy

            if target == 'V_dis':
                self.fth = self.getfve()
                axs[2,0].plot(np.arange(0,int(self.tstep*self.DT),self.DT*self.NTB),self.fth[:,0],label = 'x')
                axs[2,0].plot(np.arange(0,int(self.tstep*self.DT),self.DT*self.NTB),self.fth[:,1],label = 'y')
                axs[2,0].plot(np.arange(0,int(self.tstep*self.DT),self.DT*self.NTB),self.fth[:,2],label = 'z')
                axs[2,0].set_ylabel('V')
                axs[2,0].set_xlabel('Time')
                axs[2,0].set_title('Thermal Velocity')   
                axs[2,0].legend()



        plt.savefig("testdiag_" + str(self.IDRUN) + '.jpg')

    def TimeSliceStudy(self,timeslices):
        """This function gives me the time slice study 

        Args:
            timeslices (Array): Array of interested time slices
        """
    # @staticmethod
    # def lineplot()

    def getbk(self):
        """
        This function gives the bk1 data in the bk1 file

        return: nparray of bk1
        """

        # check the file folder
        bk1file = self.folder + '/' + 'bk1.' + str(self.IDRUN)
        if os.path.isfile(bk1file) is False:
            raise ValueError("The bk1 file is not exist")
        else:# read the file
            complex_type = np.complex64
            #nx = int(math.pow(2,self.INDX))
            bk1 = np.fromfile(bk1file,dtype = complex_type)
            self.nbrec =int( self.tstep // self.NTB)
            bk1_reshape = bk1.reshape(self.nbrec ,self.MODESXB,self.mdim)

        return bk1_reshape
    
    def getenergy(self):
        """This function gives the energy data stored in output file

        Notice that fortran version has no E field and B field (or maybe have)
        so this requries the python3 mbbeps1.py to get the new version outputfile

        Returns:
            energy matrix
        """
        energyfile = self.folder + '/' + 'output1.' + str(self.IDRUN)
        if os.path.isfile(energyfile) is False:
            raise ValueError("The output1 file is not exist")
        else:
            total_Energy = np.zeros((self.tstep,6))
            for i in range(self.tstep):#notice here might not be time step, need to check 
                line_num_1 = 5 * i + 4
                line_num_2 = 5 * i + 6 
                test_string_1 = self.get_line_context(energyfile,line_num_1)
                test_string_2 = self.get_line_context(energyfile,line_num_2)
                #print(test_string_1.split())
                #print(list(map(float, test_string_1.split(' '))))
                list_line_1 = np.array([float(n) for n in test_string_1.split()])
                list_line_2 = np.array([float(n) for n in test_string_2.split()])
                total_Energy[i,:3] = list_line_1
                total_Energy[i,3:] = list_line_2

        return total_Energy
    
    def getfve(self):
        """This function gives the fve of electrons, only read for thermal velocity for now

        Raises:
            ValueError: [description]

        Returns:
            [type]: [description]
        """
        fvefile = self.folder + '/' + 'fve1.' + str(self.IDRUN)
        if os.path.isfile(fvefile) is False:
            raise ValueError("The fve1 file is not exist")
        else:
            fve = np.fromfile(fvefile,dtype = np.float32)

            fve = fve.reshape(self.nbrec,fve.shape[0]//self.nbrec)
            print(fve.shape)
            #only check the thermal
            return fve[:,3:6]

    def calculateGrowth(self,mode = 5):
        bzkt_abs = np.absolute(self.bkt[:,:,1])
        fig,axs = plt.subplots(2,figsize = (16,16))
        axs[0].plot(np.arange(0,int(self.tstep*self.DT),self.DT*self.NTB),bzkt_abs[:,mode])
        axs[0].set_title('single mode B')
        #calculate the growth rate
        gamma_t  = np.log(bzkt_abs[1:,mode]/bzkt_abs[:-1,mode])
        axs[1].plot(np.arange(0,int(self.tstep*self.DT) -self.DT*self.NTB ,self.DT*self.NTB),gamma_t *8,label = 'growth rate')
        axs[1].set_title('Growth_rate')
        #calculate the trapping frequency
        v = self.fth[:,1].squeeze()
        B = self.bz_xt.max(axis=1)
        wtrap = np.sqrt((2 *np.pi *mode/self.nx) * v * B)
        axs[1].plot(np.arange(0,int(self.tstep*self.DT)  ,self.DT*self.NTB),wtrap,label = 'trap frequency')
        #trap_f = np.sqrt()
        axs[1].legend()
        plt.savefig('testgrowth_'+ str(self.IDRUN) + '.jpg')

    def extraplot(self):
        fig,ax = plt.subplots(2,figsize=(16,16))
        a = np.arange(1,101,1)
        norm = np.sqrt(8/(27*np.pi))*0.1
        ax[0].plot(a,norm *(a - 1)**(1.5)/a,color = 'red',label = 'theory')
        ax[0].scatter(np.array([4,9,16,25,36,49,64,81,100]),np.array([0.05,0.1,0.13,0.15,0.18,0.21,0.22,0.26,0.28]))
        plt.savefig('testtheory.jpg')
    @staticmethod
    def decompressVCfield(vfc,nx,mdim = 2):
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
        vf = np.empty((vfc.shape[0],nxh,mdim),vfc.dtype,'F')
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
    
    @staticmethod
    def get_line_context(file_path, line_number):
        return linecache.getline(file_path, line_number).strip()


    



if __name__ == '__main__':
    cl = DL_diag(16,mode = 'EM',folder = '../')
    cl.TimeRangeStudy()
    cl.calculateGrowth(mode = 5)
    cl.extraplot()

        
