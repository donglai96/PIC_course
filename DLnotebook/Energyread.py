import numpy as np 
import pandas as pd
import linecache
import matplotlib.pyplot as plt
# Read Energy file 

run_id = 21
t_num = 4000
file_name = '../mbeps1.source/output1.' + str(run_id)

def get_line_context(file_path, line_number):
    return linecache.getline(file_path, line_number).strip()

def make_energy_plot(file_name,energy_data, logscale = True):
    
    t = np.arange(energy_data.shape[0])
    plt.plot(t,energy_data[:,0],label = 'Total Field')
    plt.plot(t,energy_data[:,1],label = 'Kinetic Energy')
    plt.plot(t,energy_data[:,2],label = 'Total Energy')
    plt.legend()
    if logscale:
        plt.yscale('log')
    plt.title('Energy diagnostic')
    plt.savefig(file_name,format = 'jpg')
    



# test_string = get_line_context(file_name,4)

# df = list(map(float, test_string.split(' ')))
# print(get_line_context(file_name,4))
# print(df)
if __name__ == "__main__":
    total_Energy = np.zeros((t_num,3)) # 3 is the Total Field, Kinetic and Total energies
    for i in range(t_num):
        line_num = 5 * i +4
        test_string = get_line_context(file_name,line_num)
        list_line = list(map(float, test_string.split(' ')))
        np_line = np.array(list_line)
        total_Energy[i,:] = np_line 
    print(total_Energy)

    make_energy_plot('Energy' + str(run_id) +'.jpg',energy_data = total_Energy)

    

