import numpy as np
import matplotlib.pyplot as plt

def PlotSlumping():
    print('hi you there, im plotting')
    data2 = np.loadtxt('time_data.dat')
    #data3 = np.loadtxt('geometric_data.dat')
    #data = np.loadtxt('density_data.dat')
    #print(data)
    #print(data3)
    #for i in range(0,len(data),):
    multiple_data = ['10','100','500','1000']
    
    for i in multiple_data:
        data3 = np.loadtxt('geometric_data_'+i+'.dat')
        data = np.loadtxt('density_data_'+i+'.dat')
        plt.plot(data3,data[-1],label=i)
    
    plt.legend(title='Spatial Nodes')
    plt.title('Grid Independance at t=500s for Slumping process')
    plt.xlabel('Relative height up cylinder from bottom (x/a)')
    plt.ylabel(r'Relative density change ($\rho/\rho_{0}$)')
    #plt.show()
    plt.savefig('slumping_plot.png')

def PlotSlump():
    print('hi you there, im plotting')
    data2 = np.loadtxt('time_data.dat')
    data3 = np.loadtxt('geometric_data.dat')
    data = np.loadtxt('density_data.dat')
    print(len(data),len(data2),len(data3))
    for i in range(0,len(data),100):
        plt.plot(data3,data[i],label=data2[i])
    
    plt.legend(title='Time')
    plt.title('Progression of slump through time')
    plt.xlabel('Relative height up cylinder from bottom (x/a)')
    plt.ylabel(r'Relative density change ($\rho/\rho_{0}$)')
    plt.show()
    #plt.savefig('slumping_plot_progression.png')


PlotSlump()

