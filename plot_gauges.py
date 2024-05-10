import matplotlib.pyplot as plt
import numpy as np
import os.path
import h5py

thetis_output_dirs = [
        'outputsB_dt0.01_alpha0.01', 
        'outputsE_dt0.01_alpha0.01',
        'outputsF_dt0.01_alpha0.01'] 
hdf5_name = 'diagnostic_timeseries_gauge{}_elev.hdf5'
gauge_csv = 'raw_data/WaveGages.csv'
data = np.loadtxt(gauge_csv, delimiter=',')

for i in range(1,4):
    plt.figure()
    for dir in thetis_output_dirs:
        f = h5py.File(os.path.join(dir, hdf5_name.format(i)), 'r')
        plt.plot(f['time'], f['elev'], label=dir)
    plt.plot(data[:,0], data[:,i]/100, label=f'lab')
    plt.legend()
    plt.xlim(0, 25)
    plt.savefig(f'gauge{i}.png')
