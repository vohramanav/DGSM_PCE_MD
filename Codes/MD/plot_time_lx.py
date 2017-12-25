import argparse
import numpy as np
import matplotlib
matplotlib.use('pdf')
from matplotlib import pyplot as plt
from matplotlib import mlab
from matplotlib import rc

rc('text', usetex=True)

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('log_sw_file',help='name of sw log-file')
  args = parser.parse_args()
  f1 = open(args.log_sw_file,'r')
  f2 = open(args.log_sw_file,'r') # total number of rows in the log file

# Extract lx from the log file
  data_1 = f2.readlines()
  lx_npt = np.zeros((50,2))
  for i in range(np.array(data_1).shape[0]):
    row = f1.readline()
    if "{tnpt}" in row:
      for n_skip in range(4):
        line_skip = f1.readline()
      for k in range(lx_npt.shape[0]):
        lx_data = f1.readline()
        p  = lx_data.split()
        lx_npt[k,0] = float(p[0])
        lx_npt[k,1] = float(p[1])

  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(lx_npt[:,0],lx_npt[:,1],lw=2.0)
  ax.set_xlabel(r'$\mathrm{Number~of~Steps}$',fontsize=15)
  ax.set_ylabel(r'$\mathrm{L_x (\AA)}$',fontsize=15)
  fig.savefig('lx_npt.pdf')

if __name__ == "__main__":
  main()
