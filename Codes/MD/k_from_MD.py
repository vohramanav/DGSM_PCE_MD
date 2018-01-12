import numpy as np
import matplotlib
matplotlib.use('pdf')
from matplotlib import pyplot as plt
from matplotlib import mlab
from matplotlib import rc

import yaml
import argparse
import math

rc('text', usetex=True)

if __name__ == "__main__":

  parser = argparse.ArgumentParser()
  parser.add_argument('input_params',help='name of yml input file')
  parser.add_argument('energy_file',help='data file with exchange energies')
  parser.add_argument('Lr_file',help='data file with relaxed lengths')
  args = parser.parse_args()
  input_params = yaml.load(open(args.input_params).read())
  lc = 5.43 # Lattice constant in Ang
  dt = input_params["dt"]
  ld = input_params["ld"]
  W = input_params["w"]
  W = W*lc
  L = 50*lc
  cte = 1.5938
  eV_to_J = 1.602e-19
  ang_to_m = 1.0e-10
  ps_to_s = 1.0e-12
  energy_file = args.energy_file
  Lr_file = args.Lr_file
  f1 = open(energy_file,'r')
  data = f1.readlines()
  nr = np.array(data).shape[0]
  idx,e1,e2 = [],[],[]
  
  f2 = open(Lr_file,'r')
  data2 = f2.readlines()
  Lr = []

  k_data = np.zeros((nr,2))
  for line in data:
    p = line.split()
    idx.append(float(p[0]))
    e1.append(float(p[1]))
    e2.append(float(p[2]))
    idxv,e1v,e2v = np.array(idx),np.array(e1),np.array(e2)
  k_data[:,0] = idxv

  for line in data2:
    p = line.split()
    Lr.append(float(p[0]))
    Lrv = np.array(Lr)  
  
  e_avg,k = np.zeros((nr,1)),np.zeros((nr,1))
  e_avg[:,0] = [0.5*(abs(e1v[i])+abs(e2v[i])) for i in range(nr)]

  for i in range(nr):
   dTdx = float((cte*50.0))/float((Lrv[i]*lc))
   #dTdx = cte/lc
   dq = float((e_avg[i,0]*eV_to_J))/float(2.0*dt*ld*ps_to_s*W*W*pow(ang_to_m,2))
   k[i,0] = float(dq*ang_to_m)/float(dTdx)
  
  print (k)

  k_data[:,1] = k[:,0]
  np.savetxt('k_data.txt',k_data,fmt='%d %5.4f')

















