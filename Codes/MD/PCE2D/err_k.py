# This code is used to generate error estimates in thermal conductivity from input 
# energy exchange data. The data would be used as input to uqlab for constructing the
# response surface of the error as a function of length and thermal gradient. 

import numpy as np
import matplotlib
matplotlib.use('pdf')
from matplotlib import pyplot as plt
from matplotlib import mlab
from matplotlib import rc
from matplotlib.patches import Circle

import yaml
import argparse
import math

rc('text', usetex=True)

if __name__ == "__main__":

  parser = argparse.ArgumentParser()
  parser.add_argument('input_params',help='name of yml input file')
  parser.add_argument('energy_file',help='data file with exchange energies')
  args = parser.parse_args()
  input_params = yaml.load(open(args.input_params).read())
  lc = 5.43 # Lattice constant in Ang
  dt = input_params["dt"]
  ld = input_params["ld"]
  km = input_params["km"]
  w = input_params["w"]
  nlen = input_params["nlen"]
  ncte = input_params["ncte"]
  w = w*lc
  eV_to_J = 1.602e-19
  ang_to_m = 1.0e-10
  ps_to_s = 1.0e-12
  energy_file = args.energy_file
  
  f1 = open(energy_file,'r')
  data = f1.readlines()
  nr = np.array(data).shape[0]
  e1,e2,l,dTdx = [],[],[],[]
  samples_qoi,x_data,y_data = np.zeros((nr,3)),np.zeros((nr,2)),np.zeros((nr,1))
  for line in data:
    p = line.split()
    e1.append(float(p[0]))
    e2.append(float(p[1]))
    l.append(float(p[2])) # (actual length)/LC in Ang
    dTdx.append(float(p[3])) # K/Ang
  e1v,e2v,lv,dTdxv = np.array(e1),np.array(e2),np.array(l),np.array(dTdx)
  # reorganize e1v,e2v,lv and dTdxv
  lvo,dTdxvo = np.zeros((nr,1)),np.zeros((nr,1))
  e1vo,e2vo = np.zeros((nr,1)),np.zeros((nr,1))
  co = 0
  for i in range(ncte):
    e1vo[co:co+nlen,0] = [e1v[i+k*ncte] for k in range(nlen)]
    e2vo[co:co+nlen,0] = [e2v[i+k*ncte] for k in range(nlen)]
    lvo[co:co+nlen,0] = [lv[i+k*ncte] for k in range(nlen)]
    dTdxvo[co:co+nlen,0] = [dTdxv[i+k*ncte] for k in range(nlen)]
    co += nlen

  for i in range(nlen*ncte):
    print (e1vo[i,0],e2vo[i,0],lvo[i,0],dTdxvo[i,0])

  samples_qoi[:,0],samples_qoi[:,1] = lvo[:,0],dTdxvo[:,0]
  lvo[:,0],dTdxvo[:,0] = lvo[:,0]*lc,dTdxvo[:,0]/float(lc)
  e_avg,k,err_k = np.zeros((nr,1)),np.zeros((nr,1)),np.zeros((nr,1))
  e_avg[:,0] = [0.5*(abs(e1vo[i])+abs(e2vo[i])) for i in range(nr)]

  for i in range(nr):
    dq = float((e_avg[i,0]*eV_to_J))/float(2.0*dt*ld*ps_to_s*w*w*pow(ang_to_m,2))
    k[i,0] = float(dq*ang_to_m)/float(dTdxvo[i])
#
  #print (k)
  err_k = abs(k-km)
  samples_qoi[:,2] = err_k[:,0]
  #samples_qoi[:,2] = k[:,0]
  y_data[:,0] = samples_qoi[:,2]
  np.savetxt('samples_qoi.txt',samples_qoi,fmt='%5.4f')
  #np.savetxt('x_data.dat',x_data,fmt='%5.4f')
  xyz = np.zeros((nr,3))
  xyz[:,0] = (lvo[:,0])/float(np.amax(lvo[:,0]))
  xyz[:,1] = (dTdxvo[:,0])/float(np.amax(dTdxvo[:,0]))
  xyz[:,2] = 0.05*(err_k[:,0]-np.amin(err_k))/float(np.amax(err_k)-np.amin(err_k))
  mi = np.argmin(xyz[:,2])
  #y_data[:,0] = xyz[:,2]
  #np.savetxt('y_data.dat',y_data,fmt='%5.4f')
  #print (samples_qoi)


  fig = plt.figure()
  ax = fig.add_subplot(111)
  #plt.scatter(lv,dTdxv,s=k)
  for i in range(nr):
    if (xyz[i,2] > 0.02):
      ax.add_artist(Circle(xy=(xyz[i,0],xyz[i,1]), radius=0.7*xyz[i,2], color='b'))
      ax.text(xyz[i,0]-0.01,xyz[i,1]-0.01,r'$\mathrm{%3.2f}$' %(err_k[i,0]),color='w',fontsize=7)
    else:
      ax.add_artist(Circle(xy=(xyz[i,0],xyz[i,1]), radius=1.4*xyz[i,2], color='b'))
      
  #plt.scatter(xyz[mi,0],xyz[mi,1],s=30,marker='*',c='r')
  ax.add_artist(Circle(xy=(xyz[mi,0],xyz[mi,1]), radius=1.4*xyz[nr-1,2], color='b'))
  #ax.add_artist(Circle(xy=(xyz[mi,0],xyz[mi,1]), radius=0.015, color='r'))
  ax.set_xlim([np.amin(xyz[:,0])-0.1,np.amax(xyz[:,0])+0.1])
  ax.set_ylim([np.amin(xyz[:,1])-0.1,np.amax(xyz[:,1])+0.1])
  ax.set_xlabel(r'$\mathrm{L/L_{max}}$',fontsize=18)
  ax.set_ylabel(r'$\mathrm{\frac{dT}{dt}/\frac{dT}{dt}|_{max}}$',fontsize=18)
  #x0,x1 = ax.get_xlim()
  #y0,y1 = ax.get_ylim()
  #ax.set_aspect((x1-x0)/(y1-y0))
  #ax.grid(color='k', linestyle='--', alpha=0.5,linewidth=1)
  ax.tick_params(labelsize=18)
  plt.gca().set_aspect('equal', adjustable='box')
  plt.title (r'$\mathrm{\epsilon_k~=~k_{m} - k_{MD}~ at~300~K}$',fontsize=18)
  fig.savefig('k_plot.pdf')
