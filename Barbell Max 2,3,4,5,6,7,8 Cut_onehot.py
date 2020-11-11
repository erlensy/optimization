from qiskit import *
from qiskit.tools.monitor import job_monitor

import matplotlib as mpl
mpl.use('Agg')

import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.axes_grid1 import make_axes_locatable

import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
from qiskit.visualization import *
from scipy import optimize as opt
from qaoa import *
import os

import sys
sys.path.append('../')

from qiskit_utilities.utilities import *


from matplotlib import rc
font = {'size' : 36}
rc('font', **font);
rc('text', usetex=True)
pl.rcParams["figure.figsize"] = 12,8
pl.rcParams["axes.titlesize"] = 24
pl.rcParams["axes.labelsize"] = 36
pl.rcParams["lines.linewidth"] = 3
pl.rcParams["lines.markersize"] = 10
pl.rcParams["xtick.labelsize"] = 36
pl.rcParams["ytick.labelsize"] = 36
#pl.style.use('bmh')

V = np.arange(0,2,1)
E =[(0,1,1.0)]

G = nx.Graph()
G.add_nodes_from(V)
G.add_weighted_edges_from(E)

#pos = nx.spring_layout(G)
#nx.draw_networkx(G,apos=pos)

Aer.backends()
backend = Aer.get_backend('qasm_simulator')

beta_n = 12
gamma_n = 24
beta_max = np.pi/2
gamma_max = np.pi
optmethod='Nelder-Mead'
shots=1024*2*2*2
rerun=True

maxdepth=3

depths=range(1,maxdepth+1)

outstr=""

for k_cuts in [2,3]:
    for alpha in [0, 2*max(G.number_of_nodes()/k_cuts, k_cuts*G.number_of_edges()), None]:
        if alpha == None:
            circuit_version=2
        else:
            circuit_version=1

        outstr+="k="+str(k_cuts)
        outstr+="alpha="+str(alpha)
        outstr+="version="+str(circuit_version)
        print(" k=", k_cuts, " alpha=", alpha)
        name="Barbell"+str(gamma_n)+"x"+str(beta_n)+"_v"+str(circuit_version)+"_k"+str(k_cuts)+"_onehot_alpha_"+str(alpha)
        Elandscape, gammabetas, E, best =  runQAOA(G, k_cuts, backend, gamma_n, beta_n, gamma_max, beta_max, optmethod=optmethod, circuit_version=circuit_version, shots=shots, name=name, rerun=rerun, maxdepth=maxdepth, onehot=True, onehot_alpha=alpha)

        max_val=1

        shiftg=gamma_max/(2*gamma_n)
        shiftb=beta_max/(2*beta_n)

        pl.figure(figsize=(20,10));
        pl.clf()
        pl.imshow(Elandscape,interpolation='spline36',origin='lower'
                    ,extent=[-shiftg,gamma_max+shiftg,-shiftb,beta_max+shiftb], aspect=1)
        pl.xticks([0,gamma_max/2, gamma_max], ['0', r'$\pi$', r'$2\pi$'])
        pl.yticks([0,beta_max], ['0', r'$\pi/2$'])
        pl.xlabel('$\gamma$',loc='left')
        pl.ylabel(r'$\beta$')
        pl.colorbar(shrink=0.75, pad=0.05, orientation="horizontal")


        #pl.plot(gammabetas['x0_d1'][0], gammabetas['x0_d1'][1],'xw')
        pl.plot(gammabetas['xL_d1'][0], gammabetas['xL_d1'][1],'or')

        pl.tight_layout()
        pl.savefig("pics/E_"+name+".png")
        pl.close()

        for depth in depths:
            outstr+=" & "+"{:.3f}".format(E[str(depth)]/max_val)
        outstr+="\\\\ \n"
        print(outstr)