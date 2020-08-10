#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import os
import random
import sys
import subprocess
import shutil
from decimal import Decimal
#from scipy.linalg import orth
import time

epsilon0 = 8.854187817*10**(-12) #F/m
hbar = 6.582119514*10**(-16) #eV s
hbar2 = 1.054571800*10**(-34) #J s
mass = 9.10938356*10**(-31) # kg
c = 299792458 #m/s
e = 1.60217662*10**(-19) #C
pi = np.pi
kb = 8.6173303*(10**(-5)) #eV/K
amu = 1.660539040*10**(-27) #kg


def pega_geom(freqlog):
    atomos = []
    status = 0
    if ".log" in freqlog:
        busca = "Input orientation:"
        with open(freqlog, 'r') as f:
            for line in f:
                if "Standard orientation" in line:
                    busca = "Standard orientation:"
                    break 
        #print("\nEstou usando o", busca,"\n")           
        G = np.zeros((1,3))
        n = -1
        with open(freqlog, 'r') as f:
            for line in f:
                if "Optimized Parameters" in line:
                    status = 1
                if status == 1: 
                    if n < 0:
                        if busca in line:
                            n = 0
                        else:
                            pass
                    elif n >= 0 and n < 4:
                        n += 1
                    elif n >= 4 and "---------------------------------------------------------------------" not in line:    
                        line = line.split()
                        NG = []
                        for j in range(3,len(line)):
                            NG.append(float(line[j]))
                        atomos.append(line[1])
                        G = np.vstack((G,NG))       
                        n += 1  
                    elif "---------------------------------------------------------------------" in line and n>1:
                        break       
    else:
        G = np.zeros((1,3))
        with open(freqlog, 'r') as f:
            for line in f:
                line = line.split()
                try:
                    vetor = np.array([float(line[1]),float(line[2]), float(line[3])])
                    atomos.append(line[0])
                    G = np.vstack((G,vetor))
                except:
                    pass
    G = G[1:,:] #input geometry                 
    return G, atomos


 
def busca_input(freqlog):
    base = 'lalala'
    with open(freqlog, 'r') as f:
        for line in f:
            if "#" in line:
                line = line.split()
                for elem in line:
                    if "/" in elem:
                        base = elem
                        break
    return base                 

              
def gera_optcom(atomos,G,base,nproc,mem,omega,op):
    header = "%nproc=JJJ\n%mem=MEM\n# BASE iop(3/108=MMMMM00000) iop(3/107=MMMMM00000) "+op+"\n\nTITLE\n\n0 1\n"
    header = header.replace("JJJ",nproc).replace("MEM", mem).replace("BASE", base).replace("MMMMM",omega)
    with open("OPT_"+omega+"_.com", 'w') as f:
        f.write(header)
        for i in range(0,len(atomos)):
            texto = "%2s % 2.14f % 2.14f % 2.14f" % (atomos[i],G[i,0],G[i,1],G[i,2])
            f.write(texto+"\n")
        f.write("\n")


def gera_ioncom(atomos,G,base,nproc,mem,omega):
    header = "%nproc=JJJ\n%mem=MEM\n# BASE iop(3/108=MMMMM00000) iop(3/107=MMMMM00000)\n\nTITLE\n\n1 2\n"
    header = header.replace("JJJ",nproc).replace("MEM", mem).replace("BASE", base).replace("MMMMM",omega)
    with open("pos_"+omega+"_.com", 'w') as f:
        f.write(header)
        for i in range(0,len(atomos)):
            texto = "%2s % 2.14f % 2.14f % 2.14f" % (atomos[i],G[i,0],G[i,1],G[i,2])
            f.write(texto+"\n")
        f.write("\n")
    header = "%nproc=JJJ\n%mem=MEM\n# BASE iop(3/108=MMMMM00000) iop(3/107=MMMMM00000)\n\nTITLE\n\n-1 2\n"
    header = header.replace("JJJ",nproc).replace("MEM", mem).replace("BASE", base).replace("MMMMM",omega)
    with open("neg_"+omega+"_.com", 'w') as f:
        f.write(header)
        for i in range(0,len(atomos)):
            texto = "%2s % 2.14f % 2.14f % 2.14f" % (atomos[i],G[i,0],G[i,1],G[i,2])
            f.write(texto+"\n")
        f.write("\n")           

def watcher(file):
    nome = file[:-4]+".log"
    status = 0
    while status != 1:
        try:            
            with open(nome, 'r') as f:
                for line in f:
                    if "Normal termination of Gaussian" in line:
                        status = 1
            if status == 1:
                break
            else:
                time.sleep(60)
        except:
            pass

def pega_energia(file):
    energias = []
    with open(file, 'r') as f:
        for line in f:
            if "SCF Done:" in line:
                line = line.split()
                energias.append(float(line[4]))
    return min(energias)

def pega_homo(file):
    status = 0
    HOMOS = []
    with open(file, 'r') as f:
        for line in f:
            if ('OPT' in file and "Optimized Parameters" in line) or 'pos' in file or 'neg' in file:
                    status = 1
            if status == 1:
                if "occ. eigenvalues" in line:
                    line = line.split()
                    homos = line[4:]
                    HOMOS.extend(homos)
    if len(HOMOS) == 0:
        with open(file, 'r') as f:
            for line in f:
                if "occ. eigenvalues" in line:
                    line = line.split()
                    homos = line[4:]
                    HOMOS.extend(homos)
    HOMOS = list(map(float,HOMOS))
    return max(HOMOS)               

def rodar_omega(atomos,G,base,nproc,mem,omega,gauss,geomlog,op): 
    omega = str(omega)
    gera_optcom(atomos,G,base,nproc,mem,omega,op)
    files = [i for i in os.listdir(".") if ".com" in i and geomlog not in i]
    for file in files:
        try:
            subprocess.call(['tsp', gauss, file])  #os.system("tsp g09 "+file)
        except:
            subprocess.call(['ts', gauss, file])  #os.system("tsp g09 "+file)
    for file in files:  
        watcher(file)
    for file in files:  
        if op == 'opt':
            G, atomos = pega_geom(file[:-4]+".log")
        neutro = pega_energia(file[:-4]+".log")
        homo_neutro = pega_homo(file[:-4]+".log")
        gera_ioncom(atomos,G,base,nproc,mem,omega)
    files = [i for i in os.listdir(".") if ".com" in i and "pos" in i or "neg" in i]
    for file in files:
        try:
            subprocess.call(['tsp', gauss, file])  #os.system("tsp g09 "+file)
        except:
            subprocess.call(['ts', gauss, file])  #os.system("tsp g09 "+file)
    for file in files:
        watcher(file)
    for file in files:
        if "pos" in file:
            cation = pega_energia(file[:-4]+".log")
        elif "neg" in file:
            anion = pega_energia(file[:-4]+".log")
            homo_anion = pega_homo(file[:-4]+".log")        
    try:
        os.mkdir("Logs")
    except:
        pass
    os.system("mv *.com Logs")
    os.system("mv *.log Logs")
    J = (homo_neutro + cation - neutro)**2 + (homo_anion + neutro - anion)**2
    return J, G, atomos


try:
    geomlog = sys.argv[1]
except:
    print("#                       #     #")
    print("#        ######   ####  #  #  #")
    print("#        #       #    # #  #  #")
    print("#        #####   #    # #  #  #")
    print("#        #       #    # #  #  #")
    print("#        #       #    # #  #  #")
    print("#######  ######   ####   ## ## ")
    print("----HIDEKI AJUDOU UM POUCO----\n")
    print("Modo de usar:\n")
    print("lw geom nproc mem w_i passo g09_ou_g16 relax(y/n) &")
    print('\n')
    sys.exit()



base = busca_input(geomlog)
G, atomos  = pega_geom(geomlog)
print(geomlog,'\n')
for i in range(len(atomos)):
    print(atomos[i],G[i,:])   
nproc = sys.argv[2]
mem = sys.argv[3]
omega1 = sys.argv[4]
passo = sys.argv[5]
gauss = sys.argv[6]
relax = sys.argv[7]

if relax == 'y':
    op = 'opt'
else:
    op = ''

try:
    int(nproc)
except:
    print("nproc tem que ser um inteiro, animal!")
    sys.exit()

try:
    passo = int(float(passo)*10000)
except:
    print("Passo tem que ser um inteiro, animal!\n")
    sys.exit()

try:
    omega1 = float(omega1)*10000
except:
    print("Valor inviável de omega!\n")
    sys.exit()
omega1 = format(int(omega1), '05')
omegas, Js = [], []
oms, jotas = [], []
try:
    with open("omega.lw", 'r') as f:
        for line in f:
            line = line.split()
            if len(line) == 2:
                om = format(int(line[0]), '05')
                omegas.append(om)
                Js.append(line[1])
except:
    pass


while passo > 50:#max_iter:
    if omega1 in omegas:
        ind = omegas.index(omega1)
        J = Js[ind]
    else:
        J, G, atomos = rodar_omega(atomos,G,base,nproc,mem,omega1,gauss,geomlog,op)              
        omegas.append(omega1)
        Js.append(J)  
    oms.append(omega1)
    jotas.append(J)    
    try:
        if jotas[-1] - jotas[-2] > 0:     
            passo = int(passo/2)
            sign = -1*np.sign(int(oms[-1]) - int(oms[-2])) 
        else:           
            try:
                sinais = np.asarray(oms[1:]) - np.asarray(oms[:-1])
                sinais = [np.sign(i) for i in sinais]
                if len(set(sinais[-4:])) == 1:
                    passo = int(2*passo)
            except:
                pass
            sign = +1*np.sign(int(oms[-1]) - int(oms[-2]))
        omega1 = int(omega1) + sign*passo
        omega1 = format(int(omega1), '05')    
    except:
        omega1 = int(omega1) + passo
        omega1 = format(int(omega1), '05')
    with open("omega.lw", 'w') as f:
        list1, list2 = zip(*sorted(zip(omegas, Js))) #(omegas[t] for t in zip(*sorted(zip(omegas, Js))))
        for i in range(len(list1)):
            f.write(str(list1[i])+"    "+str(list2[i])+"\n")
        f.write("\nBest value so far: "+str(list1[list2.index(min(list2))])+"\n")    


with open("omega.lw", 'w') as f:
        list1, list2 = zip(*sorted(zip(omegas, Js))) #(omegas[t] for t in zip(*sorted(zip(omegas, Js))))
        for i in range(len(list1)):
            f.write(str(list1[i])+"    "+str(list2[i])+"\n")
        f.write("\nOptimized value: "+str(list1[list2.index(min(list2))])+"\n")
        f.write("\nC'est fini!")


        
        
        


    
        
