#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 11:28:01 2017
@author: alexandre

Argumentos:
argv[1]: modelo = {bf, bn, bfn}
argv[2]: cenário = {0: bf_classicdiffusion, 1: bf_gdiffusion, 2: bf_nonlineardiffusion, 3: bn_small_l, 4: bn_small_s_and_l,
5: bn_normal, 6:bfn_normal, 7: bfn_small_s_and_l, 8: bfn_small_k}
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.pyplot import imshow
from scipy.interpolate import interp2d

def getMatrixFromFile(filename, u):
    try:
        file_ = open(filename, 'r');
        t = 0;
        x = 0;
        y = 0;
        for content in file_.readlines():
            y = 0
            if (content.startswith('--')): #increment time step
                t += 1
                x = 0;
                y = 0;
            else:
                content = content.rstrip()
                content = content.split(' ')
                #print('content:'+str(content))
                for v in content:
                    u[t][x][y] = float(v)
                    #print('v:' + v)
                    y += 1;
                x += 1
        file_.close();
        return u
    except (OSError, IOError) as exc:
        print("Your file could not be open to, your exception details are: {}".format(exc))
        return None

def getMaxMinValueFrom3DMatrix(u):
    maxvalue = 0
    minvalue = 2
    for i in range(len(u)):
        for j in range(len(u[i])):
            for k in range(len(u[i][j])):
                if u[i][j][k] > maxvalue:
                    maxvalue = u[i][j][k]
                if u[i][j][k] < minvalue:
                    minvalue = u[i][j][k]
    return [maxvalue,minvalue]

def getArrayFromFile(filename):
    a = []
    with open(filename) as afile:
        a = afile.readlines()
    a = [content.strip() for content in a]
    a = [float(i) for i in a]
    return a

def initializeGraphParameters():
    global scalarMap, cm;
    cm = plt.get_cmap('jet')
    cNorm  = colors.Normalize(vmin=0, vmax=10)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)

def createFigure(title, x, u, maxDensity, minDensity):
    global cm
    fig = plt.figure(figsize=(10,7));
    ax = fig.gca();
    plt.title(title, fontsize=17);
    plt.xlabel('x',fontsize=17);
    plt.ylabel('y',fontsize=17);
    plt.axis([min(x), max(x), min(x), max(x)])
    pcm = plt.pcolormesh(x, x, u[t], cmap='jet', vmin=minDensity, vmax=maxDensity)
    plt.colorbar(pcm)
    return [fig,ax];

def generateGraphs(fig, filename):
    #pp = PdfPages(filename + '.pdf')
    #pp.savefig(fig,bbox_inches='tight')
    #pp.close()
    fig.savefig(filename + '.svg', format='svg', dpi=1200,bbox_inches='tight')


#Remove previous images
svgimagelist = [ f for f in os.listdir(".") if f.endswith(".svg") ]
for f in svgimagelist:
    os.remove(f)
epsimagelist = [ f for f in os.listdir(".") if f.endswith(".eps") ]
for f in epsimagelist:
    os.remove(f)
pdfimagelist = [ f for f in os.listdir(".") if f.endswith(".pdf") ]
for f in pdfimagelist:
    os.remove(f)

#Global constants and parameters
maxDensity = 1.0
deltaY = 0.1
deltaX = 0.1
#size = int(5/deltaX) + 1
size = 51
#Create space and time axis values
for files in os.listdir("."):
    if files.startswith("x.dat"):
        xfile = files
    elif files.startswith("t.dat"):
        tfile = files
    elif files.startswith("bacteria.dat"):
        bfile = files
    elif files.startswith("coa_vwbp.dat"):
        coafile = files
    elif files.startswith("fibrin.dat"):
        fibrinfile = files
    elif files.startswith("neutrophil.dat"):
        nfile = files
    elif files.startswith("deadneutrophil.dat"):
        deadnfile = files
    elif files.startswith("toxin.dat"):
        toxinfile = files
    elif files.startswith("mr.dat"):
        mrfile = files
    elif files.startswith("ma.dat"):
        mafile = files
    elif files.startswith("nd.dat"):
        ndfile = files
    elif files.startswith("dmt.dat"):
        dmtfile = files


xis = getArrayFromFile(xfile)
t = getArrayFromFile(tfile)
T = len(t)

print('x:'+ str(xis))
print('t:'+ str(t))
print('T:'+ str(T))

initializeGraphParameters()
color_index = 0
colorVal = scalarMap.to_rgba(color_index)
"""
abs = abscess
ce = complete elimination
co = colony formation
inf = disseminated infection
def = defficient response
sp = small phagocytosis
"""
scenario = str(sys.argv[1])
b_filename = 'bacteria_' + scenario + '_'
coa_filename = 'coa_vwbp_' + scenario + '_'
n_filename = 'neutrophil_' + scenario + '_'
f2_filename = 'f2_' + scenario + '_'
fib_filename = 'fibrin_' + scenario + '_'
to_filename = 'toxin_' + scenario + '_'
mr_filename = 'mr_' + scenario + '_'
ma_filename = 'ma_' + scenario + '_'
nd_filename = 'nd_' + scenario + '_'
dmt_filename = 'dmt_' + scenario + '_'

b = [[[0. for y in range(0,size)] for x in range(0,size)] for t in range(0,T)]
coa = [[[0. for y in range(0,size)] for x in range(0,size)] for t in range(0,T)]
fb = [[[0. for y in range(0,size)] for x in range(0,size)] for t in range(0,T)]
n = [[[0. for y in range(0,size)] for x in range(0,size)] for t in range(0,T)]
to = [[[0. for y in range(0,size)] for x in range(0,size)] for t in range(0,T)]
nd = [[[0. for y in range(0,size)] for x in range(0,size)] for t in range(0,T)]
mr = [[[0. for y in range(0,size)] for x in range(0,size)] for t in range(0,T)]
ma = [[[0. for y in range(0,size)] for x in range(0,size)] for t in range(0,T)]
#d = [[[0. for y in range(0,size)] for x in range(0,size)] for t in range(0,T)]
dmt = [[[0. for y in range(0,size)] for x in range(0,size)] for t in range(0,T)]

getMatrixFromFile(bfile,b)
getMatrixFromFile(coafile,coa)
getMatrixFromFile(fibrinfile,fb)
getMatrixFromFile(nfile,n)
getMatrixFromFile(toxinfile,to)
getMatrixFromFile(mrfile,mr)
getMatrixFromFile(mafile,ma)
getMatrixFromFile(ndfile,nd)
getMatrixFromFile(dmtfile,dmt)
for t in range(0,T):
    [maxv,minv] = getMaxMinValueFrom3DMatrix(b);
    [fig1,ax1] = createFigure('Bacteria',xis,b,maxv,minv);
    [maxv,minv] = getMaxMinValueFrom3DMatrix(coa);
    [fig2,ax2] = createFigure('Coa/vWbp',xis,coa,maxv,minv);
    [maxv,minv] = getMaxMinValueFrom3DMatrix(n);
    [fig3,ax3] = createFigure('Neutrophils',xis,n,maxv,minv);
    [maxv,minv] = getMaxMinValueFrom3DMatrix(fb);
    [fig4,ax4] = createFigure('Fibrin',xis,fb,maxv,minv);
    [maxv,minv] = getMaxMinValueFrom3DMatrix(nd);
    [fig5,ax5] = createFigure('Apoptotic Neutrophils',xis,nd,maxv,minv);
    [maxv,minv] = getMaxMinValueFrom3DMatrix(to);
    [fig6,ax6] = createFigure('Toxins',xis,to,maxv,minv);
    [maxv,minv] = getMaxMinValueFrom3DMatrix(mr);
    [fig7,ax7] = createFigure('Resting Macrophages',xis,mr,maxv,minv);
    [maxv,minv] = getMaxMinValueFrom3DMatrix(ma);
    [fig8,ax8] = createFigure('Activated Macrophages',xis,ma,maxv,minv);
    [maxv,minv] = getMaxMinValueFrom3DMatrix(dmt);
    [fig9,ax9] = createFigure('Tissue damage',xis,dmt,maxv,minv);
    ax1.legend(loc='upper right')
    ax2.legend(loc='upper right')
    ax3.legend(loc='upper right')
    ax4.legend(loc='upper right')
    ax5.legend(loc='upper right')
    ax6.legend(loc='upper right')
    ax7.legend(loc='upper right')
    ax8.legend(loc='upper right')
    ax9.legend(loc='upper right')
    generateGraphs(fig1,b_filename + str(t));
    generateGraphs(fig2,coa_filename + str(t));
    generateGraphs(fig3,n_filename + str(t));
    generateGraphs(fig4,fib_filename + str(t));
    generateGraphs(fig5,nd_filename + str(t));
    generateGraphs(fig6,to_filename + str(t));
    generateGraphs(fig7,mr_filename + str(t));
    generateGraphs(fig8,ma_filename + str(t));
    generateGraphs(fig9,dmt_filename + str(t));
    #Plotar bacteria n tempo calculando a soma dos elementos da matriz
