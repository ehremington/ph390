#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Exercise 17.3                                                      *
%*     filename: ch17pr03.m                                               *
%*     program listing number: 17.3                                       *
%*                                                                        *
%*     This program simulates two-dimensional percolation.                *
%*     Hoshen and Kopelman algorithm is used for cluster labeling.        *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  03/17/2017.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

# ease all previous figures
plt.close('all')

N=32
# p=0.59 is close to the transition point
p=0.59

# Graphics setting
plt.figure(figsize=(6,6))
plt.axis('equal')
plt.axes(xlim=(-1, N), ylim=(-1, N))

# place atoms at random
r=np.random.rand(N,N)
lattice = r < p

# draw clusters
for i in range(0,N):
    for j in range(0,N):
        if lattice[i,j]:
            circle=plt.Circle((i,j),0.5,fc='b')
            plt.gca().add_patch(circle)
            plt.pause(0.0001)

# Labeling: 1st pass (Making initial labels and map)
label=np.zeros([N,N],dtype=np.int)  # allocate array
remap=np.zeros(N*N,dtype=np.int)

new=0
for j in range(0,N):
    for i in range(0,N):
        
        i1 = i-1
        j1 = j-1
        
        if lattice[i,j]>0:
            if i1 < 0:
                left=0            # outside the box
            else:
                left=label[i1,j]  # left neighbor

            if j1 < 0:
                down=0            # outside the box.
            else:
                down=label[i,j1]  # down neighbor
            
            if down==0 and left==0:  # if both are unocupied
                new=new+1            # create a new cluster.
                label[i,j]=new
                remap[new]=new
                
            elif down*left>0:       # both are occupied.
                
                if down==left:      # if they belong to the same cluster.
                    label[i,j]=left # join to the cluster.
                    
                else:               # connecting two different clusters.
                    found = False
                    while not(found):
                        if remap[left]==left:
                            found = True
                        else:
                            left=remap[left]

                    found = False
                    while not(found):
                        if remap[down]==down:
                            found = True
                        else:
                            down=remap[down]
                    
                    if left==down:     # they again belong to the same
                        label[i,j]=left
                    else:
                        nmax=np.max([left,down])
                        nmin=np.min([left,down])
                        label[i,j]=nmin    # coalesce two clusters
                        remap[nmax]=nmin   # add to the chain
                
            elif down>0:            # only down neighbor is occupied
                label[i,j] = down   # join to the neighbor
                
            else:                   # only the left neighbor is occupied
                label[i,j] = left   # join to the neighbor

# Labeling: 2nd pass (Collapse the label in the same cluster)

nmax = np.max(label)
for i in range(nmax,0,-1):
    label[label==i]=remap[i]

# Labeling: 3rd pass (Make the label continuous)
#          This procedure is not essential.
j=0
for i in range(1,nmax+1):
    if remap[i]==i:
        j=j+1
        label[label==i] = j

# Find cluster size and find percolation
nmax = np.max(label)
size = np.zeros(nmax+1,dtype=np.int)
maxsize=0
largest=1
print('Cluster   Size  Percolation')
for i in range(1,nmax+1):
    percx = any(label[0,:]==i)*any(label[N-1,:]==i) >0
    percy = any(label[:,0]==i)*any(label[:,N-1]==i) >0
    size[i]= (label==i).sum()
    if size[i] > maxsize:
        largest = i
        maxsize = size[i]

    if percx or percy:
        perc='YES'
    else:
        perc=' NO'

    print('{0:5d}:  {1:5d},    {2:s}'.format(i, size[i], perc))

# Show the largest cluster
for i in range(0,N):
    for j in range(0,N):
        if label[i,j] == largest:
            circle=plt.Circle((i,j),0.5,fc='r')
            plt.gca().add_patch(circle)
            plt.pause(0.0001)

# Plot size distribution
plt.figure(figsize=(6,5))
plt.hist(size,maxsize)
plt.show()