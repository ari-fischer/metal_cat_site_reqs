#!/usr/bin/env python
# coding: utf-8

from ase import Atoms
import numpy as np
import pandas as pd
import copy

def center_ads(atoms):
    #import ASE atoms structure and export the adsorbate moiety (C, H, O) with coordinates centered to origin
    
    #extract the atoms constituting the adsorbate (H, C, O atoms)
    atoms_H = atoms[np.array(list(atoms.symbols))=='H']
    atoms_C = atoms[np.array(list(atoms.symbols))=='C']
    atoms_O = atoms[np.array(list(atoms.symbols))=='O']

    #build an object of just the hydrocarbon moiety from the H and C with periodic boundary conditions
    atoms_ads= atoms_H+atoms_C+atoms_O
    atoms_ads.set_pbc(True)
    
    #turn atoms_ads into dataframe with the atom symbol and the xyz coordinates
    #requires cartesian coordinates
    df_atoms = pd.DataFrame(list(atoms_ads.symbols),columns=['atom'])

    df_pos = pd.DataFrame(atoms_ads.get_positions(),columns=['x','y','z'])
    df = pd.concat([df_atoms,df_pos], axis=1)
    
    #build an adsorbate centered at the origin
    # get center of position of the adsorbate as calculated
    COP = [df['x'].mean(),df['y'].mean(),df['z'].mean()]

    #build two matrices of C and H positions, relative to the adsorbate COP 
    df_H = df[df['atom']=='H']
    df_C = df[df['atom']=='C']
    df_O = df[df['atom']=='O']

    #shift coordinates to the origin
    M_H = df_H[['x','y','z']].to_numpy()-COP
    M_C = df_C[['x','y','z']].to_numpy()-COP
    M_O = df_O[['x','y','z']].to_numpy()-COP

    #build dataframe combining adsorbate positions
    M_ads = np.vstack([M_H, M_C, M_O])

    return M_H, M_C, M_O, M_ads, COP



def build_ads_grid(M_ads,d,grid_ex):
    # build a 3D grid and populate it with a sphere with the vdW volume

    # get the vdw volume of C and H from Bondi reference data
    #Cs vdw volume
    vdW_R_C = 1.7 #radius, Ang
    vdW_V_C = np.pi*4/3*vdW_R_C**3 #vol of single C

    #Hs vdw volume
    vdW_R_H = 1.2 #radius, Ang
    vdW_V_H = np.pi*4/3*vdW_R_H**3 #vol of single H

    #Os vdw volume
    vdW_R_O = 1.52 #radius, Ang
    vdW_V_O = np.pi*4/3*vdW_R_O**3 #vol of single O

    # build a grid that encompasses the min and max of the adsorbate in x,y, and z directions
    # add "grid_ex" Ang above and below the maximum position of the structure extends beyond the vdW radii of C and H)

    # calculate the total length of each grid dimension
    c_mins = np.min(M_ads, axis=0)-grid_ex
    c_maxs = np.max(M_ads, axis=0)+grid_ex
    c_ranges = c_maxs-c_mins

    # build the matrix by normalizing the cell lengths by the element size (d)
    M_size = np.ceil(c_ranges/d)
    M_size = (np.ceil(M_size/2))*2 #make sure they are all even

    # build an empty 3D matrix to fill with vdW spheres with M_size dimensions
    M_fill = np.zeros([int(M_size[0]),int(M_size[1]),int(M_size[2])])
    return M_size,M_fill


def build_atom_volume(R,d):
    #construct a matrix to represent a sphere with radius R on a grid with unit size d

    #construct A_fill as the matrix to enclose the sphere
    A_size = [int(np.ceil((R+0.4)/d))*2,
              int(np.ceil((R+0.4)/d))*2,
             int(np.ceil((R+0.4)/d))*2]
    A_fill = np.zeros(A_size)

    #iterate over the dimensions of A_size and fill elements with 1 if they overlap with a sphere centered about the origin
    
    #initialize the z-array
    zs = np.array([])
    
    for i in range(int(A_size[2]/2)):
        # Initialize an empty array to store the z-values for this slice
        z_slice = np.array([])

        # Loop over y-coordinates
        for j in range(int(A_size[1]/2)):
            # Assign current x and y values
            x = i
            y = j
            # Calculate the positive z-value for the current (x, y) using the sphere equation
            z_plus = +np.sqrt((R/d)**2-x**2-y**2)
            # Calculate the negative z-value for the current (x, y)
            z_minus = -np.sqrt((R/d)**2-x**2-y**2)
            # Append the positive z-value to the z_slice array (for now, only considering the upper half of the sphere)
            z_slice = np.append(z_slice,[z_plus])

        # Create a symmetrical slice by appending the flipped z_slice (which represents the lower half of the sphere)
        z_slice = np.append(np.flip(z_slice),z_slice)

        # Loop over all x-coordinates in the grid
        for k in range(A_size[0]):
            # Calculate the maximum z index for filling
            z_max = np.ceil(z_slice[k])
            # If the z_max value is positive, fill the grid at the appropriate positions with 1s
            if z_max>0:
                # Fill the upper half of the sphere in the grid
                A_fill[i+int(A_size[0]/2),k,int(A_size[0]/2-z_max-1):int(z_max-1+A_size[0]/2)]=1
                # Fill the lower half of the sphere in the grid (by symmetry)
                A_fill[int(A_size[0]/2)-i,k,int(A_size[0]/2-z_max-1):int(z_max-1+A_size[0]/2)]=1

    #compare numerical volume with geometry calculation to assess accuracy
    print('grid volume: '+str(np.sum(A_fill)*d**3)) 
    print('vdW volume: '+str(4*np.pi/3*R**3)) 
    return A_size, A_fill

def mol_2_vol_area(M_fill,pos_A,A_size,A_fill,d):
    #fills a matrix enclosing the molecular volume with ones at locations corresponding to the van der Waals spheres surrounding atomic coordinates

    #takes in the matrix to fill with molecule volume (M_fill),
    #the coordinates of the atoms to fill (pos_A),
    #the shape of the matrix representing a single atom (A_size)
    #the matrix for the atom used to fill the molecular volume (A_fill)
    #and the unit size for grid elements (d)

    #outputs the filled matrix (M_fill) and the area projection (vdw_area)

    #loop over the atoms used to fill the volume and input them into the molecular volume matrix 
    for i in range(np.shape(pos_A)[0]):
        M_fill[int(np.ceil(pos_A[i][0])-int(A_size[0]/2)):int(np.ceil(pos_A[i][0])+int(A_size[0]/2)),
           int(np.ceil(pos_A[i][1])-int(A_size[0]/2)): int(np.ceil(pos_A[i][1])+int(A_size[0]/2)),
           int(np.ceil(pos_A[i][2])-int(A_size[0]/2)): int(np.ceil(pos_A[i][2])+int(A_size[0]/2))] += A_fill

    #replace all of the values larger than 1 (reflecting overlapping spheres) with 1
    M_fill[M_fill>1]=1
    #collapse the z-axis of the 3D geomertry down onto the xy plane
    vdW_2d = np.sum(M_fill,axis=2)
    #replace any non-zero element with 1
    vdW_2d[vdW_2d>1] = 1
    #count the filled grid elements and multiply by grid size for the area element to get the vdW area
    vdW_area=np.sum(vdW_2d)*d**2
    return M_fill,vdW_area


def mol_overlap_2D(M0,pos_A,A_size,A_fill,d):
    #takes in the filled matrix and places probe molecule at nearby to sample if there is overlap.
    #M_fill -- matrix of 0 and 1s of where the molecule fills the volume
    #pos_A -- the trial location to add probe molecule
    #add A at the location of pos_A
    M = copy.deepcopy(M0)
    
    M[int(np.ceil(pos_A[0])-int(A_size[0]/2)):int(np.ceil(pos_A[0])+int(A_size[0]/2)),
    int(np.ceil(pos_A[1])-int(A_size[0]/2)): int(np.ceil(pos_A[1])+int(A_size[0]/2))] += A_fill

    overlap = len(np.where(M==2)[0])>1
    
    return overlap