#!/usr/bin/python

import os, re
import numpy as np 
import argparse

#This needs numpy version 1.8.0 to work; indexing errors seen with version 1.16.0

''' 
takes an ascii vtk file
rotates and scales to place it approximately centered to another ascii vtk file
change filenames in main 

'''

def rotation_matrix_x(theta):
    return np.array([[1,             0,              0], 
                     [0, np.cos(theta), -np.sin(theta)],  
                     [0, np.sin(theta),  np.cos(theta)]])


def rotation_matrix_y(theta):
    return np.array([[np.cos(theta),  0, -np.sin(theta)], 
                     [            0,  1,             0 ], 
                     [np.sin(theta),  0,  np.cos(theta)]])


def rotation_matrix_z(theta):
    return np.array([[np.cos(theta), -np.sin(theta), 0],
                     [np.sin(theta),  np.cos(theta), 0],
                     [            0,              0, 1]])



def position_vertices(orig_file, new_file, ring_file, scaling_orig=1.0, R_0=np.eye(3)):
    '''
    orig_file      vtk file to position 
    new_file       name of new vtk file, will be overwritten 
    ring_file      positions to a ring specified by this vtk file 
    scaling_orig   scales vertices by scalar 
    R_0            Initial rotation matrix 
    '''

    vertices_working = read_vtk(orig_file)
    vertices_working *= scaling_orig

    # putting the data aligned with here 
    vertices_ring = read_vtk(ring_file)

    tol = 1.0e-14; 

    centroid = np.mean(vertices_ring, axis=1)

    if vertices_ring.shape[0] is not 3:
        raise ValueError('need three dimensional array for ring')

    n_pts_ring = vertices_ring.shape[1]    

    normal = np.cross(vertices_ring[:,0] - centroid, vertices_ring[:, np.floor(n_pts_ring/4)] - centroid);       

    # require up facing normal on z component 
    if normal[2] < 0.0: 
        normal *= -1

    # normalize 
    normal = normal / np.linalg.norm(normal)

    # initial rotation if requested, default is identity
    vertices_working = np.dot(R_0, vertices_working)

    # set vertices to be centered at origin 
    # first to all min zero 
    for component in range(3):
        vertices_working[component,:] -= np.min(vertices_working[component,:])
    # x,y centered on origin, z still with minimum at zero 
    for component in range(2):
        vertices_working[component,:] -= np.mean(vertices_working[component,:])    

    # check inverse rotation first 
    # rotate normal so it has zero y component 
    phi = np.arctan2(normal[1], normal[2])
    R_x = rotation_matrix_x(phi)

    normal_no_y = np.dot(R_x, normal) 
    if abs(normal_no_y[1]) > tol:
        raise ValueError('did not remove y component')

    # then again so it has zero x 
    theta = np.arctan2(normal_no_y[0], normal_no_y[2])
    R_y = rotation_matrix_y(theta)

    normal_zdir_from_rotation = np.dot(R_y,normal_no_y)

    initial_normal = (np.array([0, 0, 1])).T 

    if np.linalg.norm(normal_zdir_from_rotation - initial_normal) > tol:
        raise ValueError('inverse rotation incorrect')

    # rotation to apply is the inverse of these 
    # and note that all matrices are orthogonal 
    R = np.dot(R_x.T, R_y.T)   

    if np.linalg.norm(np.dot(R,initial_normal) - normal) > tol:
        raise ValueError('R*[0;0;1] does not give correct normal')

    if np.linalg.norm(np.dot(R,R.T) - np.eye(3)) > tol:
        raise ValueError('rotation matrix is not orthogonal') 

    # apply rotation matrix to each vertex 
    vertices_working = np.dot(R, vertices_working) 

    # finally place the vertices 
    for component in range(3):
        vertices_working[component,:] += centroid[component]

    vertices_into_vtk(orig_file, new_file, vertices_working)


def read_vtk(filename):
    ''' reads vtk and returns 3*n_vertices numpy array '''

    if filename.endswith('vtk'):
        strings = filename.split('.')
        print 'found vtk file ', filename 
    else:
        raise OSError('File with correct extension not found')

    vtk_file  = open(filename, 'r')
    vtk_file_as_string = vtk_file.read()
    vtk_file.close()

    # take everything after float or double 
    vertices_split_1 = re.split('float|double', vtk_file_as_string)[1]
    # and stops at METADATA of CELLS
    all_vertices_string = re.split('METADATA|CELLS', vertices_split_1)[0]

    # may or may not be leading white space, remove if so 
    all_vertices_string = all_vertices_string.lstrip()

    vertices = []
    for vertex_str in all_vertices_string.split():
        vertices.append(float(vertex_str))

    if (len(vertices) % 3) is not 0: 
        raise ValueError('must have three dimensional data')
    n_vertices = len(vertices) / 3

    vertices = np.array(vertices)

    # reshape continiguous vertices in rows 
    # so without transpose this is each vector in the row
    # and with transpose, vectors are in columns as desired 
    vertices = vertices.reshape((n_vertices,3)).T

    return vertices


def vertices_into_vtk(vtk_file, vtk_file_new, vertices):
    ''' 
    takes an example vtk file and writes vertices into vtk_file_new
    topology is kept constant
    '''

    with open(vtk_file, 'r') as vtk_file_obj: 

        vtk_file_as_string = vtk_file_obj.read()

        # header ends with float\n
        # header into vertices_split_1[0]
        # rest of file into vertices_split_1[1]
        vertices_split_1 = re.split('float|double', vtk_file_as_string, 1)
        header = vertices_split_1[0] + ' float\n'

        # next thing after vertices is 'METADATA'
        if 'METADATA' in vertices_split_1[1]:
            vertices_split_2 = (vertices_split_1[1]).split('METADATA', 1)
            end_prefix = 'METADATA'
        elif 'CELLS' in vertices_split_1[1]:
            vertices_split_2 = (vertices_split_1[1]).split('CELLS', 1)
            end_prefix = 'CELLS'
        else:
            raise ValueError("Vertices must end with METADATA or CELLS and otherwise parsing not implemented")

        end_string = end_prefix + vertices_split_2[1]

        f = open(vtk_file_new, 'w')
        f.write(header)

        if vertices.shape[0] is not 3: 
            raise ValueError('Must use three dimensional data')
        n_vertices = vertices.shape[1]

        for idx in range(n_vertices):
            for component in range(3):
                f.write(str(vertices[component,idx]) + ' ')
            f.write('\n')

        f.write(end_string)
        f.close()
    
        print "wrote vtk file ", vtk_file

   

if __name__== "__main__":


    parser = argparse.ArgumentParser(description = """
    Positions and places an aortic valve inside an aorta geometry    
    Developed by Alexander D. Kaiser, PhD

    Input files: 
    
    Original aortic valve: av_converted.vtk as an ASCII vtk file 
    Ring for the aortic root: edgefile.vtk as an ASCII vtk file 


    Output files: 
    aortic_valve_positioned_2.vtk as an ASCII vtk file for further processing in meshmixer""", usage='%(prog)s [OPTIONS]', formatter_class=argparse.RawTextHelpFormatter)

    args = parser.parse_args()


    orig_file = 'av_converted.vtk'
    new_file  = 'aortic_valve_positioned_2.vtk'
    
    # assumed to be desired units already
    ring_file = 'edgefile.vtk'

    # if original file is to be scaled to match ring_file units 
    # or otherwise adjusted in size
    mm_to_cm = 0.1
    scaling_orig = .7 * mm_to_cm   

    R_0 = rotation_matrix_x(np.pi/2.0)

    #Add rotation to fit shit in sinus
    R_0 = np.dot(rotation_matrix_z((np.pi/180)*15), R_0)

    #print R_0
    #print np.dot(R_0,R_0.T)
    position_vertices(orig_file, new_file, ring_file, scaling_orig, R_0)

