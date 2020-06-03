# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 11:33:51 2020

@author: giova
"""
import numpy as np
import scipy.sparse as sps

class BoundaryConditions:

    def __init__(self):
        self.disp = None
        self.neumann_elements = None
        self.neumann_nodes = None
        self.dirichlet_elements = None
        self.dirichlet_nodes = None
        self.load = None #array of distributed loads

    def find_dofs(self,mesh,nodes):
        dn = mesh.dofspernode
        dofs = np.array([nodes * dn,
                          nodes * dn + 1])
        dofs = np.concatenate(dofs.T)
        return dofs

    def find_boundary_obj(self, mesh, tag):
        if mesh.d == 1:
            boundary_elements = mesh.cell_sets_dict[tag]['vertex'] # Array of element numbers with tag "tag"
            boundary_nodes = mesh.cells_dict['vertex'][boundary_elements] # table of nodes of elements with tag "tag"
        elif mesh.d == 2:
            boundary_elements = mesh.cell_sets_dict[tag]['line'] # Array of element numbers with tag "tag"
            boundary_nodes = mesh.cells_dict['line'][boundary_elements] # table of nodes of elements with tag "tag"
        elif mesh.d == 3:
            boundary_elements = np.hstack((mesh.cell_sets_dict[tag]['triangle'], mesh.cell_sets_dict[tag]['quad'])) # Array of element numbers with tag "tag"
            boundary_nodes = np.hstack((mesh.cells_dict['triangle'][boundary_elements], mesh.cells_dict['quad'][boundary_elements])) # table of nodes of elements with tag "tag"
        # boundary_nodes = np.unique(boundary_nodes) # only consider nodes once
        return boundary_elements , boundary_nodes

    def apply_bcs(self, F, K, mesh):
        
        
        #line_element = mesh.cell_sets_dict['Neumann']['line']
        
        # Method for distributed load on the boundary of a 2D element (only constant loads on the xy plane for now)
        for i in range(self.neumann_elements.size):

            nodes = self.neumann_nodes[i]
            dofs = self.find_dofs(mesh , nodes)

            cds = mesh.points[nodes]
            length = np.sqrt((cds[1][0]-cds[0][0])**2+(cds[1][1]-cds[0][1])**2)

            c = (cds[1,0]-cds[0,0])/length
            s = (cds[1,1]-cds[0,1])/length

            R = np.array([[c , s],
                          [-s , c]])

            load_nat = R @ self.load.T
            load_nat = load_nat.reshape((1,2))

            N = np.array([.5 , .5]).reshape((1,2))
            f_nat = length/2 * 2 * (N.T @ load_nat) # det(j) * w_i * f(N @ q)
            f = f_nat @ R
            f = np.concatenate(f)
            F[dofs] += f

        # Find dofs where Dirichlet conditions are enforced
        dirichlet_dofs = self.find_dofs(mesh,np.unique(self.dirichlet_nodes))
        F[dirichlet_dofs] = 0 # zeroing forces on nodes with zero displacement

        # zeroing out zero displacement columns and rows
        K[:,dirichlet_dofs] = 0 # row slicing
        K[dirichlet_dofs,:] = 0 # column slicing
        K[dirichlet_dofs,dirichlet_dofs] = 1
