# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 11:33:51 2020

@author: giova
"""
import numpy as np

class BoundaryConditions:

    def __init__(self):
        self.disp = None
        self.neumann_elements = None
        self.dirichlet_elements = None
        self.dirichlet_nodes = None
        self.neumann_nodes = None
        self.load = None #array of distributed loads

    def find_dofs(self,mesh,nodes):
        dn = mesh.dofspernode
        dofs = np.array([nodes * dn,
                          nodes * dn + 1])
        dofs = np.concatenate(dofs.T)
        return dofs

    def find_boundary_obj(self, mesh, tag):
        boundary_elements = mesh.cell_sets_dict[tag]['line'] # Array of element numbers with tag "tag"
        boundary_nodes = mesh.cells_dict['line'][boundary_elements] # table of nodes of elements with tag "tag"
        # boundary_nodes = np.unique(boundary_nodes) # only consider nodes once
        return boundary_elements , boundary_nodes

    def apply_bcs(self, F, K, mesh):

        # Method for distributed load on the boundary of a 2D element (only constant loads on the xy plane for now)
        for i in self.neumann_elements:

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
        K[:,dirichlet_dofs] = 0
        K[dirichlet_dofs,:] = 0
        K[dirichlet_dofs,dirichlet_dofs] = 1