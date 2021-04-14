"""
    Program that takes as input either a Kirchhoff matrix or
    an affinity matrix to compute hit and commute times.
    
    Affinity matrices could be defined as:
    $A_{ij} = \frac{N_{ij}}{\sqrt{N_{i}.N_{j}}$
    
    where: 
    $N_{ij}$ refers to the number of atoms from residue i 
            in contact with residue j
    $N_{i}$ refers to the number of heavy atoms in residue i
    $N_{j}$ refers to the number of heavy atoms in residue j

    Kirchhoff matrix is defined as K = D - A
    where D is a diagonal matrix with elements that are row or 
    column sums of A.

    Kirchhoff matrix will also be referred to as graph Laplacian
    in the code below.

"""

import os, string, sys
from numpy import *
#import numpy as np
from scipy import sparse, linalg

class IT_HitCommute:
    def __init__(self, Aij):
        """
        Constructor
        
        @attention: works only on Affinity matrices defined above
        
        @type: An n x n matrix with Affinity entries computed from either a single structure or an ensemble of structures
        @param: $A_{ij}$ is a matrix
        
        @rtype: None
        @return: A new object of type hitCommuteComputer
        """
        self.size = 0
        self.size = len(Aij)
        assert self.size > 0, 'Bad argument: size: %d' %(self.size)
        self.hitTime = []
        self.comTime = []
        self.laplacian = []
        for i in range(0, self.size):
            self.hitTime.append([])
            self.comTime.append([])
            
    def makeLaplacian(self, Aij):
        """
        same as making the Kirchhoff matrix, except that 
        nodes/residues are weighed based on the number of atoms.
        @rtype: A matrix of n x n, where n is the number of residues in the matrix
        """
        D = sum(Aij, 1)
        self.laplacian = spdiags(D, 0, len(D), len(D)) - A
        
        return self.laplacian
            

    def buildHitTimes(self, A):
        """
        Hit Time computer - main function to calculate the hit times of a graph
        @param Aij: Affinity matrix
        
        @rtype: A matrix of n x n where n is the number of residues in the matrix
        @return: Hit times computed from the graph
        """
        
        #D is the degree vector if A is affinity.
        #But if A is Laplace, then D is all-zero.
        D = sum(A, 1)    

    
        #in the case the input matrix is a Laplacian        
        all_zeros = not any(D)
        #print all_zeros
        if all_zeros:
            D = diag(A)
            A = sparse.spdiags(D,0,len(D),len(D)) - A;

        # Compute stationary distribution
        st = D /sum(D)
        
        # Markov transition matrix, row normalized
        P = dot(diag(D**(-1)), A)
        W = ones((len(st),1)) * st.T
        Z = linalg.pinv(eye(P.shape[0], P.shape[1]) - P + W)
        
        Hz = ones((len(st),1)) * diag(Z).T - Z 
        Hz = Hz / (ones((len(st),1)) * st.T)
        Hz = Hz.T
        
        self.hitTime = Hz
        return self.hitTime
    
    def buildCommuteTimes(self):
        """
        main function to calculate the commute times of a graph. 
        
        @attention: Commute times are just defined to be the sum of the hitTimes(i, j) + hitTimes(j, i)
        
        @param None: No parameters are expected, since the hitTimes have to be computed
        
        @rtype: A matrix of n x n where n is the number of residues in the matrix
        @return: Hit times computed from the graph
        """
        self.comTime = self.hitTime + self.hitTime.T
        return self.comTime              
    
    def driver(self, Aij):
        """
        driver routine to compute both hit times and commute times
        @param Aij:  Affinity matrix
        
        @attention: meaningful if and only if Aij is >0 and <= n, where n is the number of residues
        
        @rtype: A list of two arrays Hz, Cz
        @return: hit Times (Hz) and commute times (Cz) of the graph computed
        """
        Hz = self.buildHitTimes(Aij)
        Cz = self.buildCommuteTimes()
        return [Hz, Cz]

      
#if __name__=='__main__':    
#    a = [[5., 3., 4., 2.], [3., 4., 2., 2.], [4., 1., 4., 2.], [2., 3., 2., 5.]]
#    A = array(a)
#    print A
#    # print A.shape
#    hc = IT_hitcommute(A)
#    H = hc.buildHitTimes(A)
#    print H
#   
        
