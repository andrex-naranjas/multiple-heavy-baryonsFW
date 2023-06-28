#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------
 authors: C. A. Vaquera Araujo (vaquera@fisica.ugto.mx)
          A. Rivero (ailierrivero@gmail.com)
          A. Ramirez Morales (andres@knu.ac.kr)
------------------------------------------------------------
"""
import sympy.physics.mechanics as mech
from sympy import *
from sympy.physics.quantum import TensorProduct
import numpy as np


class FlavorAmplitudes():
    """
    Class to get electromagnetic flavor amplitudes
    """        
    def __init__(self, baryons, m_q=299, m_s=465, m_c=1606):
        # baryons available: omegas, sigmas, lambdas, cascades_prime, cascades
        self.m_baryons = baryons 
        # Values of the magnetic moments of each quark according to Dothan 1982
        self.mu_u = (2/3) * 1/(2 * m_q)
        self.mu_d = (-1/3) * 1/(2 * m_q)
        self.mu_s = (-1/3) * 1/(2 * m_s)
        self.mu_c = (2/3) * 1/(2 * m_c)
        
        # up quark state
        self.u = Matrix([[1],[0],[0],[0]]) 
        # down quark state
        self.d = Matrix([[0],[1],[0],[0]]) 
        # strange quark state
        self.s = Matrix([[0],[0],[1],[0]]) 
        # charm quark state
        self.c = Matrix([[0],[0],[0],[1]])

    def magnetic_operators(self, index): 
        """
        Method to calculate the Magnetic Operators (mu(i))
        index == operator index (1,2,3)
        """
        return Matrix([[self.mu_u, 0, 0, 0], [0, self.mu_d, 0, 0], [0, 0, self.mu_s, 0], [0, 0, 0, self.mu_c]])

    def identity_matrix(self, dim):
        """
        Method to obtain the identity matrix
        """
        return eye(dim)

    def tensor_product_operator_1(self, index=1, dim=4): 
        """
        Method to define the direct product 2⊗2⊗2 magnetic operator for mu_1 (Mu(1))
        """
        return TensorProduct(TensorProduct(self.magnetic_operators(index), self.identity_matrix(dim)), self.identity_matrix(dim))

    def tensor_product_operator_2(self, index=2, dim=4):
        """
        Method to define the direct product 2⊗2⊗2 magnetic operator for mu_2 (Mu(2))
        """
        return TensorProduct(TensorProduct(self.identity_matrix(dim), self.magnetic_operators(index)), self.identity_matrix(dim))

    def tensor_product_operator_3(self, index=3, dim=4):
        """
        Method to define the direct product 2⊗2⊗2 magnetic operator for mu_3 (Mu(3))
        """
        return TensorProduct(TensorProduct(self.identity_matrix(dim), self.identity_matrix(dim)), self.magnetic_operators(index))        
    
    def flavor_state(self, state_a, state_b, state_c): 
        """
        Method to calculate the flavor states (Fs)
        state_a,b,c should be sympy matrices of the form: Matrix([[0],...,[1]])
        """
        return TensorProduct(TensorProduct(state_a ,state_b), state_c)
    
    #Flavor States of the triplet
    def Xi_cc_pp_p(self):
        """
        Method to define the Xi_cc_pp_p flavor state
        """
        return self.flavor_state(self.c, self.c, self.u)

    def Xi_cc_p_p(self):
        """
        Method to define the Xi_cc_p_p flavor state
        """
        return self.flavor_state(self.c, self.c, self.d)
    
    def Omega_cc_p_p(self):
        """
        Method to define the Omega_cc_p_p flavor state
        """
        return self.flavor_state(self.c, self.c, self.s)

    
    

    #Flavor State of the singlet
    def Omega_ccc_pp_p(self):
        """
        Method to define the Omega_ccc_pp_p flavor state
        """
        return (self.flavor_state(self.c, self.c, self.c))
    
    def flavor_matrix_elements(self, n, x_i, y_j):
        """
        Method to calculate the matrix element for the flavor part
        """
        if n==1:
            Mu_i=TensorProduct(TensorProduct(self.magnetic_operators(index=1), self.identity_matrix(dim=4)), self.identity_matrix(dim=4))  
        else: 
            if n==2:
                Mu_i=TensorProduct(TensorProduct(self.identity_matrix(dim=4), self.magnetic_operators(index=2)), self.identity_matrix(dim=4))
            else:
                Mu_i=TensorProduct(TensorProduct(self.identity_matrix(dim=4), self.identity_matrix(dim=4)), self.magnetic_operators(index=3))
        return np.array(conjugate(transpose(x_i)) * Mu_i * y_j).astype(np.float64).flatten()[0]
