import numpy as np
from sympy.physics.quantum import TensorProduct
from sympy.physics.paulialgebra import Pauli, evaluate_pauli_product
from binsymbols import *
from sympy import *
import itertools
import math

def alg3(B, T):
    """
    Parameters
    B : np.array(str), feasible states
    T : np.matrix, transition matrix
    returns -> pauli_string : str, hamilton in pauli basis
    """
    pauli_string = 0
    J = len(B)
    for j in range(J):
        for k in range(j + 1, J):
            if not math.isclose(T[j, k], 0, abs_tol = 1e-7):
                w, z = B[j], B[k]
                
                # compute A1
                if z[0] == '0' and w[0] == '0':
                    A_n = I + Z
                elif (z[0] == '0' and w[0] == '1') or (z[0] == '1' and w[0] == '0'):
                    A_n = X
                else: # z[0] = 1 and w[0] = 1
                    A_n = I - Z

                # compute B1
                if z[0] == '0' and w[0] == '1':
                    B_n = -Y
                elif z[0] == '1' and w[0] == '0':
                    B_n = Y
                else: # (z[0] == 0 and w[0] == 0) or (z[0] == 1 and w[0] == 1)
                    B_n = 0.

                # compute An, Bn
                for n in range(1, len(z)):
                    x, y = z[n], w[n]

                    if x == '0' and y == '0':
                        A_next = 0.5 * TensorProduct(A_n, I + Z)
                    elif x == '0' and y == '1':
                        A_next = 0.5 * (TensorProduct(A_n, X) + TensorProduct(B_n, Y))
                    elif x == '1' and y == '0':
                        A_next = 0.5 * (TensorProduct(A_n, X) - TensorProduct(B_n, Y))
                    else: # x == 1 and y == 1
                        A_next = 0.5 * TensorProduct(A_n, I - Z)

                    if x == '0' and y == '0':
                        B_next = 0.5 * (TensorProduct(B_n, I + Z))
                    elif x == '0' and y == '1':
                        B_next = 0.5 * (TensorProduct(B_n, X) - TensorProduct(A_n, Y))
                    elif x == '1' and y == '0':
                        B_next = 0.5 * (TensorProduct(B_n, X) + TensorProduct(A_n, Y))
                    else: # x == 1 and y == 1
                        B_next = 0.5 * (TensorProduct(B_n, I - Z))

                    A_n = A_next
                    B_n = B_next
                    
                pauli_string += A_n * T[j, k]
    return pauli_string
