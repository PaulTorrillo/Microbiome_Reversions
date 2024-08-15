
from sympy import symbols, Matrix, diag, eye, N, simplify, lambdify, exp
import numpy as np
import matplotlib.pyplot as plt
# Note p is stand in for 1/tau
# Define the variables
r, f, n, p= symbols('r f n p', real=True, positive=True)
t = symbols('t', integer=True, nonnegative=True)  # t is discrete and non-negative

# Define the transition matrix for a Markov chain
A = Matrix([
    [1-p, 0, 0],
    [p, 1-p, p],
    [0, p, 1-p]
])

# Define the transition matrices for a Markov chain in continuous time
A2 = Matrix([
    [-p, 0, 0],
    [p, -p, p],
    [0, p, -p]
])
# Perform eigendecomposition
P, D = A.diagonalize()

# Formulate the state transition matrix for discrete time using eigendecomposition
state_transition_matrix_discrete = P * D**t * P.inv()


# Compute the matrix exponential for continuous time state transitions
state_transition_matrix_continuous = exp(A2 * t)
# Print the algebraic formula for the state transition matrix
print("Algebraic formula for the state transition matrix in continuous time:")
print(state_transition_matrix_continuous)
# Simplify the expression for better readability
state_transition_matrix_discrete_simplified = simplify(state_transition_matrix_discrete)

# Print the algebraic formula for the state transition matrix
print("Algebraic formula for the state transition matrix in discrete time:")
print(state_transition_matrix_discrete_simplified)

# Define the initial vector

A=state_transition_matrix_continuous

part1=simplify(n*(p/(p+p))*(A[1,1]*A[2,1]+A[2,1]*A[1,1]+2*A[2,1]*A[2,1]))
part2=simplify(n*(p/(p+p)))*(A[0,0]*A[1,0]+A[1,0]*A[0,0]+2*A[0,0]*A[2,0]+2*A[2,0]*A[0,0]+A[1,0]*A[2,0]+A[2,0]*A[1,0]+2*A[2,0]*A[2,0])

initial_equilibrium=part1+part2
initial_WT=simplify(n*(A[1,1]*A[2,1]+A[2,1]*A[1,1]+2*A[2,1]*A[2,1]))
print('dN assuming equilibrum')
print(simplify(initial_equilibrium))
print('dN assuming no initial transient')
print(simplify(initial_WT))
