import numpy as np
import os

# Read from a file.
with open("deneme.txt", ) as file:
    data = file.read()  # This reads file as string.

# Print data for validation.
print(">>>>>>>>>> input <<<<<<<<<<<")
print(data)

# Split lines to get a list of lines.
data = data.splitlines()
# print(data,type(data))

# Split every character in a single line.
LP_Problem = []
for line in data:
    LP_Problem.append(line.split())
# print(LP_Problem,type(LP_Problem))

# Convert LP Problem from a list of strings to a list of arrays.
for i in range(len(LP_Problem)):
    LP_Problem[i] = np.asarray(LP_Problem[i]).astype(float)
# print(LP_Problem,type(LP_Problem))
print(">>>>>>>>>> output <<<<<<<<<<")

# m constaints,n variables
m = int(LP_Problem[0][0])
n = int(LP_Problem[0][1])
print("Constraints m =", m, "\tVariables n =", n)
LP_Problem.pop(0)  # drop these values to make it easier.

# c = objective function
c = LP_Problem[0]
c = np.append(c, [0] * m)
print("c = ", *c)
LP_Problem.pop(0)  # Drop objective function.

LP_Problem = np.asarray(LP_Problem)  # Convert LP prob. into array.
b = LP_Problem[:, -1]  # Take the last column.
print("b = ", *b)

A = np.delete(LP_Problem, -1, axis=1)  # Delete last column to get coefficient matrix.
A = np.concatenate((A, np.identity(m)), axis=1)  # Identity matrix joins coefficient matrix.
print("A =\n", A)
##########################################################################################
#Create initial basic feasible solution.
N = []
BFS = []
for i in range(1, n + 1):
    N.append("x{}".format(i)) #Create non-basic(decision) variables.
    BFS.append("x{}={}".format(i, float(0)))
B = []
for i in range(m):
    B.append("x{}".format(n+i+1)) #Create basic(slack) variables.
    BFS.append("x{}={}".format(n+i+1, b[i]))
x = N + B
print("x = {", *x, "}")
print("N = {", *N, "}","\tB = {", *B, "}")
print("\nInitial BFS:", *BFS[:n],"\n\t    ",*BFS[n:])
#########################################################################################

Basic_coef = np.identity(m)
Basic_coef_inv = np.linalg.inv(Basic_coef)
cb = np.array([[0]]*m)
cb_transpose = np.transpose(cb)

kapa = np.dot(cb_transpose,Basic_coef_inv)
M_inverse = np.identity(m+1)
initial_resulting_matrix = np.dot(M_inverse,np.concatenate(([0],b),0))
z = np.array([0]*len(c))
y = [0]*m
a=0
E=[]

xb = b
while(a<2):

    I_m = np.identity(m)
    print("y = ", *y)
    for i in range(n):
        print("x{}:".format(i + 1), c[i], "\t", end=" ")
    # Define entering variable.
    entering_var_canditates = []
    carpan = np.concatenate(([1],kapa[0]))
    print("\n",carpan)

    for i in range(len(c)):
        zj_cj = np.transpose(np.concatenate(([[-c[i]]], [A[:, i]]), axis=1))
        entering_var_canditates.append(np.dot(carpan, zj_cj))
    print(entering_var_canditates)
    entering_var = min(entering_var_canditates)
    entering_var_idx = np.argmin(entering_var_canditates)
    p = entering_var_idx
    print("p=", p+1)
    print("\nMinimum Value:", *entering_var, "\nEntering var:", "x{}".format(entering_var_idx + 1))
    # Extract coefficients of entering variable.
    print("Basic_coef_inv:",Basic_coef_inv)
    t = np.dot(Basic_coef_inv, A[:, p])
    print(A)
    print("Coefficients of entering variable in a matrix A:", *t)
    # Take the positive coefficients.
    positive_values = t[t>0]
    print("Positive values:",*positive_values)
    idx_of_positive_values = np.array(np.where(t > 0)).reshape((-1, 1))
    print("Index of positive values:", idx_of_positive_values)
    # Calculate ratio which defines the leaving variable.

    print(xb)
    ratio = xb[idx_of_positive_values]
    ratio = ratio / t[idx_of_positive_values]
    print("Ratio", end=" ")
    for i in range(len(ratio)):
        print("x{}:".format(*(n+1 + idx_of_positive_values[i])), "{:3.4f}\t".format(*ratio[i]), end=" ")

    # Define leaving variable.
    leaving_var_idx = idx_of_positive_values[np.argmin(ratio)]
    leaving_var = (n+1) + leaving_var_idx
    print("\nLeaving variable: x{}".format(*leaving_var))
    q = leaving_var_idx
    print("q =", *q+1)
    eta_vector = -t
    eta_vector[q] = 1
    print("t:", t)
    eta_vector= eta_vector/t[q]
    print("eta vector:", eta_vector)

    I_m[:, q] = eta_vector.reshape((-1, 1))
    eta_matrix = I_m
    print("eta matrix:", eta_matrix)

    Basic_coef_inv = np.dot(eta_matrix, Basic_coef_inv)
    print("Basic coefficient inverse:", Basic_coef_inv)
    cb[q] = c[p]
    cb_transpose = np.transpose(cb)
    print("cb:", cb)
    kapa = np.dot(cb_transpose, Basic_coef_inv)
    print("cb_transpose*B_inverse:", kapa)
    M_inverse = np.concatenate(([[1]], kapa), axis=1)
    dim = len(M_inverse[0])-1
    dim = np.array([[0]]*dim)
    expanded_B_inv = np.concatenate((dim, Basic_coef_inv) , axis=1)
    M_inverse = np.concatenate((M_inverse, expanded_B_inv))
    print("M_inverse:", M_inverse)

    #The current solution is:
    current_sol = np.dot(M_inverse, np.concatenate(([0], b)))
    print(current_sol)
    xb = np.dot(eta_matrix,xb)






    a+=1



