import numpy as np
import sys
# Read from a file.
with open("problem8.txt", ) as file:
    data = file.read()  # This reads file as string.

file = open('output.txt', 'w')
sys.stdout = file

# Print data for validation.
print(">>>>>>>>>> input <<<<<<<<<<<")
print(data)

# Split lines to get a list of lines.
data = data.splitlines()
# print(data, type(data))

# Split every character in a single line.
LP_Problem = []
for line in data:
    LP_Problem.append(line.split())
# print(LP_Problem, type(LP_Problem))


# Convert LP Problem from a list of strings to a list of arrays.
for i in range(len(LP_Problem)):
    LP_Problem[i] = np.asarray(LP_Problem[i]).astype(float)
# print(LP_Problem, type(LP_Problem))
# print(LP_Problem[0], type(LP_Problem[0]))

print(">>>>>>>>>> output <<<<<<<<<<")

# m constraints, n variables
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

# Create initial basic feasible solution.
N = []
BFS = []
for i in range(1, n + 1):
    N.append("x{}".format(i))  # Create non-basic(decision) variables.
    BFS.append("x{}={}".format(i, float(0)))
B = []
for i in range(m):
    B.append("x{}".format(n + i + 1))  # Create basic(slack) variables.
    BFS.append("x{}={}".format(n + i + 1, b[i]))
x = N + B
N_idx = list(range(1, n + 1))
B_idx = list(range(n + 1, n + m + 1))
N_idx_orijinal = np.copy(N_idx)
B_idx_original = np.copy(B_idx)
print("x = {", *x, "}")
print("N = {", *N, "}", "\tB = {", *B, "}")
print("\nInitial BFS:", *BFS[:n], "\n\t    ", *BFS[n:])
# print(N)
# print(B)

#########################################################################################


# Initial values
Basic_coef = np.identity(m)
Basic_coef_inv = np.linalg.inv(Basic_coef)
cb = np.array([[0]] * m)
cb_transpose = np.transpose(cb)
kapa = np.dot(cb_transpose, Basic_coef_inv)  # Here kapa has no meaning. It is just a dummy variable.
M_inverse = np.identity(m + 1)
initial_resulting_matrix = np.dot(M_inverse, np.concatenate(([0], b), 0))

z = np.array([0] * len(c))
y = [0] * m
a = 1
E = []
xb = b
entering_var_candidates = [-1]
# Print option to get rid of scientific notation.
# any(i<0 for i in entering_var_canditates)
np.set_printoptions(precision=4, suppress=True)

while 1:

    print("#################################################################################")
    print("{}.loop".format(a))
    print("y =", *y)
    I_m = np.identity(m)
    if a == 1:
        print("cbar\t", end=" ")
        for i in range(n):
            print("x{}:".format(N_idx[i]), c[i], "\t", end=" ")
        print()

    # Define entering variable.
    entering_var_candidates = []
    carpan = np.concatenate(([1], kapa[0]))
    # print("\nKapa:", carpan)

    for i in range(len(c)):
        zj_cj = np.transpose(np.concatenate(([[-c[i]]], [A[:, i]]), axis=1))
        entering_var_candidates.append(np.dot(carpan, zj_cj))
    entering_var_candidates = np.asarray(entering_var_candidates)
    # print("z-c", np.round(entering_var_candidates), type(entering_var_candidates))
    # print("ennnnn", entering_var_candidates)

    if a > 1:
        print("cbar\t", end=" ")
        for i in range(n):
            print("x{}:{}".format(N_idx[i], -entering_var_candidates[N_idx[i] - 1]), "\t", end=" ")
        print()
    """
    entering_var_candidates = np.around(entering_var_candidates, decimals=4)
    entering_var_candidates[entering_var_candidates==0] = 0
    """
    ##############################################################################################
    # For optimal solution:
    if all(i >= 0 for i in np.around(entering_var_candidates, decimals=4)):

        print()
        print("Congrats! Optimal Solution has been achieved. \nNumber of loops:{}".format(a))
        print("Optimal value of", current_sol[0], "has been reached.")
        print("Original:", end=" ")
        B_dict = {}
        original_variables = [0] * n
        slack_variables = [0] * m
        for i in range(len(B_idx)):
            B_dict[B_idx[i]] = (bbar[i])
            # np.round(bbar[i])
        for var in N_idx_orijinal:
            if var in B_dict:
                original_variables[var - 1] = (B_dict[var])

        for var in B_idx_original:
            # print(var, m, B_idx_original, B_dict, slack_variables)
            if var in B_dict:
                slack_variables[var - B_idx_original[0]] = B_dict[var]
        for i in range(n):
            print("x{} = {:.2f}".format((i + 1), original_variables[i]), end=" ")
        print()
        print("Slack:", end=" ")
        for i in range(m):
            print("x{} = {:.2f}".format((i+n+1), slack_variables[i]), end=" ")
        print()
        # print(Basic_coef_inv)
        # print(current_sol)
        break
    ##############################################################################################

    entering_var = min(entering_var_candidates)
    entering_var_idx = np.argmin(entering_var_candidates)
    p = entering_var_idx
    # print("p=", p)
    # print("\nMinimum Value:", *entering_var, "\nEntering var:", "x{}".format(entering_var_idx + 1))
    print("Entering variable is x{}".format(entering_var_idx + 1))

    # Extract coefficients of entering variable.
    # print("Basic_coef_inv:",Basic_coef_inv)
    t = np.dot(Basic_coef_inv, A[:, p])
    # print("Matrix A:",A)
    # print("Coefficients of entering variable in a matrix A:", *t)
    print("abarj ={}".format(t))

    # Check if the solution is unbounded.Condition: if every elements of t lower than zero,
    #  solution is unbounded.
    """    
    t = np.around(t, decimals=4)  # To solve decimal point issues and negative zeros.
    t[t == 0] = 0
    """
    ########################################################################################

    # For unbounded solution:
    if all(i <= 0 for i in np.round(t, decimals=4)):
        print("This solution is unbounded. Number of loops:{}".format(a))
        # print(t)
        # print(B_idx, B_idx_original)
        # print(current_sol)
        print("Original:", end=" ")
        B_dict = {}
        original_variables = [0] * n
        slack_variables = [0] * m
        for i in range(len(B_idx)):
            B_dict[B_idx[i]] = np.float(bbar[i])

        for var in N_idx_orijinal:
            if var in B_dict:
                original_variables[var - 1] = (B_dict[var])

        for var in B_idx_original:
            if var in B_dict:
                slack_variables[var - B_idx_original[0]] = B_dict[var]
            # print(var, m, B_idx_original, B_dict, slack_variables)
        for i in range(n):
            print("x{} = {:.2f}".format((i + 1), original_variables[i]), end=" ")
        print()
        print("Slack:", end=" ")
        for i in range(m):
            print("x{} = {:.2f}".format((i + n + 1), slack_variables[i]), end=" ")
        print()

        print("z =", current_sol[0], end=" + ")
        for i in range(len(N_idx)-1):
            print("({:.2f})x{}".format(*-entering_var_candidates[N_idx[i]-1], N_idx[i]), end=" + ")
        print("({:.2f})x{}".format(*-entering_var_candidates[N_idx[-1]-1], N_idx[-1]))
        #print(M_inverse[1:, 1:])
        sol = M_inverse[1:, 1:]
        """for i in range(len(B_idx)):
            print("x{} = ".format(B_idx[i]), end=" ")
            print("{:.2f}".format(original_variables[i]), end=" ")
            #print("")
            #print(eta_matrix)
        """
        print()
        # print(Basic_coef_inv)
        break

    ########################################################################################

    # Take the positive coefficients.
    positive_values = t[t > 0]
    # print("Positive values:", *positive_values)

    idx_of_positive_values = np.array(np.where(t > 0)).reshape((-1, 1))
    # print("Index of positive values:", idx_of_positive_values)

    # Calculate ratio which defines the leaving variable.
    # print("xb:",xb)
    ratio = xb[idx_of_positive_values]
    ratio = ratio / t[idx_of_positive_values]
    print("Ratio", end=" ")
    for i in range(len(ratio)):
        print("x{}:".format(B_idx[idx_of_positive_values[i][0]]), "{:3.4f}\t".format(*ratio[i]), end="")
    print()

    # Define leaving variable.
    leaving_var_idx = idx_of_positive_values[np.argmin(ratio)]
    # print(leaving_var_idx)
    leaving_var = (n + 1) + leaving_var_idx
    print("Leaving variable is x{}".format(*leaving_var))
    q = leaving_var_idx
    # print("q =", *q)

    eta_vector = -t
    # print(eta_vector)
    eta_vector[q] = 1
    # print(eta_vector)
    # print("t:", *t)
    print("E{} = column {}:{}".format(a, *leaving_var - n, t))
    eta_vector = eta_vector / t[q]
    # print("eta vector:", eta_vector)

    I_m[:, q] = eta_vector.reshape((-1, 1))
    eta_matrix = I_m
    # print("eta matrix:", eta_matrix)

    Basic_coef_inv = np.dot(eta_matrix, Basic_coef_inv)
    # print("Basic coefficient inverse:", Basic_coef_inv)

    cb[q] = c[p]
    cb_transpose = np.transpose(cb)
    # print("cb:", cb)

    kapa = np.dot(cb_transpose, Basic_coef_inv)
    # print("cb_transpose*B_inverse:", kapa)

    M_inverse = np.concatenate(([[1]], kapa), axis=1)
    dim = len(M_inverse[0]) - 1
    dim = np.array([[0]] * dim)
    expanded_B_inv = np.concatenate((dim, Basic_coef_inv), axis=1)
    M_inverse = np.concatenate((M_inverse, expanded_B_inv))
    # print("M_inverse:", M_inverse)

    # The current solution is:
    current_sol = np.dot(M_inverse, np.concatenate(([0], b)))
    # print(current_sol)
    xb = np.dot(eta_matrix, xb)
    # print("Condition for continuing:", entering_var_candidates)
    """
    # Organize basic and non-basic variables.
    N_copy = np.copy(N)
    B_copy = np.copy(B)
    N[int(p)] = B_copy[int(q)]
    B[int(q)] = N_copy[int(p)]
    print("N = {", *N, "}", "\tB = {", *B, "}")
    """
    # Organize basic and non-basic variables.
    N_copy = np.copy(N)
    B_copy = np.copy(B)
    N[int(p)] = B_copy[int(q)]
    B[int(q)] = N_copy[int(p)]
    print("N = {", *N, "}", "\tB = {", *B, "}")
    # To calculate cbar we have to use index of variables.
    N_idx_copy = np.copy(N_idx)
    B_idx_copy = np.copy(B_idx)
    N_idx[int(p)] = B_idx_copy[int(q)]
    B_idx[int(q)] = N_idx_copy[int(p)]
    # print(N_idx, B_idx)
    # Print current solution
    # print(current_sol, B)

    bbar = current_sol[1:]
    print("bbar =", bbar)

    y = kapa
    a += 1

file.close()
