import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


def create_u_vector(L, num_points_lattice, spatial_step):
    # Create array with the given number of points
    u_vector = np.zeros(num_points_lattice + 1) # Add 1 to ensure symmetry
    # Each element in u should be separated by the spatial stepsize
    for i in range(len(u_vector)):
        u_vector[i] = i * spatial_step
        # Subtract L so that the first value of u is -L and the final value is L
        u_vector = u_vector - L
    return u_vector


def initial_conditions(a, A, u, v):
    # Initial cancer cells
    nC_0 = (10 ** 4) * (1 + a * np.cos(A * u))
    # Initial T-cells
    nT_0 = (10 ** 4) * (2 + a * np.cos(A * v))
    return nC_0, nT_0
def g(x_i, y, chi):
    # Find difference between point x_i in lattice and each element y_i of attice
    d = abs(x_i - y)
    # We look for the points y_i close to x_i, so for the elements of d smaller than chi
    close_y_index = np.where(d <= chi)[0]
    # Find size of interval containing elements of lattice within chi of x_i
    l_max = min(x_i + chi, y[np.max(close_y_index)])
    l_min = max(x_i - chi, y[np.min(close_y_index)])
    47
    l = l_max-l_min
    # Create vector of values of g at each y_i
    g_vals = np.zeros(len(y))
    # If l>0 we define the value of g for points withing chi of x_i as 1/l
    if l != 0:
        for index in close_y_index:
            g_vals[index] = 1 / l
    return g_vals
def calculate_competition_vector(x, y, selectivity_range, n_vec, spatial_step):
    # Create vector with zeros of length y
    competition_vector = np.zeros(len(y))
    # Iterate through each element x_i in x
    for i in range(len(x)):
    # Apply the function g given the selectivity range
    g_vec = g(x[i], y, selectivity_range)
    # Multiply each element in g_vec by the corresponding element in n_vec
    and the spatial step, then sum over all elements
    competition_vector[i] = np.sum(spatial_step * np.multiply(g_vec, n_vec))
    return competition_vector


def calculate_gamma_times_j(gamma_matrix, x, y, selectivity_range, n_vec, spatial_step):
    # Create vector with zeros of length y
    vector = np.zeros(len(y))
    # Iterate through each element x_i in x
    for i in range(len(x)):
    # Apply the function g given the selectivity range
    g_vec = g(x[i], y, selectivity_range)
    # Multiply each element in g_vec by the corresponding element in n_vec
    and the spatial step, then sum over all elements
    vector[i] = np.sum(spatial_step *
    np.multiply.reduce([g_vec, n_vec, gamma_matrix[i,:]]))
    return vector


def calculate_Rc_interaction_matrix(alpha, mu, K, zeta, gamma_times_j):
    R = alpha - mu * K - zeta * gamma_times_j
    return R
def calculate_Rt_interaction_matrix(alpha, mu, K, zeta, gamma_times_j):
    R = alpha - mu * K + zeta * gamma_times_j
    return R
def create_matrix(nC_k):
    # Obtain length of nC_k
    n = len(nC_k)
    48
    # Create central diagonal array
    diagonal = -2 * np.ones(n)
    # Set first and last element to 0
    diagonal[0] = 0
    diagonal[-1] = 0
    # Create upper diagonal array
    upper_diagonal = 1 * np.ones(n - 1)
    # Set first element to 0
    upper_diagonal[0] = 0
    # Create lower diagonal array
    lower_diagonal = 1 * np.ones(n - 1)
    # Set last element to 0
    lower_diagonal[-1] = 0
    # Create Sparse matrix
    diagonals = [lower_diagonal, diagonal, upper_diagonal]
    offset = [-1, 0, 1]
    A = sp.sparse.diags(diagonals, offset)
    return A
def update_parameters(nC_k, nT_k, time_step, spatial_step, R_c, R_t, beta_c,
matrix, lambda_c):
    # Cancer cells
    nC_k_plus_one_half = np.multiply(nC_k, np.divide(1 + time_step * np.maximum(0,
    R_c),
    1 + time_step * abs(np.minimum(0, R_c))))
    nC_k_plus_one = np.zeros(len(nC_k))
    nC_k_plus_one[1:-1] = nC_k_plus_one_half[1:-1] + beta_c * (time_step / spatial_step
    ** 2) * matrix.dot(nC_k_plus_one_half)[1:-1]
    if lambda_c > 0:
    # Neumann Boundary Conditions
    nC_k_plus_one[0] = nC_k_plus_one[1]
    nC_k_plus_one[-1] = nC_k_plus_one[-2]
    # T-cells
    nT_k_plus_one = np.multiply(nT_k, np.divide(1 + time_step * np.maximum(0,
    R_t),
    1 + time_step * abs(np.minimum(0, R_t))))
    return nC_k_plus_one, nT_k_plus_one
def plot_total_cells_evolution(time_vector, total_cancer_cells, total_t_cells):
    # Plot results of total cells over time
    plt.title(’Evolution of Cancer cells and T-cells over time’, fontsize=20)
    plt.plot(time_vector[:-1], total_cancer_cells, label=’Cancer cells’)
    plt.plot(time_vector[:-1], total_t_cells, label=’T-cells’)
    plt.xlabel(’Time’, fontsize=18)
    49
    plt.ylabel(’Total number of cells’, fontsize=18)
    plt.legend(fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.show()
def interaction_matrix_model(L, num_points_lattice, time_step, final_time, a, A, theta_c, theta_t, eta, alpha_c, mu_c, zeta_c, alpha_t, mu_t, zeta_t, gamma_matrix, lambda_c, plot_cell_evolution=True, print_results=True):
    # Calculate spatial step
    spatial_step = 2 * L / num_points_lattice
    # Calculate beta_c
    beta_c = lambda_c * (spatial_step**2) / (2 * time_step)
    # Create empty lists to store cell density at each point over time
    nC_matrix = []
    nT_matrix = []
    # Define lattice points
    u_vector = create_u_vector(L, num_points_lattice, spatial_step)
    v_vector = u_vector.copy()
    # Define initial conditions
    nC_k, nT_k = initial_conditions(a, A, u_vector, v_vector)
    nC_matrix.append(nC_k)
    nT_matrix.append(nT_k)
    # Define matrix for spatial discretisation
    matrix = create_matrix(nC_k)
    # Define total cell lists for each time
    total_cancer_cells = []
    total_t_cells = []
    time = 0
    time_vector = [0]
    while time <= final_time:
    if print_results:
    print(’time =’, round(time, 2))
    # Count cells of each type across lattice
    total_cancer_cells.append(np.sum(nC_k) * spatial_step)
    total_t_cells.append(np.sum(nT_k) * spatial_step)
    # Calculate Kc, Kt, Jc, Jt with the corresponding vectors nC_k and nT_k
    Kc = calculate_competition_vector(u_vector, u_vector, theta_c, nC_k,
    spatial_step)
    Kt = calculate_competition_vector(v_vector, v_vector, theta_t, nT_k,
    spatial_step)
    gamma_Jc = calculate_gamma_times_j(gamma_matrix, u_vector, v_vector,
    50
    eta, nT_k, spatial_step)
    gamma_Jt = calculate_gamma_times_j(gamma_matrix.T, v_vector, u_vector,
    eta, nC_k, spatial_step)
    # Use Kc, Kt, Jc, Jt and the matrix gamma to calculate R_c and R_t
    R_c = calculate_Rc_interaction_matrix(alpha_c, mu_c, Kc, zeta_c, gamma_Jc)
    R_t = calculate_Rt_interaction_matrix(alpha_t, mu_t, Kt, zeta_t, gamma_Jt)
    # Update vectors
    nC_k, nT_k = update_parameters(nC_k, nT_k, time_step, spatial_step, R_c,
    R_t, beta_c, matrix, lambda_c)
    nC_matrix.append(nC_k)
    nT_matrix.append(nT_k)
    # Update time
    time += time_step
    time_vector.append(time)
    if print_results:
    # Print results at each timestep
    print(’Total cancer cells:’, total_cancer_cells[-1])
    print(’Total T-cells:’, total_t_cells[-1])
    # Plot results of total cells over time
    if plot_cell_evolution:
    plot_total_cells_evolution(time_vector, total_cancer_cells, total_t_cells)
    return nC_matrix, nT_matrix, u_vector, time_vector