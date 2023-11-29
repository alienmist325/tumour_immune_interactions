from marta_model import matrix_exponential_dist, interaction_matrix_model
import numpy as np
from inputs import get_sim_configuration


print("blom")

cf = get_sim_configuration()
mat_size = int(cf.no_possible_phenotypes + 1)
const_affinity_matrix = cf.binding_affinity * np.ones((mat_size, mat_size))


nC, nT, u_vector, time_vector = interaction_matrix_model(
    L=1,
    num_subintervals_lattice=int(cf.no_possible_phenotypes),
    time_step=cf.time_step,
    final_time=cf.final_time,
    A=cf.A,
    a=cf.a,
    theta_c=cf.tumour_selectivity,
    theta_t=cf.CTL_selectivity,
    eta=cf.affinity_range,
    alpha_c=cf.tumour_natural_prolif_rate,
    mu_c=cf.tumour_natural_death_rate,
    zeta_c=cf.tumour_interaction_induced_rate,
    alpha_t=cf.CTL_natural_prolif_rate,
    mu_t=cf.CTL_natural_death_rate,
    zeta_t=cf.CTL_interaction_induced_rate,
    gamma_matrix=const_affinity_matrix,
    lambda_c=cf.tumour_phenotypic_variation_probability,
    print_results=False,
)

"""
nC, nT, u_vector, time_vector = interaction_matrix_model(
    L=1,
    num_subintervals_lattice=100,
    time_step=0.05,
    final_time=100,
    a=1,
    A=5,
    theta_c=1.8,
    theta_t=1.8,
    eta=2,
    alpha_c=1.5,
    mu_c=5 * 10 ** (-6),
    zeta_c=5 * 10 ** (-6),
    alpha_t=0.05,
    mu_t=5 * 10 ** (-6),
    zeta_t=3 * 10 ** (-5),
    gamma_matrix=const_affinity_matrix,
    lambda_c=0.01,
    print_results=False,
)
"""
