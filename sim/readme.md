## Configurations
This section describes what parameters need to be specified to run the various simulations.

### Marta's Continuous Model

Required:

time_step
final_time

no_possible_phenotypes
A
a

tumour_natural_prolif_rate
tumour_natural_death_rate
tumour_interaction_induced_rate
tumour_selectivity

CTL_natural_prolif_rate
CTL_natural_death_rate
CTL_interaction_induced_rate
CTL_selectivity

affinity_range
tumour_phenotypic_variation_probability

Optional:
name

Not required (could be made optional):
binding_affinity_matrix_config
print_results?


## Discrete Model

### Lattice

Required:

name

time_step
final_time

no_possible_phenotypes
no_init_tumour_cells
no_init_CTL_cells

affinity_range
binding_affinity
tumour_phenotypic_variation_probability

tumour_natural_prolif_rate
tumour_natural_death_rate
tumour_interaction_induced_rate
tumour_selectivity

CTL_natural_prolif_rate
CTL_natural_death_rate
CTL_interaction_induced_rate
CTL_selectivity