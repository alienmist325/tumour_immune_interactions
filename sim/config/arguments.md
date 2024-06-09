# Which arguments are required, optional and ignored for each simulation type?

## Continuous
### Marta's Continuous Model

--required--

time_step
final_time
no_possible_phenotypes
tumour_selectivity
CTL_selectivity
affinity_range
tumour_natural_prolif_rate
tumour_natural_death_rate
tumour_interaction_induced_rate
CTL_natural_prolif_rate
CTL_natural_death_rate
CTL_interaction_induced_rate
tumour_phenotypic_variation_probability
A
a

--optional--


--ignored--

binding_affinity
no_init_tumour_cells
no_init_CTL_cells
sequence_matrix_config
affinity_matrix_config
CTL_sequence_path
tumour_sequence_path

## Discrete Model

### Lattice

--required--

time_step
final_time
no_possible_phenotypes
tumour_selectivity
CTL_selectivity
binding_affinity
affinity_range
no_init_tumour_cells
no_init_CTL_cells
tumour_natural_prolif_rate
tumour_natural_death_rate
tumour_interaction_induced_rate
CTL_natural_prolif_rate
CTL_natural_death_rate
CTL_interaction_induced_rate
tumour_phenotypic_variation_probability

--optional--


--ignored--

A
a
sequence_matrix_config
affinity_matrix_config
CTL_sequence_path
tumour_sequence_path

### Sequence



--required--

time_step
final_time
no_possible_phenotypes
tumour_selectivity
CTL_selectivity
affinity_range
no_init_tumour_cells
no_init_CTL_cells
tumour_natural_prolif_rate
tumour_natural_death_rate
tumour_interaction_induced_rate
CTL_natural_prolif_rate
CTL_natural_death_rate
CTL_interaction_induced_rate
tumour_phenotypic_variation_probability
sequence_matrix_config
affinity_matrix_config
CTL_sequence_path
tumour_sequence_path

--optional--


--ignored--

binding_affinity
A
a
