## Configurations
This section describes what parameters need to be specified to run the various simulations.

The configuration loaded is dependent on:
- `simtype` (this will come from the file you run: `continuous_run.py` or `run.py` i.e. whether you run Marta's model, or the discrete model)
`discrete`
`continuous`

defaults to `discrete`

- `subtype` (specified in the config file, and defaulting to `lattice`)
`lattice`
`sequence`

`name` must be unique.


All arguments ar either:
- required (0)
- optional (1)
- ignored (2)

Below, we list (excluding `simtype`, `subtype` and `name`) all the arguments in each category for all the currently supported simulations.

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

## Matrix configuration
This will be done as a json file. 

### Path
- "from" : "path"
- "path" : from the source directory (simulation)
- "delimiter" : If empty "" or unspecified, then we assume a space " " delimiting

Internally, this will load the matrix, have it stored as the output of a function, and pass this function into the simulation.

### Function
- "from" : "function"
- "where" : if empty "", then we assume globally in the source file (`Simulation`); otherwise we import that relative reference
- "function" : the name of the function (cannot be inside a class). This function takes in the `Simulation` object after all initialisations have been done (so you may need to inspect the init file), and returns a matrix 


### Aside
Affinity matrix needs to be width CTL, height tumour