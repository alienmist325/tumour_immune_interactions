## Configurations
This section describes what parameters need to be specified to run the various simulations.

The configuration loaded is dependent on:
- `simtype` (this will come from the file you run: `continuous_run.py` or `run.py` i.e. whether you run Marta's model, or the discrete model)
`discrete`
`continuous`
- `subtype` (specified in the config file, and defaulting to `lattice`)

`name` must be unique.


All arguments ar either:
- required (0)
- optional (1)
- ignored (2)

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
affinity_matrix_config
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

> subtype

### Sequence


Required:

name

time_step
final_time

no_possible_phenotypes
no_init_tumour_cells
no_init_CTL_cells

affinity_range
// binding_affinity
tumour_phenotypic_variation_probability

tumour_natural_prolif_rate
tumour_natural_death_rate
tumour_interaction_induced_rate
tumour_selectivity

CTL_natural_prolif_rate
CTL_natural_death_rate
CTL_interaction_induced_rate
CTL_selectivity

> subtype
> sequence_matrix_config
> affinity_matrix_config


## Matrix configuration
This will be done as a json file. 

### Path
- "from" : "path"
- "path" : from the source directory (simulation)
- "delimiter" : If empty "" or unspecified, then we assume a space " " delimiting

Internally, this will load the matrix, have it stored as the output of a function, and pass this function into the simulation.

# Function
- "from" : "function"
- "where" : if empty "", then we assume globally in the source file (`Simulation`); otherwise we import that relative reference
- "function" : the name of the function (cannot be inside a class). This function takes in the `Simulation` object after all initialisations have been done (so you may need to inspect the init file), and returns a matrix 