from simulation import Simulation, UniversalCellParams

time_step = 0.05
final_time = 30
no_possible_phenotypes = 1500

tumour_selectivity = 0.1  # in [0.1,2]
CTL_selectivity = 0.1  # in [0.1,2]
binding_affinity = 0.1  # in [0.1, 3.5]
affinity_range = 0.1  # in [0.1,2]


sim = Simulation(
    time_step,
    final_time,
    no_possible_phenotypes,
    1,
    no_init_tumour_cells=10000,
    no_init_CTL_cells=20000,
    tumour_universal_params=UniversalCellParams(
        1.5, 1.5 * (10**-6), 5 * (10**-6), tumour_selectivity
    ),
    CTL_universal_params=UniversalCellParams(
        5 * (10**-2), 5 * (10**-6), 3 * (10**-5), CTL_selectivity
    ),
    TCR_affinity_range=affinity_range,
    TCR_binding_affinity=binding_affinity,
    tumour_phenotypic_variation_probability=0.01,
)


print("begin")
sim.run()
print("done")
