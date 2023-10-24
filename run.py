from simulation import Simulation, UniversalCellParams
from conf import path_to_data, get_sim_configuration


def create_simulation():
    cf = get_sim_configuration()
    """
    time_step = 0.05
    final_time = 30
    no_possible_phenotypes = 101

    no_init_tumour_cells = 10000
    no_init_CTL_cells = 10000

    cf.tumour_natural_prolif_rate = 1.5
    cf.tumour_natural_death_rate = 1.5e-6
    cf.tumour_interaction_induced_rate = 5e-6

    cf.CTL_natural_prolif_rate = 1.5
    cf.CTL_natural_death_rate = 1.5e-6
    cf.CTL_interaction_induced_rate = 5e-6

    tumour_selectivity = 0.1  # in [0.1,2]
    CTL_selectivity = 0.1  # in [0.1,2]
    binding_affinity = 0.1  # in [0.1, 3.5]
    affinity_range = 0.1  # in [0.1,2]
    """

    cf.no_possible_phenotypes = int(cf.no_possible_phenotypes)

    tumour_universal_params = UniversalCellParams(
        cf.tumour_natural_prolif_rate,
        cf.tumour_natural_death_rate,
        cf.tumour_interaction_induced_rate,
        cf.tumour_selectivity,
    )

    CTL_universal_params = UniversalCellParams(
        cf.CTL_natural_prolif_rate,
        cf.CTL_natural_death_rate,
        cf.CTL_interaction_induced_rate,
        cf.CTL_selectivity,
    )

    sim = Simulation(
        cf.time_step,
        cf.final_time,
        int(cf.no_possible_phenotypes),
        1,
        no_init_tumour_cells=int(cf.no_init_tumour_cells),
        no_init_CTL_cells=int(cf.no_init_CTL_cells),
        tumour_universal_params=tumour_universal_params,
        CTL_universal_params=CTL_universal_params,
        TCR_affinity_range=cf.affinity_range,
        TCR_binding_affinity=cf.binding_affinity,
        tumour_phenotypic_variation_probability=cf.tumour_phenotypic_variation_probability,
    )
    return sim


def run():
    try:
        sim = Simulation.load_simulation(path_to_data)
    except IOError:
        print("Creating a new simulation.")
        sim = create_simulation()

    sim.run()

    Simulation.save_simulation(path_to_data, sim)
