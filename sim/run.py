"""
Running my discrete model simulation (and setting it up).
"""

from discrete_model import Simulation, UniversalCellParams
from config.conf import path_to_data
from inputs import get_sim_configuration, read_phenotypes, get_matrix_function_from_config


def create_simulation(config_name=None):
    cf = get_sim_configuration("discrete", config_name)
    print(config_name)

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

    if cf.subtype == "lattice":
        sim = Simulation(
            cf.time_step,
            cf.final_time,
            no_init_tumour_cells=int(cf.no_init_tumour_cells),
            no_init_CTL_cells=int(cf.no_init_CTL_cells),
            tumour_universal_params=tumour_universal_params,
            CTL_universal_params=CTL_universal_params,
            TCR_affinity_range=cf.affinity_range,
            tumour_phenotypic_variation_probability=cf.tumour_phenotypic_variation_probability,
            config_name=cf.name,
            no_possible_phenotypes=int(cf.no_possible_phenotypes),
            absolute_max_phenotype=1,
            TCR_binding_affinity=cf.binding_affinity
        )
    elif cf.subtype == "sequence":
        CTL_sequences = read_phenotypes(cf.CTL_sequence_path)
        tumour_sequences = read_phenotypes(cf.tumour_sequence_path)
        get_sequence_matrix = get_matrix_function_from_config(cf.sequence_matrix_config)
        get_affinity_matrix = get_matrix_function_from_config(cf.affinity_matrix_config)

        sim = Simulation(
            cf.time_step,
            cf.final_time,
            no_init_tumour_cells=int(cf.no_init_tumour_cells),
            no_init_CTL_cells=int(cf.no_init_CTL_cells),
            tumour_universal_params=tumour_universal_params,
            CTL_universal_params=CTL_universal_params,
            TCR_affinity_range=cf.affinity_range,
            tumour_phenotypic_variation_probability=cf.tumour_phenotypic_variation_probability,
            config_name=cf.name,
            CTL_sequences=CTL_sequences,
            tumour_sequences=tumour_sequences,
            get_sequence_matrix=get_sequence_matrix,
            get_affinity_matrix=get_affinity_matrix
        )

    else:
        raise NotImplementedError(f"Subtype {cf.subtype} has not been implemented into run.py. Add this functionality, or check your spelling.")
    return sim


def run(overwrite=None, config_name=None):
    create_new = False
    try:
        sim = Simulation.load_simulation(path_to_data)
        if overwrite is None:
            overwrite = input("Would you like to overwrite this simulation?")

        if overwrite == "y":
            create_new = True
    except IOError:
        create_new = True

    if create_new:
        print("Creating a new simulation.")
        sim = create_simulation(config_name)
    sim.run()

    Simulation.save_simulation(path_to_data, sim)


def extend(additional_time):
    """
    Extend the end time of the simulation.
    """
    try:
        sim = Simulation.load_simulation(path_to_data)
        sim.extend(additional_time)
        Simulation.save_simulation(path_to_data, sim)
        print("Simulation extended.")
    except IOError:
        print("Simuation to extend could not be found")
