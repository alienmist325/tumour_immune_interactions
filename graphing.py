from simulation import Simulation
from conf import path_to_data
import numpy as np


def get_sim(path=path_to_data) -> Simulation:
    try:
        sim = Simulation.load_simulation(path)
        return sim
    except IOError:
        print("Simulation not found.")
        return


def get_pops(sim: Simulation):
    tumour_cells_pop = []
    CTL_cells_pop = []
    for state in sim.history.history:
        tumour_cells_pop.append(state.tumour_cells_pop)
        CTL_cells_pop.append(state.CTL_cells_pop)
    return tumour_cells_pop, CTL_cells_pop


def graph(sim: Simulation):
    import matplotlib.pyplot as plt

    tumour_cells_pop, CTL_cells_pop = get_pops(sim)
    times = np.linspace(0, sim.time_step * sim.time_step_size, sim.time_step)
    plt.plot(times, CTL_cells_pop, label="CTL Cells")
    plt.plot(times, tumour_cells_pop, label="Tumour Cells")
    plt.legend()
    plt.show()


def graph_from_path(path=path_to_data):
    sim = get_sim(path)
    graph(sim)


def flatten_dict(dict: dict):
    return [key for key, val in dict.items() for _ in range(val)]


def hist(sim: Simulation):
    import matplotlib.pyplot as plt

    tumour_cell_phenotypes = flatten_dict(sim.tumour_cells.cells_at_phenotype)
    CTL_cell_phenotypes = flatten_dict(sim.CTL_cells.cells_at_phenotype)

    plt.hist(
        CTL_cell_phenotypes,
        label="CTL Cells",
        bins=sim.phen_struct.no_possible_values,
        alpha=0.6,
    )
    plt.hist(
        tumour_cell_phenotypes,
        label="Tumour Cells",
        bins=sim.phen_struct.no_possible_values,
        alpha=0.6,
    )
    plt.legend()
    plt.show()


def hist_from_path(path=path_to_data):
    sim = get_sim(path)
    hist(sim)
