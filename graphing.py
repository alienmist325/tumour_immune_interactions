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


def graph(path=path_to_data):
    import matplotlib.pyplot as plt

    sim = get_sim(path)
    tumour_cells_pop, CTL_cells_pop = get_pops(sim)
    times = np.linspace(0, sim.final_time, sim.final_time_step)
    plt.plot(times, tumour_cells_pop, label="Tumour Cells")
    plt.plot(times, CTL_cells_pop, label="CTL Cells")
    plt.legend()
    plt.show()
