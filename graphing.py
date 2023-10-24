from simulation import Simulation
from conf import path_to_data


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
    plt.plot(tumour_cells_pop)
    plt.plot(CTL_cells_pop)
    plt.show()
