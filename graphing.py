from simulation import Simulation
from conf import path_to_data


def get_sim(path=path_to_data) -> Simulation:
    try:
        sim = Simulation.load_simulation(path_to_data)
        return sim
    except IOError:
        print("Simulation not found.")
        return
