"""
Graphing functionality, including line graphs, histograms and fish plots. Used inside the simulation, and can also be utilised externally.
"""


from discrete_model import Simulation, SimulationState, CellBundle
from config.conf import path_to_data, path_to_output
import numpy as np
import pandas as pd
from typing import Callable
import os


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
    return plt


def fish_tumour(sim: Simulation):
    return fish(sim, "tumour")


def fish_CTL(sim: Simulation):
    return fish(sim, "CTL")


def fish(sim: Simulation, bundle_name):
    from pyfish import fish_plot, process_data, setup_figure
    import matplotlib.pyplot as plt

    if sim.history.state_init_type == "detailed":
        bundles = [get_bundle(state, bundle_name) for state in sim.history]
        populations = []
        phenotypes = set()
        for step in range(len(bundles)):
            bundle = bundles[step]
            for phenotype, pop in bundle.cells_at_phenotype.items():
                populations.append([phenotype, step, pop])
                phenotypes.add(phenotype)
        populations_df = pd.DataFrame(
            np.array(populations), columns=["Id", "Step", "Pop"]
        )
        parent_tree = []
        """
        We need to construct the parent tree. I'm not sure quite how this works with mutations, since the parent could be any from the left/ right.
        This means no phenotypes have no parents.. And hmm, one phenotype can have multiple phenotypes.
        To generate this, I'd want to collect all present phenotypes, then order chronologically, and assign parents in descending order.
        """
        ordered_phenotypes = list(phenotypes)
        ordered_phenotypes.sort()
        parent = None
        child = None
        for i in range(len(ordered_phenotypes)):
            child = ordered_phenotypes[i]
            if parent is not None:
                parent_tree.append([parent, child])
            parent = child
        parent_tree_df = pd.DataFrame(
            np.array(parent_tree), columns=["ParentId", "ChildId"]
        )
        data = process_data(populations_df, parent_tree_df)
        setup_figure()
        fish_plot(*data)
        return plt


def get_bundle(state: SimulationState, bundle_name) -> CellBundle:
    if bundle_name == "tumour":
        return state.tumour_cells
    elif bundle_name == "CTL":
        return state.CTL_cells
    else:
        return None


def graph_from_path(path=path_to_data):
    sim = get_sim(path)
    graph(sim).show()


def flatten_dict(dict: dict):
    return [key for key, val in dict.items() for _ in range(val)]


def hist_base(tumour_cells: CellBundle, CTL_cells: CellBundle, phen_struct):
    import matplotlib.pyplot as plt

    tumour_cell_phenotypes = flatten_dict(tumour_cells.cells_at_phenotype)
    CTL_cell_phenotypes = flatten_dict(CTL_cells.cells_at_phenotype)

    plt.hist(
        CTL_cell_phenotypes,
        label="CTL Cells",
        bins=phen_struct.no_possible_values,
        alpha=0.6,
    )
    plt.hist(
        tumour_cell_phenotypes,
        label="Tumour Cells",
        bins=phen_struct.no_possible_values,
        alpha=0.6,
    )
    plt.legend()
    return plt


def hist_from_state(state: SimulationState):
    return hist_base(
        state.tumour_cells, state.CTL_cells, state.tumour_cells.phen_struct
    )


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
    return plt


plt_fn_label = {
    graph: "graph",
    hist: "hist",
    fish_tumour: "fish_tumour",
    fish_CTL: "fish_CTL",
}


def savefig(
    sim: Simulation = None, plt_fn: Callable = graph, path_to_output=path_to_output
):
    if sim is None:
        sim = get_sim()
    plt = plt_fn(sim)
    out_path = f"{path_to_output}output_{sim.config_name}_{plt_fn_label[plt_fn]}.png"
    # I'm not going to worry about removing the inputs here, because this should never have a file in it, because it's a completely empty directory!
    while os.path.exists(out_path):
        overwrite = input("Would you like to overwrite the existing file?")
        if overwrite == "y":
            break
        else:
            still_save = input("Would you still like to save the file?")
            if still_save == "y":
                out_path = input("Enter the new path:")
                continue
            else:
                print("The plot has not been saved.")
                return
    plt.savefig(out_path)


def hist_from_path(path=path_to_data):
    sim = get_sim(path)
    hist(sim).show()
