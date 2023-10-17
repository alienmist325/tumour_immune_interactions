import random
import numpy as np
from copy import deepcopy
import conf
import importlib
import pickle


class PhenotypeStructure:
    def __init__(self, abs_max_value, no_possible_values):
        self.abs_max_value = abs_max_value
        self.range = np.linspace(
            -self.abs_max_value, self.abs_max_value, no_possible_values, True
        )
        self.id_range = range(no_possible_values)
        print(len(self.range))
        self.step_size = 2 * abs_max_value / no_possible_values  # check this
        self.no_possible_values = no_possible_values

    """
    def shift(self, phen, no_steps, direction):
        print(phen)
        phen += no_steps * direction * self.step_size
        if phen > self.abs_max_value:
            phen = self.abs_max_value
        if phen < -self.abs_max_value:
            phen = -self.abs_max_value
        return phen
    """

    def shift(self, phen_id, no_steps, direction):
        phen_id += no_steps * direction
        if phen_id > self.no_possible_values - 1:
            phen_id = self.no_possible_values - 1
        elif phen_id < 0:
            phen_id = 0
        return phen_id

    def get_phenotype_by_id(self, id):
        return (id * self.step_size) - self.abs_max_value

    def get_random_phenotype(self):
        return random.randint(0, self.no_possible_values - 1)

    """
    def get_random_phenotype(self):
        return (
            random.randint(0, self.no_possible_values - 1) * self.step_size
        ) - self.abs_max_value
    """


class Cell:
    def __init__(self, phenotype_id, phen_struct: PhenotypeStructure):
        self.phenotype_id = phenotype_id
        self.phen_struct = phen_struct

    def mutate_int(self, no_steps, direction):
        self.phenotype_id = self.phen_struct.shift(
            self.phenotype_id, no_steps, direction
        )

    def mutate_left(self):
        self.mutate_int(1, -1)

    def mutate_right(self):
        self.mutate_int(1, 1)

    @classmethod
    def random(self, phenotype_structure: PhenotypeStructure):
        return Cell(
            phenotype_structure.get_random_phenotype(), phenotype_structure
        )


class UniversalCellParams:
    def __init__(
        self,
        natural_prolif_rate,
        natural_death_rate,
        interaction_induced_rate,
        selectivity,
    ):
        # Death rate will be negative for consistency
        self.natural_prolif_base_rate = natural_prolif_rate
        self.natural_death_base_rate = natural_death_rate
        self.interaction_induced_base_rate = interaction_induced_rate
        self.selectivity = selectivity


class Cells:
    def __init__(
        self,
        cells: set[Cell],
        universal_params: UniversalCellParams,
    ):
        self.cells = cells
        self.no_cells_at_phenotype = {}
        self.universal_params = universal_params

    def compute_cells_at_each_phenotype(self):
        self.no_cells_at_phenotype = (
            {}
        )  # Reset dictionary (inefficient, but clear)

        for cell in self.cells:
            if cell.phenotype_id in self.no_cells_at_phenotype:
                self.no_cells_at_phenotype[cell.phenotype_id] += 1
            else:
                self.no_cells_at_phenotype[cell.phenotype_id] = 1

    def get_no_cells_at_phenotype(self, phenotype_id):
        if phenotype_id in self.no_cells_at_phenotype:
            return self.no_cells_at_phenotype[phenotype_id]
        else:
            return 0

    def __len__(self):
        return len(self.cells)

    @classmethod
    def random(
        self,
        number,
        universal_params: UniversalCellParams,
        phen_struct: PhenotypeStructure,
    ):
        # Create some number of random cells
        cells = set()
        for i in range(number):
            cell = Cell.random(phen_struct)
            cells.add(cell)
        return Cells(cells, universal_params)

    @classmethod
    def evolve_population(
        self,
        cells,
        get_phenotype_probabilities,
    ):
        """
        phenotype_probabilities has the birth, death and quiescence probabilities of the population
        """
        new_cells = set()
        dead_cells = set()
        for cell in cells.cells:
            action_name = random.choices(
                population=["birth", "death", "quiescence"],
                weights=get_phenotype_probabilities(cell.phenotype_id),
            )[0]
            """
            if action_name != "quiescence":
                print(action_name)
            """
            if action_name == "birth":
                new_cell = Cell(cell.phenotype_id, cell.phen_struct)
                new_cells.add(new_cell)
            elif action_name == "death":
                dead_cells.discard(cell)
        for cell in dead_cells:
            cells.cells.discard(cell)
        for cell in new_cells:
            if cell in cells.cells:
                print("Already here")
            cells.cells.add(cell)


class SimulationState:
    def __init__(self, CTL_cells: Cells, tumour_cells: Cells):
        self.CTL_cells_pop = len(CTL_cells)
        self.tumour_cells_pop = len(tumour_cells)


class SimulationHistory:
    """For recording history."""

    def __init__(self, history: list[SimulationState] = []):
        # Do I need to copy the simulation?
        self.history = history

    def update(self, sim_state: SimulationState):
        self.history.append(sim_state)

        # We pickle and unpickle this object to do stuff.


class Simulation:
    def __init__(
        self,
        time_step_size,
        final_time,
        no_possible_phenotypes,
        absolute_max_phenotype,
        no_init_tumour_cells,
        no_init_CTL_cells,
        tumour_universal_params,
        CTL_universal_params,
        TCR_affinity_range,
        TCR_binding_affinity,
        tumour_phenotypic_variation_probability,
    ):
        self.time_step_size = time_step_size
        self.time_step = 0  # An integer describing which time step we're on
        self.final_time_step = int(final_time / time_step_size)

        self.phen_struct = PhenotypeStructure(
            absolute_max_phenotype, no_possible_phenotypes
        )
        self.tumour_cells = Cells.random(
            no_init_tumour_cells, tumour_universal_params, self.phen_struct
        )
        self.CTL_cells = Cells.random(
            no_init_CTL_cells, CTL_universal_params, self.phen_struct
        )

        self.TCR_affinity_range = TCR_affinity_range
        self.TCR_binding_affinity = TCR_binding_affinity
        self.tumour_phenotypic_variation_probability = (
            tumour_phenotypic_variation_probability
        )

        self.phenotype_tumour_probabilities = {}
        self.phenotype_CTL_probabilities = {}

        self.history = SimulationHistory()

    def get_immune_score(self):
        return len(self.CTL_cells.cells) / len(self.tumour_cells.cells)

    def get_average_immune_score(self):
        pass

    def get_phenotypic_separation_scaling(
        self,
        phenotype_1_id,
        phenotype_2_id,
        range,
    ):
        phenotype_1 = self.phen_struct.get_phenotype_by_id(phenotype_1_id)
        phenotype_2 = self.phen_struct.get_phenotype_by_id(phenotype_2_id)
        if abs(phenotype_1 - phenotype_2) <= range:
            return 1 / (
                min(phenotype_1 + range, self.phen_struct.abs_max_value)
                - max(phenotype_1 - range, -self.phen_struct.abs_max_value)
            )
        else:
            return 0

    def get_phenotype_natural_death_rate(self, cells: Cells, phenotype_id):
        # Based on death base rate, and a weighted sum of the competition from "close species"
        return cells.universal_params.natural_death_base_rate * sum(
            [
                self.get_phenotypic_separation_scaling(
                    phenotype_id,
                    other_phenotype_id,
                    cells.universal_params.selectivity,
                )
                * cells.get_no_cells_at_phenotype(phenotype_id)
                for other_phenotype_id in self.phen_struct.id_range
            ]
        )

    def get_phenotype_interaction_induced_rate(
        self,
        cells: Cells,
        other_cells: Cells,
        phenotype_id,
    ):
        # The rate of growth/ death resulting from the interaction of two sets of cells (tumour and CTL)
        return (
            cells.universal_params.interaction_induced_base_rate
            * self.TCR_binding_affinity
            * sum(
                [
                    self.get_phenotypic_separation_scaling(
                        phenotype_id,
                        other_phenotype_id,
                        self.TCR_affinity_range,
                    )
                    * other_cells.get_no_cells_at_phenotype(phenotype_id)
                    for other_phenotype_id in self.phen_struct.id_range
                ]
            )
        )

    def run(self):
        self.print("The simulation is starting.")
        while self.time_step < self.final_time_step:
            importlib.reload(conf)

            if conf.interrupt:
                print(
                    "The simulation has been interrupted and will now safely save."
                )
                return

            self.time_step += 1

            # Pre-calculation
            self.tumour_cells.compute_cells_at_each_phenotype()
            self.CTL_cells.compute_cells_at_each_phenotype()

            # Simulating effects
            for tumour_cell in self.tumour_cells.cells:
                r_1 = random.randrange(0, 1)

                if r_1 < self.tumour_phenotypic_variation_probability:
                    action = random.choice(
                        [tumour_cell.mutate_left, tumour_cell.mutate_right]
                    )
                    action()

            self.phenotype_tumour_probabilities = {}
            self.phenotype_CTL_probabilities = {}

            Cells.evolve_population(
                self.tumour_cells, self.get_phenotype_tumour_probabilities
            )

            Cells.evolve_population(
                self.CTL_cells, self.get_phenotype_CTL_probabilities
            )

            self.print(
                "C: ", len(self.tumour_cells), " | T:", len(self.CTL_cells)
            )
            self.print("Iteration done.")
            self.print(
                "Time step: ", self.time_step, "/", self.final_time_step
            )
            # Post-calculation

            self.history.update(
                SimulationState(self.CTL_cells, self.tumour_cells)
            )

            # End it
        self.print(
            "The final time has been reached, so the simulation is over."
        )

    def print(self, *string):
        if conf.debug:
            print(*string)

    def get_phenotype_tumour_probabilities(self, phenotype_id):
        if phenotype_id not in self.phenotype_tumour_probabilities:
            self.phenotype_tumour_probabilities[
                phenotype_id
            ] = self.compute_phenotype_tumour_probabilities(phenotype_id)

        return self.phenotype_tumour_probabilities[phenotype_id]

    def get_phenotype_CTL_probabilities(self, phenotype_id):
        if phenotype_id not in self.phenotype_CTL_probabilities:
            self.phenotype_CTL_probabilities[
                phenotype_id
            ] = self.compute_phenotype_CTL_probabilities(phenotype_id)

        return self.phenotype_CTL_probabilities[phenotype_id]

    def compute_phenotype_tumour_probabilities(self, phenotype_id):
        birth = (
            self.time_step_size
            * self.tumour_cells.universal_params.natural_prolif_base_rate
        )
        death = self.time_step_size * (
            self.get_phenotype_natural_death_rate(
                self.tumour_cells, phenotype_id
            )
            + self.get_phenotype_interaction_induced_rate(
                self.tumour_cells, self.CTL_cells, phenotype_id
            )
        )
        return birth, death, 1 - (birth + death)

    def compute_phenotype_CTL_probabilities(self, phenotype_id):
        birth = self.time_step_size * (
            self.CTL_cells.universal_params.natural_prolif_base_rate
            + self.get_phenotype_interaction_induced_rate(
                self.CTL_cells, self.tumour_cells, phenotype_id
            )
        )
        death = self.time_step_size * self.get_phenotype_natural_death_rate(
            self.CTL_cells, phenotype_id
        )
        return birth, death, 1 - (birth + death)

    @classmethod
    def load_simulation(self, path_to_data):
        with open(path_to_data, "rb") as f:
            sim = pickle.load(f)
            print("Successfully opened the previous simulation.")

        return sim

    @classmethod
    def save_simulation(self, path_to_data, sim):
        with open(path_to_data, "wb") as f:
            print("Pickling....")
            pickle.dump(sim, f, pickle.HIGHEST_PROTOCOL)
            print("Pickling done.")
