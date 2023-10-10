import random
import numpy as np
from copy import deepcopy


class PhenotypeStructure:
    def __init__(self, abs_max_value, no_possible_values):
        self.abs_max_value = abs_max_value
        self.range = np.linspace(
            -self.abs_max_value, self.abs_max_value, no_possible_values, True
        )
        print(len(self.range))
        self.step_size = 2 * abs_max_value / no_possible_values  # check this
        self.no_possible_values = no_possible_values

    def shift(self, phen, no_steps, direction):
        print(phen)
        phen += no_steps * direction * self.step_size
        if phen > self.abs_max_value:
            phen = self.abs_max_value
        if phen < -self.abs_max_value:
            phen = -self.abs_max_value
        return phen

    def get_random_phenotype(self):
        return (
            random.randint(0, self.no_possible_values - 1) * self.step_size
        ) - self.abs_max_value


class Cell:
    def __init__(self, phenotype, phen_struct: PhenotypeStructure):
        self.phenotype = phenotype
        self.phen_struct = phen_struct

    def mutate_int(self, no_steps, direction):
        self.phenotype = self.phen_struct.shift(
            self.phenotype, no_steps, direction
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
        cells: list[Cell],
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
            if cell.phenotype in self.no_cells_at_phenotype:
                self.no_cells_at_phenotype[cell.phenotype] += 1
            else:
                self.no_cells_at_phenotype[cell.phenotype] = 0

    @classmethod
    def random(
        self,
        number,
        universal_params: UniversalCellParams,
        phen_struct: PhenotypeStructure,
    ):
        # Create some number of random cells
        cells = []
        for i in range(number):
            cells.append(Cell.random(phen_struct))
        return Cells(cells, universal_params)

    @classmethod
    def evolve_population(self, cells, get_phenotype_probabilities):
        """
        phenotype_probabilities has the birth, death and quiescence probabilities of the population
        """
        for cell in cells.cells:
            action_name = random.choices(
                population=["birth", "death", "quiescence"],
                weights=get_phenotype_probabilities(cell.phenotype),
            )
            if action_name == "birth":
                new_cell = deepcopy(cell)
                cells.cells.append(new_cell)
            elif action_name == "death":
                cells.cells.remove(cell)


class SimulationState:
    """For recording history."""

    def __init__(self):
        pass


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
        self.final_time_step = final_time / time_step_size

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

        self.history = []

    def get_immune_score(self):
        return len(self.CTL_cells.cells) / len(self.tumour_cells.cells)

    def get_average_immune_score(self):
        pass

    def get_phenotypic_separation_scaling(
        self, phenotype_1, phenotype_2, range
    ):
        if abs(phenotype_1 - phenotype_2) <= range:
            return 1 / (
                min(phenotype_1 + range, self.phen_struct.abs_max_value)
                - max(phenotype_1 - range, -self.phen_struct.abs_max_value)
            )
        else:
            return 0

    def get_phenotype_natural_death_rate(self, cells: Cells, phenotype):
        # Based on death base rate, and a weighted sum of the competition from "close species"
        return cells.universal_params.natural_death_base_rate * sum(
            [
                self.get_phenotypic_separation_scaling(
                    phenotype,
                    other_phenotype,
                    cells.universal_params.selectivity,
                )
                * cells.no_cells_at_phenotype[phenotype]
                for other_phenotype in self.phen_struct.range
            ]
        )

    def get_phenotype_interaction_induced_rate(
        self, cells: Cells, other_cells: Cells, phenotype
    ):
        # The rate of growth/ death resulting from the interaction of two sets of cells (tumour and CTL)
        return (
            cells.universal_params.interaction_induced_base_rate
            * self.TCR_binding_affinity
            * sum(
                [
                    self.get_phenotypic_separation_scaling(
                        phenotype, other_phenotype, self.TCR_affinity_range
                    )
                    * other_cells.no_cells_at_phenotype[phenotype]
                    for other_phenotype in self.phen_struct.range
                ]
            )
        )

    def run(self):
        while self.time_step < self.final_time_step:
            self.time_step += 1

            # Pre-calculation
            self.tumour_cells.compute_cells_at_each_phenotype()
            self.CTL_cells.compute_cells_at_each_phenotype()

            # Simulating effects
            for tumour_cell in self.tumour_cells.cells:
                r_1 = random.randrange(0, 1)
                r_2 = random.randrange(0, 1)

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

            # Post-calculation

            # End it

    def get_phenotype_tumour_probabilities(self, phenotype):
        if phenotype not in self.phenotype_tumour_probabilities:
            self.phenotype_tumour_probabilities[
                phenotype
            ] = self.compute_phenotype_tumour_probabilities(phenotype)

        return self.phenotype_tumour_probabilities[phenotype]

    def get_phenotype_CTL_probabilities(self, phenotype):
        if phenotype not in self.phenotype_CTL_probabilities:
            self.phenotype_CTL_probabilities[
                phenotype
            ] = self.compute_phenotype_CTL_probabilities(phenotype)

        return self.phenotype_tumour_probabilities[phenotype]

    def compute_phenotype_tumour_probabilities(self, phenotype):
        birth = (
            self.time_step_size
            * self.tumour_cells.universal_params.natural_prolif_base_rate
        )
        death = self.time_step_size * (
            self.get_phenotype_natural_death_rate(self.tumour_cells, phenotype)
            + self.get_phenotype_interaction_induced_rate(
                self.tumour_cells, self.CTL_cells
            ),
            phenotype,
        )
        return birth, death, 1 - (birth + death)

    def compute_phenotype_CTL_probabilities(self, phenotype):
        birth = self.time_step_size * (
            self.CTL_cells.universal_params.natural_prolif_base_rate
            + self.get_phenotype_interaction_induced_rate(
                self.CTL_cells, self.tumour_cells
            )
        )
        death = self.time_step_size * self.get_phenotype_natural_death_rate(
            self.CTL_cells, phenotype
        )
        return birth, death, 1 - (birth + death)
