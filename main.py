import random

time_step = 0.05
final_time = 30 * 24 * 60 * 60
no_possible_phenotypes = 1500


class Cell:
    def __init__(self, phenotype, simulation):
        self.phenotype = phenotype
        self.simulation = (
            simulation  # Just point to the simulation which owns this
        )

    @classmethod
    def random(self, simulation):
        return Cell(simulation.get_random_phenotype, simulation)


class Cells:
    def __init__(self, cells):
        self.cells = cells


class SimulationState:
    def __init__(self):
        pass


class UniversalCellParams:
    def __init__(
        self,
        natural_prolif_rate,
        natural_death_rate,
        interaction_induced_rate,
        selectivity,
    ):
        # Death rate will be negative for consistency
        self.natural_prolif_rate = natural_prolif_rate
        self.natural_death_rate = natural_death_rate
        self.interaction_induced_rate = interaction_induced_rate
        self.selectivity = selectivity


class Simulation:
    def __init__(
        self,
        time_step_size,
        final_time,
        no_possible_phenotypes,
        max_absolute_phenotype,
        no_tumour_cells,
        no_CTL_cells,
        tumour_universal_params,
        CTL_universal_params,
        TCR_affinity_range,
        TCR_binding_affinity,
        tumour_phenotypic_variation_probability,
    ):
        self.time_step_size = time_step_size
        self.time_step = 0  # An integer describing which time step we're on
        self.final_time_step = final_time / time_step_size

        self.max_absolute_phenotype = max_absolute_phenotype
        self.no_possible_phenotypes = no_possible_phenotypes
        self.phenotype_step_size = (
            2 * max_absolute_phenotype / no_possible_phenotypes  # check this
        )
        self.phenotype_range = range(
            -self.max_absolute_phenotype,
            self.max_absolute_phenotype,
            self.phenotype_step_size,
        )  # check
        self.history = []
        self.cells = []

        self.no_tumour_cells = no_tumour_cells
        self.no_CTL_cells = no_CTL_cells

    def get_immune_score(self):
        return self.no_CTL_cells / self.no_tumour_cells

    def get_average_immune_score(self):
        pass

    def get_random_phenotype(self):
        return (
            random.randint(0, self.no_possible_phenotypes - 1)
            * self.phenotype_step_size
        ) - self.max_absolute_phenotype

    def run(self):
        self.time_step += 1

        if self.time_step > self.final_time_step:
            # End it
            pass
        pass
