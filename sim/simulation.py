import random
import numpy as np
from copy import deepcopy
import conf
import importlib
import pickle
from typing import Self
from abc import ABC, abstractmethod
from functools import singledispatch


class PhenotypeStructure(ABC):
    @abstractmethod
    def get_random_phenotype(self):
        """
        Generate a valid phenotype id.
        """
        pass

    @abstractmethod
    def get_random_mutation(self, phenotype):
        """
        Get a possible mutation of the phenotype.
        """
        pass

    @abstractmethod
    def get_value(self, phenotype):
        """
        Get the value associated with a phenotype id.
        """
        return phenotype.id  # By default, these are the same

    # IDEA: We want to have a general function here which will take in the two ids and get the scaling. If they're the same type, we should redirect back down to the
    # relevant class, since e.g. LatticePhenotypeStructure _can definitely_ compare between two lattice phenotypes. If they're different, we'll cover it here and overload.

    @classmethod
    @abstractmethod
    def get_interaction_scaling(self, phen_1, phen_2):
        if type(phen_1.struct) is type(phen_2.struct):
            type(phen_1.struct).get_interaction_scaling(phen_1, phen_2)
        else:
            return -1


class Phenotype:
    def __init__(self, struct: PhenotypeStructure, id):
        # Notice we will have a problem if our value isn't one that compares nicely (e.g. float), so we have to use a "comparable" id here instead
        self.struct = struct
        self.id = id

    def __eq__(self, other):
        return (self.struct == other.struct) and (self.id == other.id)

    def get_value(self):
        return self.struct.get_value(self)

    def get_random_mutation(self):
        return self.struct.get_random_mutation(self)

    def __hash__(self):
        return hash((self.struct, self.id))


class SequencePhenotypeStructure(PhenotypeStructure):
    def __init__(self, sequences):
        self.seqeuences = sequences

    def get_random_phenotype(self):
        """
        Generate a valid phenotype id.
        """
        pass

    def get_random_mutation(self, phenotype):
        """
        Get a possible mutation of the phenotype.
        """
        pass

    def get_value(self, phenotype):
        """
        Get the value associated with a phenotype id.
        """
        return phenotype.id  # By default, these are the same

    # IDEA: We want to have a general function here which will take in the two ids and get the scaling. If they're the same type, we should redirect back down to the
    # relevant class, since e.g. LatticePhenotypeStructure _can definitely_ compare between two lattice phenotypes. If they're different, we'll cover it here and overload.

    @classmethod
    def get_affinity_distance(phen_1, phen_2, binding_affinity_matrix):
        # For comparing different seqeuence phenotypes (i.e. tumour and T-cell), using binding affinities
        pass

    @classmethod
    def get_sequence_difference(phen_1, phen_2, sequence_matrix):
        pass

    @classmethod
    def get_interaction_function(self):
        return (
            Self.get_sequence_distance
        )  # or is it better to do self.etc rather than the class method? Which works?

    @classmethod
    def get_standard_interaction_data(self):
        # An example of some data
        # Is this the best way to do this? In truth, we just need a specification of _what_ we need to provide
        sequence_matrix = np.identity(10)  # TODO: fix this
        return sequence_matrix

    @classmethod
    def get_cross_interaction_function(self):
        return Self.get_affinity_distance

    @classmethod
    def get_standard_cross_interaction_data(self):
        # An example of some data
        # Is this the best way to do this? In truth, we just need a specification of _what_ we need to provide
        affinity_matrix = np.identity(10)  # TODO: fix this
        return affinity_matrix


class LatticePhenotypeStructure(PhenotypeStructure):
    """
    Phenotypes are in the range [0,1], and are floats
    Phenotype IDs are in the range [0, no_possible_values - 1], and are integers
    """

    def __init__(self, abs_max_value, no_possible_values):
        self.abs_max_value = abs_max_value
        self.range = np.linspace(
            -self.abs_max_value, self.abs_max_value, no_possible_values, True
        )
        self.id_range = range(no_possible_values)
        self.step_size = 2 * abs_max_value / no_possible_values  # check this
        self.no_possible_values = no_possible_values  # For an id

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

    def shift(self, phen: Phenotype, no_steps, direction):
        phen_id = phen.id
        phen_id += no_steps * direction
        if phen_id > self.no_possible_values - 1:
            phen_id = self.no_possible_values - 1
        elif phen_id < 0:
            phen_id = 0
        return Phenotype(phen.struct, phen_id)

    def get_value(self, phen: Phenotype):
        id = phen.id
        return (id * self.step_size) - self.abs_max_value

    def get_random_phenotype(self):
        return Phenotype(self, random.randint(0, self.no_possible_values - 1))

    def get_random_mutation(self, phenotype: Phenotype):
        no_steps = 1
        direction = random.choice([-1, 1])
        return self.shift(phenotype, no_steps, direction)

    @classmethod
    def is_excluded_phenotype(self, phenotype: Phenotype, exclude_percent=0.1):
        phen_struct = phenotype.struct
        phenotype_id = phenotype.id
        no_values = phen_struct.no_possible_values
        return (
            phenotype_id < no_values * exclude_percent
            or phenotype_id > no_values * (1 - exclude_percent)
        )

    @classmethod
    def get_circular_distance(self, a, b):
        """
        Get the circular distance on [0,1] between a and b.
        """
        if b < a:
            # Switch so a < b
            temp = a
            a = b
            b = temp
        return min(abs(b - a), abs(a + 1 - b))

    @classmethod
    def get_interaction_function(self):
        return (
            Self.compute_interaction_scaling
        )  # or is it better to do self.etc rather than the class method? Which works?

    @classmethod
    def get_standard_interaction_data(self):
        # An example of some data
        # Is this the best way to do this? In truth, we just need a specification of _what_ we need to provide
        distance_type = "length"
        return distance_type

    @classmethod
    def compute_interaction_scaling(
        self, phenotype_1_, phenotype_2_, range, distance_type="line"
    ):
        """
        Return 0 if phenotypes are out of the selectivity range; otherwise return a constant, modified to account for the boundary.
        """
        phenotype_1 = phenotype_1_.get_value()
        phenotype_2 = phenotype_2_.get_value()
        if distance_type == "line":
            distance = abs(phenotype_1 - phenotype_2)
        if distance_type == "circular":
            distance = LatticePhenotypeStructure.get_circular_distance(
                phenotype_1, phenotype_2
            )
        if distance <= range:
            return 1 / (
                min(phenotype_1 + range, self.abs_max_value)
                - max(phenotype_1 - range, -self.abs_max_value)
            )
        else:
            return 0

    """
    def get_random_phenotype(self):
        return (
            random.randint(0, self.no_possible_values - 1) * self.step_size
        ) - self.abs_max_value
    """


class PhenotypeInteractions:
    """
    Outlines how the different PhenotypeStructure objects in the program interact with each other.

    interaction_data : dict from (struct1, struct2) to data
    """

    def __init__(self, interaction_data: dict, interaction_functions: dict):
        self.interaction_data = interaction_data
        self.interaction_functions = interaction_functions
        self.phenotype_separation_scaling = {}

    def compute_interaction_scaling(self, phen_1: Phenotype, phen_2: Phenotype, range):
        struct_tuple = (phen_1.struct, phen_2.struct)
        data = self.interaction_data[struct_tuple]
        function = self.interaction_functions[struct_tuple]
        function(phen_1, phen_2, range, data)

    def get_interaction_scaling(
        self,
        phenotype_1,
        phenotype_2,
        range,
    ):
        if (
            phenotype_1,
            phenotype_2,
        ) not in self.phenotype_separation_scaling:
            self.phenotype_separation_scaling[
                (phenotype_1, phenotype_2)
            ] = self.compute_interaction_scaling(phenotype_1, phenotype_2, range)
        return self.phenotype_separation_scaling[(phenotype_1, phenotype_2)]


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


class CellBundle:
    def __init__(
        self,
        universal_params: UniversalCellParams,
        phen_struct: PhenotypeStructure,
        cells_at_phenotype: dict,
    ):
        self.cells_at_phenotype = cells_at_phenotype
        self.phen_struct = phen_struct
        self.universal_params = universal_params

    def __len__(self):
        return sum(self.cells_at_phenotype.values())

    def create_cells(self, phenotype: Phenotype, number):
        """
        if PhenotypeStructure.is_excluded_phenotype(self.phen_struct, phenotype):
            return
        """

        if phenotype not in self.cells_at_phenotype:
            self.cells_at_phenotype[phenotype] = number
        else:
            self.cells_at_phenotype[phenotype] += number

    def kill_cells(self, phenotype: Phenotype, number):
        if phenotype not in self.cells_at_phenotype:
            raise ValueError("No cells of this phenotype exist. Cannot kill cells.")
        else:
            if self.cells_at_phenotype[phenotype] < number:
                raise ValueError(
                    "Not enough cells of this phenotype exist. Cannot kill cells."
                )
            else:
                self.cells_at_phenotype[phenotype] -= number

    def mutate(self, phenotype: Phenotype, number):
        for i in range(number):
            new_phenotype = phenotype.get_random_mutation()
            self.kill_cells(phenotype, 1)
            self.create_cells(new_phenotype, 1)

    @classmethod
    def random(
        self,
        number,
        universal_params: UniversalCellParams,
        phen_struct: PhenotypeStructure,
    ):
        cell_bundle = CellBundle(universal_params, phen_struct, {})
        for i in range(number):
            cell_bundle.create_cells(phen_struct.get_random_phenotype(), 1)
        return cell_bundle

    @classmethod
    def evolve_population(
        self,
        cells: Self,
        get_phenotype_probabilities,
    ):
        new_cells = deepcopy(cells)
        for phenotype, number in cells.cells_at_phenotype.items():
            # print(number)
            # number = 100 * number
            """
            if PhenotypeStructure.is_excluded_phenotype(
                cells.phen_struct, phenotype
            ):
                continue
            """
            weights = get_phenotype_probabilities(phenotype)
            # print(weights)
            rng = np.random.default_rng()
            births, deaths, quiescences = rng.multinomial(number, weights)
            # print(births, "|", deaths, "|", quiescences)
            new_cells.create_cells(phenotype, births)
            new_cells.kill_cells(phenotype, deaths)
            # print(len(new_cells))
            # Could just subtract and do this in one step
        return new_cells


class SimulationStateTypes:
    @classmethod
    def populations_only(
        self, state: Self, CTL_cells: CellBundle, tumour_cells: CellBundle
    ):
        state.CTL_cells_pop = len(CTL_cells)
        state.tumour_cells_pop = len(tumour_cells)
        return state

    @classmethod
    def whole_cell_bundles(
        self, state: Self, CTL_cells: CellBundle, tumour_cells: CellBundle
    ):
        state.CTL_cells_pop = len(CTL_cells)
        state.tumour_cells_pop = len(tumour_cells)
        state.CTL_cells = CTL_cells
        state.tumour_cells = tumour_cells
        return state


class SimulationState:
    type_to_init_dict = {
        "default": SimulationStateTypes.populations_only,
        "detailed": SimulationStateTypes.whole_cell_bundles,
    }

    def __init__(self, CTL_cells: CellBundle, tumour_cells: CellBundle):
        from conf import sim_state_init_type

        initialiser = SimulationState.type_to_init_dict[sim_state_init_type]
        self = initialiser(self, CTL_cells, tumour_cells)
        self.init_type = sim_state_init_type


class SimulationHistory:
    """For recording history."""

    def __init__(self, history: list[SimulationState] = []):
        # Do I need to copy the simulation?
        from conf import sim_state_init_type

        self.history = history
        self.state_init_type = sim_state_init_type  # We hope this is the same as each sim state, but technically, it could not be

    def update(self, sim_state: SimulationState):
        self.history.append(sim_state)

        # We pickle and unpickle this object to do stuff.

    def __iter__(self):
        return self.history.__iter__()


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
        config_name="Unspecified",
    ):
        self.config_name = config_name
        self.final_time = final_time
        if conf.m_adjustment:
            scalar = 1  # 1.15, 2
            time_step_size *= scalar
            final_time *= scalar
        self.time_step_size = time_step_size
        self.time_step = 0  # An integer describing which time step we're on
        self.final_time_step = int(final_time / time_step_size)

        self.phen_struct = PhenotypeStructure(
            absolute_max_phenotype, no_possible_phenotypes
        )
        self.tumour_cells = CellBundle.random(
            no_init_tumour_cells, tumour_universal_params, self.phen_struct
        )
        self.CTL_cells = CellBundle.random(
            no_init_CTL_cells, CTL_universal_params, self.phen_struct
        )

        self.TCR_affinity_range = TCR_affinity_range
        self.TCR_binding_affinity = TCR_binding_affinity
        self.tumour_phenotypic_variation_probability = (
            tumour_phenotypic_variation_probability
        )

        self.phenotype_tumour_probabilities = {}
        self.phenotype_CTL_probabilities = {}

        self.phenotype_separation_scaling = {}

        self.history = SimulationHistory()
        self.temp_scalar = 1

        # Something like this (implement properly)
        sequence_matrix = np.identity(10)
        affinity_matrix = np.identity(10)

        self.tumour_struct = SequencePhenotypeStructure([])
        self.CTL_struct = SequencePhenotypeStructure([])

        ts = self.tumour_struct
        cs = self.CTL_struct
        interaction_data = {}

        interaction_data[(ts, ts)] = sequence_matrix
        interaction_data[(cs, cs)] = sequence_matrix
        interaction_data[(ts, cs)] = affinity_matrix
        interaction_data[(cs, ts)] = np.transpose(affinity_matrix)

        interaction_function = {}
        interaction_function[
            (ts, ts)
        ] = SequencePhenotypeStructure.get_interaction_function()
        interaction_function[
            (cs, cs)
        ] = SequencePhenotypeStructure.get_interaction_function()
        interaction_function[
            (ts, cs)
        ] = SequencePhenotypeStructure.get_cross_interaction_function()
        interaction_function[
            (cs, ts)
        ] = SequencePhenotypeStructure.get_cross_interaction_function()

        self.phen_int = PhenotypeInteractions(
            interaction_functions=interaction_function,
            interaction_data=interaction_data,
        )

    def get_immune_score(self):
        return len(self.CTL_cells.cells) / len(self.tumour_cells.cells)

    def get_average_immune_score(self):
        pass

    def get_phenotype_natural_death_rate(self, cells: CellBundle, phenotype: Phenotype):
        # Based on death base rate, and a weighted sum of the competition from "close species"
        return cells.universal_params.natural_death_base_rate * sum(
            [
                self.phen_int.get_interaction_scaling(
                    phenotype,
                    other_phenotype,
                    cells.universal_params.selectivity,
                )
                * cells_at_phenotype
                for other_phenotype, cells_at_phenotype in cells.cells_at_phenotype.items()
            ]
        )

    def get_phenotype_interaction_induced_rate(
        self,
        cells: CellBundle,
        other_cells: CellBundle,
        phenotype: Phenotype,
    ):
        # The rate of growth/ death resulting from the interaction of two sets of cells (tumour and CTL)
        return (
            cells.universal_params.interaction_induced_base_rate
            * self.TCR_binding_affinity
            * sum(
                [
                    self.phen_int.get_interaction_scaling(
                        phenotype,
                        other_phenotype,
                        self.TCR_affinity_range,
                    )
                    * other_cells_at_phenotype
                    for other_phenotype, other_cells_at_phenotype in other_cells.cells_at_phenotype.items()
                ]
            )
        )

    def mutate(self, cells: CellBundle):
        new_cells = deepcopy(cells)
        for phenotype, number in cells.cells_at_phenotype.items():
            rng = np.random.default_rng()
            no_mutations = rng.binomial(
                number, self.tumour_phenotypic_variation_probability
            )
            new_cells.mutate(phenotype, no_mutations)
        return new_cells

    def run(self):
        self.print("The simulation is starting.")
        while self.time_step < self.final_time_step:
            importlib.reload(conf)

            if conf.interrupt:
                print("The simulation has been interrupted and will now safely save.")
                return

            self.time_step += 1

            # Pre-calculation

            # Simulating effects
            self.tumour_cells = self.mutate(self.tumour_cells)

            self.phenotype_tumour_probabilities = {}
            self.phenotype_CTL_probabilities = {}

            self.tumour_cells = CellBundle.evolve_population(
                self.tumour_cells, self.get_phenotype_tumour_probabilities
            )
            self.CTL_cells = CellBundle.evolve_population(
                self.CTL_cells, self.get_phenotype_CTL_probabilities
            )

            self.print("C: ", len(self.tumour_cells), " | T:", len(self.CTL_cells))
            self.print("Iteration done.")
            self.print("Time step: ", self.time_step, "/", self.final_time_step)
            # Post-calculation

            self.history.update(SimulationState(self.CTL_cells, self.tumour_cells))

            # End it
        self.print("The final time has been reached, so the simulation is over.")

    def print(self, *string):
        if conf.debug:
            print(*string)

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

        return self.phenotype_CTL_probabilities[phenotype]

    def compute_phenotype_tumour_probabilities(self, phenotype):
        birth = (
            self.time_step_size
            * self.tumour_cells.universal_params.natural_prolif_base_rate
        )
        death = self.time_step_size * (
            self.get_phenotype_natural_death_rate(self.tumour_cells, phenotype)
            + self.get_phenotype_interaction_induced_rate(
                self.tumour_cells, self.CTL_cells, phenotype
            )
        )
        birth = self.temp_scalar * birth
        death = self.temp_scalar * death
        return birth, death, 1 - (birth + death)

    def compute_phenotype_CTL_probabilities(self, phenotype):
        birth = self.time_step_size * (
            self.CTL_cells.universal_params.natural_prolif_base_rate
            + self.get_phenotype_interaction_induced_rate(
                self.CTL_cells, self.tumour_cells, phenotype
            )
        )
        death = self.time_step_size * self.get_phenotype_natural_death_rate(
            self.CTL_cells, phenotype
        )
        birth = self.temp_scalar * birth
        death = self.temp_scalar * death
        return birth, death, 1 - (birth + death)

    def extend(self, additional_time):
        self.final_time += additional_time
        self.final_time_step = int(self.final_time / self.time_step_size)

    @classmethod
    def load_simulation(self, path_to_data) -> Self:
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
