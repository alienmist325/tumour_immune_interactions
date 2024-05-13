"""
A generic Phenotype class, and the PhenotypeInteractions class that describe how different governing structures can interact/ how similar we can claim they are.
All the PhenotypeStructure classes, describing how to build and mutate phenotypes.
"""

from dataclasses import dataclass
import numpy as np
import random
from abc import ABC, abstractmethod
import Levenshtein


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

    @abstractmethod
    def __hash__(self):
        pass

    @abstractmethod
    def __eq__(self, other):
        pass

    # (Not anymore) IDEA: We want to have a general function here which will take in the two ids and get the scaling. If they're the same type, we should redirect back down to the
    # relevant class, since e.g. LatticePhenotypeStructure _can definitely_ compare between two lattice phenotypes. If they're different, we'll cover it here and overload.


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

    def __str__(self):
        return f"{self.id} inside {str(self.struct)}"


@dataclass
class SequencePhenotypeInteractionData:
    interaction_matrix: object
    scaling: float


class SequencePhenotypeStructure(PhenotypeStructure):
    def __init__(self, sequences):
        self.sequences = tuple(sequences)
        self.ids = range(len(self.sequences))
        self.distance_matrix = None

    def get_random_phenotype(self):
        """
        Generate a valid phenotype id. ==Or an actual phenotype?==
        """
        return Phenotype(self, random.choice(self.ids))

    def get_random_mutation(self, phenotype):
        """
        Get a possible mutation of the phenotype.
        """
        return Phenotype(
            self,
            random.choices(self.ids, self.get_mutation_weightings(phenotype))[0],
        )

    def get_mutation_weightings(self, phenotype: Phenotype):
        # TODO: Implement this properly
        weightings = np.ones(len(self.sequences))
        weightings[phenotype.id] = (
            0  # So that you cannot become the same phenotype you began as
        )
        return weightings

    def get_value(self, phenotype):
        """
        Get the value associated with a phenotype id.
        """
        return self.sequences[phenotype.id]  # By default, these are the same

    @classmethod
    def get_affinity_scaling(
        self, phen_1, phen_2, range, data: SequencePhenotypeInteractionData
    ):
        # For comparing different seqeuence phenotypes (i.e. tumour and T-cell), using binding affinities
        binding_affinity_matrix = data.interaction_matrix
        binding_scaling = data.scaling
        affinity = binding_affinity_matrix[phen_1.id, phen_2.id]
        # distance = 1 / binding_affinity_matrix[phen_1.id, phen_2.id]
        # TODO: Improve this implementation; it's quite crude as it is and not calibrated
        # print(distance)
        # print(range)
        return binding_scaling * affinity
        # if distance < range:
        #    return 0
        # else:
        #    return affinity

    @classmethod
    def get_sequence_scaling(
        self,
        phen_1: Phenotype,
        phen_2: Phenotype,
        range,
        data: SequencePhenotypeInteractionData,
        distance_scaling=0.1,
    ):

        sequence_scaling = data.scaling  # 0.1 by default
        return sequence_scaling / (
            1
            + distance_scaling
            * SequencePhenotypeStructure.get_sequence_distance(
                phen_1, phen_2, range, data
            )
        )

    @classmethod
    def get_sequence_distance(
        self,
        phen_1: Phenotype,
        phen_2: Phenotype,
        range,
        data: SequencePhenotypeInteractionData,
    ):
        # TODO: implement a sequence matrix
        sequence_matrix = data.interaction_matrix
        # TODO: should the scaling be here
        return Levenshtein.distance(phen_1.get_value(), phen_2.get_value())

    @classmethod
    def get_interaction_function(self):
        """
        Return a function that handles interactions between two phenotypes with this structure of the form:
        `f(phen_1 : Phenotype, phen_2 : Phenotype, range : float, *data)`
        """
        return (
            self.get_sequence_scaling
        )  # or is it better to do self.etc rather than the class method? Which works?

    @classmethod
    def get_standard_interaction_data(self):
        """
        An unused function that illustrates usage.
        """
        # An example of some data
        # Is this the best way to do this? In truth, we just need a specification of _what_ we need to provide
        sequence_matrix = np.identity(10)  # TODO: fix this
        return sequence_matrix

    @classmethod
    def get_cross_interaction_function(self):
        return self.get_affinity_scaling

    @classmethod
    def get_standard_cross_interaction_data(self):
        """
        An unused function that illustrates usage.
        """
        # An example of some data
        # Is this the best way to do this? In truth, we just need a specification of _what_ we need to provide
        affinity_matrix = np.identity(10)  # TODO: fix this
        return affinity_matrix

    def __hash__(self):
        return hash(self.sequences)

    def __eq__(self, other):
        return self.sequences == other.sequences

    def get_distance_matrix(self, data: SequencePhenotypeInteractionData):
        if self.distance_matrix is None:

            self.compute_distance_matrix(data)

        return self.distance_matrix

    def compute_distance_matrix(self, data: SequencePhenotypeInteractionData):
        # Create a list of all possible phenotypes
        print("computing")
        all_phenotypes = [Phenotype(self, id) for id in self.ids]
        num = len(all_phenotypes)
        all_phenotype_pair_matrix = np.transpose(
            np.meshgrid(all_phenotypes, all_phenotypes), (2, 1, 0)
        )
        all_phenotype_pair_tuple_matrix = np.array(np.ones((num, num)), dtype=object)

        for idx, row in enumerate(all_phenotype_pair_matrix):
            for idy, pair in enumerate(row):
                all_phenotype_pair_tuple_matrix[idx][idy] = tuple(pair)

        def get_sequence_distance_by_pair(pair: tuple):
            return SequencePhenotypeStructure.get_sequence_distance(
                pair[0], pair[1], 0, data
            )

        vec_get_sequence_distance = np.vectorize(get_sequence_distance_by_pair)

        distance_matrix = vec_get_sequence_distance(all_phenotype_pair_tuple_matrix)

        self.distance_matrix = distance_matrix

    @classmethod
    def sequence_matrix_adjustment(self, sequence_matrix):
        """
        Scale down the weights in the sequence matrix so all values are in [0,1]. Assumes all values are positive.
        """
        sequence_matrix = np.array(sequence_matrix)
        max = sequence_matrix.max()
        return sequence_matrix / max


@dataclass
class LatticeInteractionData:
    distance_type: str


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
        self.ids = range(no_possible_values)
        self.step_size = 2 * abs_max_value / no_possible_values  # check this
        self.no_possible_values = no_possible_values  # For an id
        self.distance_matrix = None

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
        return self.compute_interaction_scaling

    @classmethod
    def get_standard_interaction_data(self):
        """
        TODO: Remove as defunct
        """
        # An example of some data
        # Is this the best way to do this? In truth, we just need a specification of _what_ we need to provide
        distance_type = "line"
        return distance_type

    @classmethod
    def compute_interaction_scaling(
        self,
        phenotype_1_: Phenotype,
        phenotype_2_: Phenotype,
        range,
        data: LatticeInteractionData,
    ):
        """
        Return 0 if phenotypes are out of the selectivity range; otherwise return a constant, modified to account for the boundary.
        """
        phenotype_1 = phenotype_1_.get_value()
        phenotype_2 = phenotype_2_.get_value()
        if data.distance_type == "line":
            distance = abs(phenotype_1 - phenotype_2)
        if data.distance_type == "circular":
            distance = LatticePhenotypeStructure.get_circular_distance(
                phenotype_1, phenotype_2
            )

        struct = phenotype_1_.struct

        if distance <= range:
            return 1 / (
                min(phenotype_1 + range, struct.abs_max_value)
                - max(phenotype_1 - range, -struct.abs_max_value)
            )
        else:
            return 0

    def __hash__(self):
        return hash((self.abs_max_value, self.no_possible_values))

    def __eq__(self, other):
        return (self.abs_max_value == other.abs_max_value) and (
            self.no_possible_values == other.no_possible_values
        )

    def get_distance_matrix(self, data):
        """
        TODO: Remove or make more general, since this was copied from `SequencePhenotypeStructure`.
        This doesn't work fully correctly.
        """
        if self.distance_matrix is None:

            self.compute_distance_matrix(data)

        return self.distance_matrix

    def compute_distance_matrix(self, data):
        """
        TODO: Remove or make more general, since this was copied from `SequencePhenotypeStructure`
        This certainly doesn't work fully correctly.
        """
        # Create a list of all possible phenotypes
        print("computing")
        all_phenotypes = [Phenotype(self, id) for id in self.ids]
        num = len(all_phenotypes)
        all_phenotype_pair_matrix = np.transpose(
            np.meshgrid(all_phenotypes, all_phenotypes), (2, 1, 0)
        )
        all_phenotype_pair_tuple_matrix = np.array(np.ones((num, num)), dtype=object)

        for idx, row in enumerate(all_phenotype_pair_matrix):
            for idy, pair in enumerate(row):
                all_phenotype_pair_tuple_matrix[idx][idy] = tuple(pair)

        def get_sequence_distance_by_pair(pair: tuple):
            return self.compute_interaction_scaling(pair[0], pair[1], 0, data)

        vec_get_sequence_distance = np.vectorize(get_sequence_distance_by_pair)

        distance_matrix = vec_get_sequence_distance(all_phenotype_pair_tuple_matrix)

        self.distance_matrix = distance_matrix

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

    def __init__(
        self,
        interaction_data: dict = {},
        interaction_functions: dict = {},
        cross_term_multipliers: dict = {},
    ):
        self.interaction_data = interaction_data
        self.interaction_functions = interaction_functions
        self.cross_term_multiplers = cross_term_multipliers
        self.phenotype_separation_scaling = {}

    def compute_interaction_scaling(self, phen_1: Phenotype, phen_2: Phenotype, range):
        struct_tuple = (phen_1.struct, phen_2.struct)
        data = self.interaction_data[struct_tuple]
        function = self.interaction_functions[struct_tuple]
        return function(phen_1, phen_2, range, data)

    def set_up_interactions(
        self, struct_tuple, data, function, cross_term_multiplier=1
    ):
        """
        Assign an interaction function and its associated data to a particular pair of phenotypes.
        `function` is of the form `f(phen_1 : Phenotype, phen_2 : Phenotype, range : float, *data)`
        `data` is the parameter that is passed into the function.
        `cross_term_multiplier` is a scaling that is applied only when phenotypes of different CellBundles interact
        """
        self.interaction_data[struct_tuple] = data
        self.interaction_functions[struct_tuple] = function
        self.cross_term_multiplers[struct_tuple] = cross_term_multiplier

    def get_interaction_scaling(
        self,
        phenotype_1,
        phenotype_2,
        range,
    ):
        if (
            phenotype_1.id,
            phenotype_2.id,
        ) not in self.phenotype_separation_scaling:
            self.phenotype_separation_scaling[(phenotype_1.id, phenotype_2.id)] = (
                self.compute_interaction_scaling(phenotype_1, phenotype_2, range)
            )

        return self.phenotype_separation_scaling[(phenotype_1.id, phenotype_2.id)]

    def get_cross_term_multipler(
        self, phen_1_struct: PhenotypeStructure, phen_2_struct: PhenotypeStructure
    ):
        struct_tuple = (phen_1_struct, phen_2_struct)
        return self.cross_term_multiplers[struct_tuple]

    @classmethod
    def get_default_sequence_interactions(
        self,
        CTL_struct: SequencePhenotypeStructure,
        tumour_struct: SequencePhenotypeStructure,
        sequence_matrix,
        affinity_matrix,
        binding_scaling,
    ):
        """
        Create a PhenotypeInteractions object which handles interactions between two sequence-based phenotype structures.
        """
        ts = tumour_struct
        cs = CTL_struct
        phen_int = PhenotypeInteractions()

        sequence_matrix = SequencePhenotypeStructure.sequence_matrix_adjustment(
            sequence_matrix
        )
        homogeneous_data = SequencePhenotypeInteractionData(sequence_matrix, 0.1)

        phen_int.set_up_interactions(
            (ts, ts),
            homogeneous_data,
            SequencePhenotypeStructure.get_interaction_function(),
        )

        phen_int.set_up_interactions(
            (cs, cs),
            homogeneous_data,
            SequencePhenotypeStructure.get_interaction_function(),
        )

        cross_data = SequencePhenotypeInteractionData(affinity_matrix, binding_scaling)

        phen_int.set_up_interactions(
            (ts, cs),
            cross_data,
            SequencePhenotypeStructure.get_cross_interaction_function(),
        )

        transposed_cross_data = SequencePhenotypeInteractionData(
            np.transpose(affinity_matrix), binding_scaling
        )

        phen_int.set_up_interactions(
            (cs, ts),
            transposed_cross_data,
            SequencePhenotypeStructure.get_cross_interaction_function(),
        )

        tumour_struct.compute_distance_matrix(homogeneous_data)
        CTL_struct.compute_distance_matrix(homogeneous_data)

        return phen_int

    @classmethod
    def get_default_lattice_interactions(
        self,
        lattice_struct: LatticePhenotypeStructure,
        TCR_binding_affinity,
        distance_type="line",
    ):
        """
        Create a PhenotypeInteractions object which handles interactions where all cells have the same lattice phenotype structure.
        """
        ls = lattice_struct
        phen_int = PhenotypeInteractions()
        data = LatticeInteractionData(distance_type)

        phen_int.set_up_interactions(
            (ls, ls), data, ls.get_interaction_function(), TCR_binding_affinity
        )

        return phen_int
