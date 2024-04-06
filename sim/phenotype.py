"""
A generic Phenotype class, and the PhenotypeInteractions class that describe how different governing structures can interact/ how similar we can claim they are.
All the PhenotypeStructure classes, describing how to build and mutate phenotypes.
"""

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


class SequencePhenotypeStructure(PhenotypeStructure):
    def __init__(self, sequences):
        self.sequences = tuple(sequences)
        self.ids = range(len(self.sequences))

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

    def get_mutation_weightings(self, phenotype):
        # TODO: Implement this properly
        return np.ones(len(self.sequences))

    def get_value(self, phenotype):
        """
        Get the value associated with a phenotype id.
        """
        return self.sequences[phenotype.id]  # By default, these are the same

    @classmethod
    def get_affinity_distance(self, phen_1, phen_2, range, binding_affinity_matrix):
        # For comparing different seqeuence phenotypes (i.e. tumour and T-cell), using binding affinities
        affinity = binding_affinity_matrix[phen_1.id, phen_2.id]
        distance = 1 / binding_affinity_matrix[phen_1.id, phen_2.id]
        # TODO: Improve this implementation; it's quite crude as it is and not calibrated
        print(distance)
        print(range)
        if distance < range:
            return 0
        else:
            return affinity

    @classmethod
    def get_sequence_distance(
        self, phen_1: Phenotype, phen_2: Phenotype, range, sequence_matrix
    ):
        # TODO: implement a sequence matrix
        print("henlo")
        return Levenshtein.distance(phen_1.get_value(), phen_2.get_value())

    @classmethod
    def get_interaction_function(self):
        """
        Return a function that handles interactions between two phenotypes with this structure of the form:
        `f(phen_1 : Phenotype, phen_2 : Phenotype, range : float, *data)`
        """
        return (
            self.get_sequence_distance
        )  # or is it better to do self.etc rather than the class method? Which works?

    @classmethod
    def get_standard_interaction_data(self):
        # An example of some data
        # Is this the best way to do this? In truth, we just need a specification of _what_ we need to provide
        sequence_matrix = np.identity(10)  # TODO: fix this
        return sequence_matrix

    @classmethod
    def get_cross_interaction_function(self):
        return self.get_affinity_distance

    @classmethod
    def get_standard_cross_interaction_data(self):
        # An example of some data
        # Is this the best way to do this? In truth, we just need a specification of _what_ we need to provide
        affinity_matrix = np.identity(10)  # TODO: fix this
        return affinity_matrix

    def __hash__(self):
        return hash(self.sequences)

    def __eq__(self, other):
        return self.sequences == other.sequences


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
        distance_type="line",
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

        struct = phenotype_1_.struct

        if distance <= range:
            return (
                self.TCR_binding_affinity
                * 1
                / (
                    min(phenotype_1 + range, struct.abs_max_value)
                    - max(phenotype_1 - range, -struct.abs_max_value)
                )
            )
        else:
            return 0

    def __hash__(self):
        return hash((self.abs_max_value, self.no_possible_values))

    def __eq__(self, other):
        return (self.abs_max_value == other.abs_max_value) and (
            self.no_possible_values == other.no_possible_values
        )

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

    def __init__(self, interaction_data: dict = {}, interaction_functions: dict = {}):
        self.interaction_data = interaction_data
        self.interaction_functions = interaction_functions
        self.phenotype_separation_scaling = {}

    def compute_interaction_scaling(self, phen_1: Phenotype, phen_2: Phenotype, range):
        struct_tuple = (phen_1.struct, phen_2.struct)
        data = self.interaction_data[struct_tuple]
        function = self.interaction_functions[struct_tuple]
        return function(phen_1, phen_2, range, data)

    def set_up_interactions(self, struct_tuple, data, function):
        """
        Assign an interaction function and its associated data to a particular pair of phenotypes.
        `function` is of the form `f(phen_1 : Phenotype, phen_2 : Phenotype, range : float, *data)`
        `data` is the parameter that is passed into the function.
        """
        self.interaction_data[struct_tuple] = data
        self.interaction_functions[struct_tuple] = function

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
            self.phenotype_separation_scaling[(phenotype_1, phenotype_2)] = (
                self.compute_interaction_scaling(phenotype_1, phenotype_2, range)
            )

        return self.phenotype_separation_scaling[(phenotype_1, phenotype_2)]

    @classmethod
    def get_default_sequence_interactions(
        self,
        CTL_struct: SequencePhenotypeStructure,
        tumour_struct: SequencePhenotypeStructure,
        sequence_matrix,
        affinity_matrix,
    ):
        """
        Create a PhenotypeInteractions object which handles interactions between two sequence-based phenotype structures.
        """
        ts = tumour_struct
        cs = CTL_struct
        phen_int = PhenotypeInteractions()

        phen_int.set_up_interactions(
            (ts, ts),
            sequence_matrix,
            SequencePhenotypeStructure.get_interaction_function(),
        )

        phen_int.set_up_interactions(
            (cs, cs),
            sequence_matrix,
            SequencePhenotypeStructure.get_interaction_function(),
        )

        phen_int.set_up_interactions(
            (ts, cs),
            affinity_matrix,
            SequencePhenotypeStructure.get_cross_interaction_function(),
        )

        phen_int.set_up_interactions(
            (cs, ts),
            np.transpose(affinity_matrix),
            SequencePhenotypeStructure.get_cross_interaction_function(),
        )

        return phen_int

    @classmethod
    def get_default_lattice_interactions(
        self, lattice_struct: LatticePhenotypeStructure
    ):
        """
        Create a PhenotypeInteractions object which handles interactions where all cells have the same lattice phenotype structure.
        """
        ls = lattice_struct
        phen_int = PhenotypeInteractions()

        phen_int.set_up_interactions(
            (ls, ls), ls.get_standard_interaction_data(), ls.get_interaction_function()
        )

        return phen_int
