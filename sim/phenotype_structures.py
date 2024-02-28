"""
All the PhenotypeStructure classes, describing how to build and mutate phenotypes.
"""


import random
from abc import ABC, abstractmethod
import Levenshtein
from phenotype import Phenotype
import numpy as np

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


class SequencePhenotypeStructure(PhenotypeStructure):
    def __init__(self, sequences):
        self.sequences = sequences

    def get_random_phenotype(self):
        """
        Generate a valid phenotype id.
        """
        return random.choice(self.sequences)

    def get_random_mutation(self, phenotype):
        """
        Get a possible mutation of the phenotype.
        """
        return random.choice(self.sequences, self.get_mutation_weightings())
    
    def get_mutation_weightings(self, phenotype):
        return np.ones(len(self.sequences))

    def get_value(self, phenotype):
        """
        Get the value associated with a phenotype id.
        """
        return phenotype.id  # By default, these are the same

    @classmethod
    def get_affinity_distance(self, phen_1, phen_2, binding_affinity_matrix):
        # For comparing different seqeuence phenotypes (i.e. tumour and T-cell), using binding affinities
        pass

    @classmethod
    def get_sequence_distance(self, phen_1: Phenotype, phen_2: Phenotype, sequence_matrix):
        # TODO: implement a sequence matrix
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
        return (
            self.compute_interaction_scaling
        ) 

    @classmethod
    def get_standard_interaction_data(self):
        # An example of some data
        # Is this the best way to do this? In truth, we just need a specification of _what_ we need to provide
        distance_type = "line"
        return distance_type

    @classmethod
    def compute_interaction_scaling(
        self, phenotype_1_ : Phenotype, phenotype_2_ : Phenotype, range, distance_type="line"
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
            return 1 / (
                min(phenotype_1 + range, struct.abs_max_value)
                - max(phenotype_1 - range, -struct.abs_max_value)
            )
        else:
            return 0

    def __hash__(self):
        return hash((self.abs_max_value, self.no_possible_values))
    
    def __eq__(self, other):
        return (self.abs_max_value == other.abs_max_value) and (self.no_possible_values == other.no_possible_values)

    """
    def get_random_phenotype(self):
        return (
            random.randint(0, self.no_possible_values - 1) * self.step_size
        ) - self.abs_max_value
    """

