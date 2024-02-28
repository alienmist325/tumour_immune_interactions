"""
A generic Phenotype class, and the PhenotypeInteractions class that describe how different governing structures can interact/ how similar we can claim they are.
"""

import numpy as np
from phenotype_structures import PhenotypeStructure, SequencePhenotypeStructure, LatticePhenotypeStructure


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
            self.phenotype_separation_scaling[
                (phenotype_1, phenotype_2)
            ] = self.compute_interaction_scaling(phenotype_1, phenotype_2, range)

        return self.phenotype_separation_scaling[(phenotype_1, phenotype_2)]

    @classmethod
    def get_default_sequence_interactions(
        self,
        tumour_struct: SequencePhenotypeStructure,
        CTL_struct: SequencePhenotypeStructure,
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
