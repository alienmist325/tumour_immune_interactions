from simulation import Phenotype, SequencePhenotypeStructure

struct = SequencePhenotypeStructure([])
p = Phenotype(struct, 2)
p2 = Phenotype(struct, 2)

a = {}


a[p] = 2
a[p2] = 4

print(a[p])
print(a[p2])
