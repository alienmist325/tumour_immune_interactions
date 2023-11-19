from run import run
from cProfile import Profile
from pstats import SortKey, Stats
from graphing import graph_from_path, savefig

"""
with Profile() as profile:
    run()
    Stats(profile).strip_dirs().sort_stats(SortKey.CALLS).print_stats()
"""
run()

if input("Would you like to save the graph?") == "y":
    savefig()
