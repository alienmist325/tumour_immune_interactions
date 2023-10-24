from run import run
from cProfile import Profile
from pstats import SortKey, Stats
from graphing import graph

"""
with Profile() as profile:
    run()
    Stats(profile).strip_dirs().sort_stats(SortKey.CALLS).print_stats()
"""
run()

if input("Would you like to see the result? \n") == "y":
    graph()
