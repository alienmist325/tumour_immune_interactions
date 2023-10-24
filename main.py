from run import run
from cProfile import Profile
from pstats import SortKey, Stats

"""
with Profile() as profile:
    run()
    Stats(profile).strip_dirs().sort_stats(SortKey.CALLS).print_stats()
"""
run()
