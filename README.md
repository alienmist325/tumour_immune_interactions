A set of simulations modelling the behaviour and interaction of T-cells and tumour cells in the body.


## Setting up

Begin by installing the `tumour_immune_interactions` package. 

```py
python -m pip install -e .
```

You will then for example be able to do:

```py
from sim.discrete_model import Simulation
```

## Configuring

### First create or find the configuration you would like within `configurations.csv`
See `sim/config/readme.md`.

### The primary simulation configuration file `conf.py`

Item | Effect
--- | ---
interrupt | Set to true to safely interrupt at the next time step. This will save the simulation as a `.pickle`, and end the execution. (The file will also be reset so `interrupt = False` once more)
debug | Set to true to see output triggered by `self.print` inside `Simulation`
path_to_data | Where you read from and save to for simulation `sim.pickle` files
path_to_output | Where graph outputs and more go.
sim_state_init_type | `detailed` or ==?==. Specifies how much data is stored inside each `SimulationState`. `detailed` will be more performance heavy, but is certainly necessary for e.g. fish plotting

### Caveats

You can specify a `sequence_matrix_config` and `Simulation` will load this in, but the simulation logic does not currently support using it (currently it remains unused and an "identity" matrix is used instead).

## Running a simulation

Running `main.py` will run the discrete simulation. 
Running `continuous_run.py` will run the continuous simulation (not developed by me, but adapted by me)