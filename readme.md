# k-MST through MILP optimization

In this project, we provide different MILP formulations for the k-Minimum Spanning Tree Problem ([k-MST](https://en.wikipedia.org/wiki/Minimum_spanning_tree#k-minimum_spanning_tree)). In particular, we define
- a sequential formulation (inspired by Miller, Tucker and Zemlin, MTZ)
- a Single-Commodity Flow formulation, SCF
- a Multi-Commodity Flow formulation, MCF
- a Cycle Elimination Constraints formulation, CEC
- a Directed Cutset Constraints formulation, DCC

We model these formulations using ortools, and we compare their performance on a set of 10 graph instances.
We compare the licensed solver Gurobi 10.0.1 vs SCIP v803.
The last two formulations require the usage of a solver that supports Lazy Constraints. Only Gurobi (through the Python API `gurobipy`) was used for these formulations, even though SCIP could also have been explored. 

This project is part of the lecture [186.835 Mathematical Programming](https://tiss.tuwien.ac.at/course/courseDetails.xhtml?dswid=3722&dsrid=607&courseNr=186835&semester=2021S) at the Technical University of Vienna.

## Table of contents
- [Installation](#installation)
- [Configuration file](#configuration-file)
- [Execute](#execute)
- [Results](#results)

## Installation
Clone the repository and make sure that the folder structure is as follows:
```
mathematical_programming
├── kmst_instances
├──── g01.dat
├────  ....
├──── g10.dat
├── kmst_output
├── src_kmst
config.yaml
readme.md
requirements.txt
```

Install the required packages using pip or conda: 
```
pip install -r requirements.txt
conda install --file requirements.txt
```

The directory `problems` contains the LP or MILP solutions of some other simple problems: TSP, Minimum Steiner Tree, a League-Court scheduling problem, Bin Packing, Knapsack problem... Feel free to take a look :)

## Configuration file
The configuration file `config.yaml` is split into 3 sections:
- `basic`: contains the basic parameters of the program
- `single_execution`: contains the parameters for a single execution of the program
- `multiple_execution`: contains the list parameters for multiple executions of the program

The basic parameters are the following:
- `data_path`: path to the data folder
- `output_path`: path to the output folder
- `validate_solution`: can be [False, easy, hard]. If False, the solution is not validated. If easy, the solution is checked to equal the known optimal. If hard, the solution is validated checking graph properties: number of nodes selected, connectivity, is_tree...
- `instances`: list of the instances to be solved. Should be a list of integers. E.g.: if '3' is included, both '03_0' and '03_1' will be executed. 
- `time_limit`: time limit for the solver in seconds
- `execution_type`: can be [single, multiple]. If single, the program will execute the instances for a single configuration specified in the dictionary 'single_execution'. If multiple, the program will execute the instances trying all parameters specified in the dictionary 'multiple_execution'.

The parameters for a single execution and multiple executions are the same. In a single execution, just one parameter has to be specified. In a multiple execution, a list of parameters has to be specified. The parameters are the following:
- `solver`: can be [gurobi, scip]
- `formulation`: can be [MTZ, SCF, MCF, CEC, DCC]
- `hint_solution`: can be [True, False]. WARNING: only works for CEC and DCC. If True, a MIP start is provided.
- `tighten`: can be [True, False]. If True, more constraints are added to the minimal formulation. Does not work for CEC and DCC.
- `cuts`: can be [integral, fractional, both]. Only applies to CEC and DCC. If integral, cuts are applied to MIP solutions. If fractional, cuts are applied to LP solutions. If both, cuts are applied to both MIP and LP solutions.


## Execute
Just run the src_kmst/main.py file. The config is taken from the config.yaml file. The results will be stored in the output folder.
```
python src_kmst/main.py
```

## Results
The results are stored in the kmst_output folder.
The results are stored in a csv file named %Y%m%d-%H%M%S.csv. The columns are the following:
- `instance`: name of the instance
- `n`: number of nodes
- `m`: number of edges
- `k`: number of nodes to be selected
- `formulation`: formulation used
- `define_hints`: True if hints are defined, False otherwise
- `cuts`: cuts strategy applied
- `solver`: solver used
- `objective`: objective value
- `time`: time to formulate and solve the problem
- `nodes`: number of nodes explored
- `iterations`: number of simplex iterations
- `time_limit`: time limit for the solver in seconds
- `best_bound`: best bound found by the solver
- `status`: status of the solver
- `theoretical_optimal`: theoretical optimal value
- `opt_gap`: optimality gap
- `tighten`: True if the formulation is tightened, False otherwise
- `solve_time`: time spent in the solve method

