# Failure Response to Process Module in Cluster Tools with Residency Time Constraints Based on Two-Stage Mixed Integer Linear Programming

This repository provides some demos for [Failure Response to Process Module in Cluster Tools with Residency Time Constraints Based on Two-Stage Mixed Integer Linear Programming]().

## Requirements and Installation

- Matlab 2018b
- [YALMIP](https://yalmip.github.io/?n=Main.HomePage)

- [Gurobi Optimization](https://www.gurobi.com/)

## Getting Started

We release 3 files, named `cyclic_demo.m`, `FR_first_stage_demo.m` and `FR_second_stage_demo.m`, which represent MIP demos for cyclic scheduling MIP model, first-stage IP model and second-stage MILP model. We will give more details of each file in the following. You can just change the input to run your own experiments.

### Cyclic Scheduling

The meaning of each input in `cyclic_demo.m` is as the following:

- `robot_mode`: single or dual
- `n`: the number of PMs in the CT
- `w`: the robot loading or unloading time
- `v`: the robot rotating time
- `r`: $r_i$ denotes the processing time in ${\rm PM}_i$

- `residency`: the $i$-th element denotes the longest residency time allowed in ${\rm PM}_i$

- `K_PM`: wafer flow pattern (WFP) in the original paper

### Failure Response

#### First-stage Model

For the input in `FR_first_stage_demo.m`:

- `robot_mode`: single or dual
- `w`: the robot loading or unloading time
- `v`: the robot rotating time
- `WFP`: wafer flow pattern in the original paper

- `MT0_PM`: initial PM state
- `MTn_PM`: target PM state

### Second-stage Model

For the input in `FR_second_stage_demo.m`:

- `robot_mode`: single or dual
- `w`: the robot loading or unloading time
- `v`: the robot rotating time
- `WFP`: wafer flow pattern in the original paper
- `rho_raw`: the $i$-th element denotes the processing time in ${\rm PM}_i$
- `rho`: eliminate elements in `rho_raw` that correspond to PMs that have failed.
- `varthetaI`: initial processing time state
- `varthetaE`: target processing time state
- `MT0_PM`: initial PM state
- `MTn_PM`: target PM state
- `varpi_0`: initial robot state
- `varpi_1`: target robot state
- `delta`: longest wafer residency time vector
- `A_PM_idx`: used for extracting $A_{\rm PM}$ from matrix $A$
- `A_R_idx`: used for extracting `A_{\rm R}` from matrix $A$
- `Zl_idx`: used for extracting loading action from `Z`
- `Zu_idx`: used for extracting unloading action from `Z`
- `step`: number of tasks that the robot will perform $k$

## Contact
Please feel free to submit a Github issue if you have any questions or find any bugs. We do not guarantee any support, but will do our best if we can help.