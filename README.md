# Asynchronous Sensitivity-Based Distributed NMPC

This repository contains the code accompanying the paper:

**“Asynchronous Sensitivity-Based Distributed NMPC”** by Maximilian Pierer von Esch, Andreas Völz, and Knut Graichen.

---

## Estimate Constants

To reproduce the numerical values and parameter estimates used in the paper:

```
cd EstimateConstants
```
and run the script estimate_constants.m to reproduce the numerical values and estimates in Section V.i) of the paper. 

## Distributed Execution 
This folder contains the C++ for the distributed (a)synchronous execution of SBDP within a DMPC framework for the example of three agents considered in the paper. 

```
cd DistributedExecution
```
## Dependencies
The code uses lcm (https://lcm-proj.github.io/lcm) and CMake. To install lcm on Ubuntu run:
```
sudo apt install liblcm-dev
```

Setup LCM Communication between raspberry pis
- Configure lcm to attach to certain url e.g. "udpm://239.255.76.67:7667"
- Add multicast route on eth0 via sudo ip route add 224.0.0.0/4 dev eth0 or sudo ip route replace 224.0.0.0/4 dev eth0 if a route already exists
- switch off wlan to avoid misrouting of data 
 - For an execution on a centralized PC remove the attachment of lcm to a certain adress and leave it blank as then the loopback adress is used

# Installation 
To build the project on Linux, clone this repository and run:
```

mkdir build
cd build
cmake ..
make 
```

This generates agent executables and a simulator executable in the folder build. 

First, run the agents in the terminal (on the individual Raspberry Pis) as
```
./agent <id> <mpc_iterations> <algorithm_iterations> <num_agents>
```
and then the simulator (on an additional Raspberry Pi) as 
```
./simulator <mpc_iterations> <initial_condition1> <initial_condition2> <initial_condition3>
```



