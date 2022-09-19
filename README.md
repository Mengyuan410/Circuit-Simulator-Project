# Circuit-Simulator-Project
:star:This repository is the end of year project of 1st Year Imperial College London EIE department. 

The project is selected as the best Circuit Simulator Project (overall Runner-Up Award) among the cohort. 

## Overview
The circuit simulator in the repository written in C++ simulates the basic functions of LTSpice. It can perform AC and DC analysis on the given circuit. The circuit should be provided to the circuit simulator in a netlist format (samples are listed in */test* directory). The project requirement is specified in , and the project report is in *Report.pdf*, where you can find the logic and techniques used for the circuit simulator. The circuit simulator will find the transfer functions of a selected node and output a netlist (i.e. simdata.txt) containing the magnitude and phase of that node under different frequencies. Then MATLAB code *plotsim.m* is used to plot a magnitude-frequency graph and a phase-frequency graph.

## Running Guide
Use code <code> ./run.sh </code> to run the circuit simulator. You will notice that simdata.txt will be updated by the information of the new circuit.
Then run *plotsim.m* using MATLAB to obtain the magnitude-frequency and phase-frequency graphs.
