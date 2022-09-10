# Implementation of the central schemes for networks of scalar conservation laws
Code that implements the schemes of the paper
> Central schemes for networked scalar conservation laws

by M. Herty, N. Kolbe, S. Müller

A preprint of the paper will be available soon.

Julia codes were written by N. Kolbe

## Contents
The repository contains a Julia project environment, which includes the package `CentralNetworkScheme.jl` implementating the schemes, and a directory (`experiments`) with scripts for various numerical experiments. The former implements schemes for 1-to-1 and 2-to-1 networks. Flux functions and boundary conditions can be freely chosen. Users can either use the first or the second order central scheme. For the LWR model on 2-to-1 networks the central coupling approach can be replaced by the flow maximization approach. 

Experiments on 1-to-1 networks for the Burgers' equation (`experiments/Burgers11.jl`) and on 2-to-1 networks for the LWR model (`experiments/LWR21.jl`) and the Buckley-Leverett equation (`experiments/BL21`) are included in the `experiments` folder.

## Usage 
Clone the repository on your local machine. To run the scripts activate the project environment, for example, by starting `julia` from the command line within the root folder of the repository, then changing to pkg mode typing `]` and afterwards running `activate .`. Then the experiments can be run by including the corresponding scripts, e.g., `include("experiments/LWR21.jl")`.  You can get started by going through the scripts and modifying them as you like. Scripts are commented and further documentation  will be added in the future. 

In case of questions or problems please contact the authors of the paper or file a GitHub issue in this repository.
