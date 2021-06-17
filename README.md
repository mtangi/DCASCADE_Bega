# DCASCADE_Bega
The repository contains the code, input and output data used in the paper M. Tangi, S. Bizzi, K. Fryirs and A. Castelletti, "A dynamic, network scale sediment (dis)connectivity model to reconstruct historical sediment transfer and river reach sediment budgets"

The main folder contains the most important scripts and function used in the paper, as well as two workspaces containing the input data used and the main outputs:
- The workspace DCASCADE_Bega_input contains all the main input necessary to run D-CASCADE and reproduce the results shown in the paper;
- the workspace DCASCADE_Bega_results contains the main outputs of the D-CASCADE model, used to produce the figures and conclusions described in the paper;
- the script script_DCASCADEBega_main contains all the operations to simulate historic sediment (dis)connectivity trajectories on the Bega river network for the different hydrological scenarios;
- the script script_DCASCADEBega_futuresim contains the code used to generate future hydrological scenarios, as well as the future sediment (disc)connecitivity trajectories in D-CASCADE;
- the script script_DCASCADEBega_figures cointains all the code used to plot the figures shown in the paper;
- the function DCASCADEmodel_Bega is the main D-CASCADE function. It contains all the main step constituiting the D-CASCADE model described in the paper, as well as the model add-ons. This version of the model is specifically designed for the Bega river network case study, and includes spcifications of the model temporal steps and other simulations parameters described in the paper;
- the function DCASCADEmodel_Bega_future is an alternative version of the D-CASCADE model specifically tailored to simulate future sediment (dis)connectivity trajectories in the Bega river network To do so, it requires the outputs of a previous D-CASCADE historic simulation as a baseline;
- the folder D_CASCADE_functions contain all supporting functions and plot functions used in the code. Description of the functions purpose and functioning are provided as comments to the code.
