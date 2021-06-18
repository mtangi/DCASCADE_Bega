# DCASCADE_Bega
The repository contains the code, input and output data used in the paper M. Tangi, S. Bizzi, K. Fryirs and A. Castelletti, "A dynamic, network scale sediment (dis)connectivity model to reconstruct historical sediment transfer and river reach sediment budgets"

---
The main folder contains the most important scripts and function used in the paper, as well as two workspaces containing the input data used and the main outputs:
- The workspace DCASCADE_Bega_input contains all the main input necessary to run D-CASCADE and reproduce the results shown in the paper;
- the workspace DCASCADE_Bega_results contains the main outputs of the D-CASCADE model, used to produce the figures and conclusions described in the paper;
- the script script_DCASCADEBega_main contains all the operations to simulate historic sediment (dis)connectivity trajectories on the Bega river network for the different hydrological scenarios;
- the script script_DCASCADEBega_futuresim contains the code used to generate future hydrological scenarios, as well as the future sediment (disc)connectivity trajectories in D-CASCADE;
- the script script_DCASCADEBega_figures contains all the code used to plot the figures shown in the paper;
- the function DCASCADEmodel_Bega is the main D-CASCADE function. It contains all the main step constituting the D-CASCADE model described in the paper, as well as the model add-ons. This version of the model is specifically designed for the Bega river network case study, and includes specifications of the model temporal steps and other simulations parameters described in the paper;
- the function DCASCADEmodel_Bega_future is an alternative version of the D-CASCADE model specifically tailored to simulate future sediment (dis)connectivity trajectories in the Bega river network To do so, it requires the outputs of a previous D-CASCADE historic simulation as a baseline;
- the folder D_CASCADE_functions contain all supporting functions and plot functions used in the code. Description of the functions purpose and functioning are provided as comments to the code.

---
The workspace DCASCADE_Bega_input contains:
- depth_hist, a matrix in each row reporting the year, discharge and water depth of the major flood events in the Bega river network. The discharge values refers to the Morans Crossing section, while water depth refers to the Bega Township gauging station. Discharge before 1940 is generated via depth/discharge correlation described in the supporting informations S1 in the paper. Discharge after 1940 is measured in the field. Not all flood events after 1940 are correlated to a specific discharge;
- n_Man is a txn matrix containing the Manning's roughness dynamic parameter for each reach n in each timestep t used in the historic simulations;
- Q_ARR_reach is a nx7 matrix reporting the discharge values for each reach relative the the flood classes with 1- 2- 5- 10- 20- 50- 100-year return period, estimated with the ARR model in Fryirs and Brierley (2001);
- Q_future_scenario is a 100x1 cell struct containing the 100 randomly-generated future discharge scenarios used in the future trajectories simulations;
- Q_Moran_18002020_flood contains the flood discharge referred to the Morans crossing gauging station from 1800 to 2020 (including 50 timesteps of Q1 for the initialisation) already reconstructed using height/discharge correlations;
- Q_Moran_4420 contains the daily discharge dataset from the Morans crossing gauging station from 1944 to 2020, used to derive flood frequency;
- reach_compact_group is a cell struct containing the IDs of the reaches which are compacted to display the results in Figure 5, 6 and 7 in the paper;
- ReachData is the main D-CASCADE input. It's a n x 1 Struct defining the features of the n network reaches. The network features refers to pre-ES conditions;
- SedDelRt_field contains the reach ID and the corresponding sediment delivery ratio as reported in Fryirs and Brierley (2001);
- Width_change is a 3xn matrix which contains the change in channel width in each reach in each of the 3 temporal phases where width change is expected, compared to the channel width reported in ReachData. Change can be both user-modified as in Phase 2 and driver by the channel expansion-add on as in Phase 1 and 3.

---
The workspace DCASCADE_Bega_results contains, together with the input data:
- data_plot and extended_outputs, which contains all the outputs of the historic simulation of the D-CASCADE model for the MixQ scenario, used as a baseline for the future simulations;
- data_plot_multiple contains all the major outputs of the historic D-CASCADE simulation for all 4 discharge scenarios. Each output is composed by a txn matrix, reporting the value of the output parameter for each reach n in each timestep t
- data_plot_multiple_compacted is a matrix structured identical to data_plot_multiple, only the outputs of the reaches included in reach_compact_group are compacted to better visualize the data, either by adding the value in each reach or averaging it.
- data_plot_scenarios and data_plot_scenarios_nores contains the results of the future hydrological scenarios, for the present-day vegetation and the no-restoration scenarios respectively;
- data_plot_scenarios_compacted and data_plot_scenarios_compacted_nores are structured identical to data_plot_scenarios, only the outputs of the reaches included in reach_compact_group are compacted to better visualize the data, either by adding the value in each reach or averaging it.
- SimVal contains the values of the validation parameters for the historical scenarios, together with the present-day and the pre-ES values;
- SimVal_future and SimVal_future_nores contains the values of the validation parameters for all the future scenarios, for the present-day vegetation and the no-restoration scenarios respectively.
