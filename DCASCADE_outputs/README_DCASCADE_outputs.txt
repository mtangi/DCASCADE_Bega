The folder DCASCADE_outputs contains the folder with all the main outputs for the Bega case study, and the script to plot the model's results:
- the script script_DCASCADEBega_figures contains all the code used to plot the figures shown in the paper;
- the folder DCASCADE_Bega_results contains all the main outputs of the D-CASCADE model, used to produce the figures and conclusions described in the paper.

---

The folder DCASCADE_Bega_results contains several workspaces with all the main outputs for the Bega case divided:
- data_plot_Mix contains data_plot and extended_outputs, with all the outputs of the historic simulation of the D-CASCADE model for the MixQ scenario, used as a baseline for the future simulations;
- data_plot_multiple contains all the major outputs of the historic D-CASCADE simulation for all 4 discharge scenarios. Each output is composed by a txn matrix, reporting the value of the output parameter for each reach n in each timestep t
- data_plot_multiple_compacted is a matrix structured identical to data_plot_multiple, only the outputs of the reaches included in reach_compact_group are compacted to better visualize the data, either by adding the value in each reach or averaging it.
- data_plot_scenarios and data_plot_scenarios_nores contains the results of the future hydrological scenarios, for the present-day vegetation and the no-restoration scenarios respectively;
- data_plot_scenarios_compacted and data_plot_scenarios_compacted_nores are structured identical to data_plot_scenarios, only the outputs of the reaches included in reach_compact_group are compacted to better visualize the data, either by adding the value in each reach or averaging it.
- SimVal contains all the validation matrixes: SimVal contains the values of the validation parameters for the historical scenarios, together with the present-day and the pre-ES values, SimVal SimVal_future and SimVal_future_nores contain the values of the validation parameters for all the future scenarios, for the present-day vegetation and the no-restoration scenarios respectively;
