
%% define future scenarios length and number

future_time = 100;
n_scenarios = 100;

%% generate future Q using Q with different percentiles

reach_type = [ReachData.reach_type]; %extract the reach type from ReachData

swamp_example = find(reach_type == 4,1);
other_example = find(reach_type ~= 4,1);

%calculate frequencies for each RP class based on historic data from Morans
%Crossing gaugin station
%NB since we are simulating for the MixQ scenario, the frequencies will be
%different for type 4 reaches

% matrix Q_Moran_18002020_flood_prc can be obtained in the Q_scenario
% identification section in DCASCADE_Begamainscript
prc_values = [1 2 5 10 20 50 100];
boundaries = [ [[2 5 10 20 50 100] - 0.00001] ; [[1 2 5 10 20 50]+ 0.00001]];

prob_prc_swamp = sum(Q_Moran_18002020_flood_prc(50:end,swamp_example) == prc_values,1)/sum(Q_Moran_18002020_flood_prc(50:end,swamp_example) == prc_values,'all');
prob_prc_other = sum(Q_Moran_18002020_flood_prc(50:end,other_example) == prc_values,1)/sum(Q_Moran_18002020_flood_prc(50:end,other_example) == prc_values,'all');
%initialize cell of future scenarios

% define RP classes, boundaries between classes for flood classification
% and a vector defining if a reach is type 4 or others
reach_sel = ones(1,length(reach_type)); reach_sel(reach_type ==4) = 2;

% initialize future discharge cell structure
Q_future_scenarios = cell(1,n_scenarios);

%loop for all scenarios
for s=1:n_scenarios
    
    Q_future_scenarios{s} = zeros(future_time,length(reach_type));
        
    Q_Moran_future_flood_prc = zeros(future_time,size(Q_Moran_18002020_flood_prc,2));

    for i=1:future_time

        %for each future timestep, it extracts random RP flood event and
        %classify it according to hist probability. 
        %Since we are simulating for the MixQ scenario, the
        %classification will be different for type 4 reaches
        r = rand;
        prc_type4 = find( [r > cumsum(prob_prc_swamp)] ,1, 'last') + 1;
        prc_other = find( [r > cumsum(prob_prc_other)] ,1, 'last') + 1;
        
        if isempty(prc_type4); prc_type4 = 1; end
        if isempty(prc_other); prc_other = 1; end
        
        % The randomly generated RP is attributed to each reach for the
        % timestep. Type 4 reaches may receive a different classification
        Q_Moran_future_flood_prc(i,reach_type == 4) = repmat(prc_values(prc_type4),[1 sum(reach_type == 4)] );
        Q_Moran_future_flood_prc(i,reach_type ~= 4) = repmat(prc_values(prc_other),[1 sum(reach_type ~= 4)] );
        
        idx_prc = find(Q_Moran_future_flood_prc(i,1) < prc_values,1)-1;
        if isempty(idx_prc); idx_prc= length(prc_values)-1; end    

        % Convert the RP class with the actual discharge obtained from ARR
        % model
        for n = 1: length(reach_type)
            
            if Q_Moran_future_flood_prc(i,n) < boundaries(reach_sel(n), idx_prc)
                Q_Moran_future_flood_prc(i,n) = prc_values(idx_prc);
                Q_future_scenarios{s}(i,n) =  Q_ARR_reach(n,idx_prc);

            else
                Q_Moran_future_flood_prc(i,n) = prc_values(idx_prc+1);
                Q_future_scenarios{s}(i,n) =  Q_ARR_reach(n,idx_prc+1);

            end
         end

    end


end

clear Q_Moran_future_flood_prc idx_prc prc_swamp prc_other r prc_values  swamp_example other_example

%% Rerun the MixQ scenario
%Important : To generate the historic baseline for the future simulation,
%we need all the outputs of the D-CASCADE model for the historic
%simulations with the MixQ scenario. These are not just in the data_plot
%matrix, but also  in the extended_output retuened in the extended_output
%cell struct by function DCASCADEmodel_Bega;  

id_sim = 4;
Q = Q_scenarios{id_sim};

% Find river network struct, containing 
Network = graph_preprocessing(ReachData);

%run D-CASCADE
[data_plot,extended_output] = DCASCADEmodel_Bega( ReachData , Network , Q , Q_ARR_reach, n_Man , Width_change );

%% run CASCADE for the future scenarios
%Important : To generate the historic baseline for the future simulation,
%we need all the outputs of the D-CASCADE model for the historic
%simulations with the MixQ scenario. These are not just the data_plot
%matrix, but all storage matrix defined in the function DCASCADEmodel_Bega;

ID_scenario = 1; % ID_scenario defines which future scenario we are simulating, either the present-day vegetation scenario (ID_scenario = 1) or the pre-recovery one

if ID_scenario == 1
    [data_plot_scenarios] = DCASCADEmodel_Bega_future( ReachData , Network , Q_future_scenarios , Q_ARR_reach , data_plot, extended_output , ID_scenario );
else 
    [data_plot_scenarios_nores] = DCASCADEmodel_Bega_future( ReachData , Network , Q_future_scenarios , Q_ARR_reach , data_plot, extended_output , ID_scenario );
end

%% compact the D-CASCADE future outputs for the reaches in the lowland sections displayed in the figures
% replace data_plot_scenarios with data_plot_scenarios_nores to obtain data
% for the pre-restoration scenario

%we named data_plot_scenarios_compacted_nores the struct containing the
%compacted outputs for the pre-restoration scenario  
data_plot_scenarios_compacted = data_plot_scenarios;

for s=1:size(data_plot_scenarios_compacted,1)
    for i=1:length(reach_compact_group)
    reaches_compact = reach_compact_group{i};

        for j=1:9

            if any(j==[2,3])
                data_plot_scenarios_compacted{s,1}{j,2}(:,reaches_compact) = repmat(sum(data_plot_scenarios_compacted{s,1}{j,2}(:, reaches_compact),2) , [1, length(reaches_compact)] ) ;
            else
                data_plot_scenarios_compacted{s,1}{j,2}(:,reaches_compact) = repmat(mean(data_plot_scenarios_compacted{s,1}{j,2}(:, reaches_compact),2) , [1, length(reaches_compact)] ) ;       
            end

        end
    data_plot_scenarios_compacted{s,1}{10,2}{1}(:,reaches_compact) = repmat(sum(data_plot_scenarios_compacted{s,1}{10,2}{2}(:, reaches_compact),2) , [1, length(reaches_compact)] ) ;
    data_plot_scenarios_compacted{s,1}{10,2}{2}(:,reaches_compact) = repmat(sum(data_plot_scenarios_compacted{s,1}{10,2}{2}(:, reaches_compact),2) , [1, length(reaches_compact)] ) ;
    end
end
   
%% define validation struct 
% replace data_plot_scenarios with data_plot_scenarios_nores in data_val to
% obtain the validation parameters for the pre-restoration scenario

data_val = data_plot_scenarios;
reach_type = [ReachData.reach_type]; %extract the reach type from ReachData

total_dep = [ReachData.deposit].*[ReachData.Length];

limt_future = 100;

%we named SimVal_future_nores the validation struct for the pre-restoration
%scenario 
SimVal_future = struct('simname',num2cell(zeros(1,size(data_val,1))));

[SimVal_future(1:size(data_val,1)).simname] = data_val{:,2} ;

% SedDelRt
% the matrix SedDelRt_field contains all the available field measurements for
% the sediment delivery ratio 

for i=1:length(SedDelRt_field)
    for j = 1:size(data_val,1)
        SimVal_future(j).([ 'SedDelRt_' num2str(SedDelRt_field(i,1))]) = data_val{j,1} {9,2} (limt_future,SedDelRt_field(i,1));
    end
end

% Type 4 - % incised
prc_eroded_sw = num2cell(1 - cell2mat(cellfun( @(x)sum(x{3,2}( limt_future , reach_type == 4  )) , data_val(:,1),'UniformOutput',0)) ./ sum(total_dep(reach_type == 4)) )     ;

[SimVal_future(1:size(data_val,1)).prc_eroded_sw] = prc_eroded_sw{:};

% Type 3 - Width
Wac_ll =  cellfun( @(x)mean(x{1,2}( limt_future , 94:99  )) , data_val(:,1),'UniformOutput',0 );
[SimVal_future(1:size(data_val,1)).Wac_ll] = Wac_ll{:} ;

% Type 3 - D50
D50_ll =  cellfun( @(x)mean(x{5,2}( limt_future-6:limt_future , 94:99  ),'all') , data_val(:,1),'UniformOutput',0 );
[SimVal_future(1:size(data_val,1)).D50_ll] = D50_ll{:} ;

clear  data_val
