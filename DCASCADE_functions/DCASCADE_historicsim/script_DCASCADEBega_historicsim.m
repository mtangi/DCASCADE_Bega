% this script contaions the operations necessary to run DCASCADE simulation
% on the Bega river network, for the 

%% create Q scenarios using Q_ARR

% Q_Moran_4420 contains the daily Q from the morans crossing gauging
% station from 1944 to 2020, used to derive flood frequency

% Q_Moran_18002020_flood contains the flood Q from the morans crossing
% gauging station from 1800 to 2020 (including 50 timesteps of Q1 for the
% initialization) already reconstructed using heigth/discharge correlations

%id_sim
%1 - highQ
%2 - medQ
%3 - lowQ
%4 - mixQ

reach_type = [ReachData.reach_type]; %extract the reach type from ReachData

%initialize matrix
Q_scenarios = cell(1,4);

prc_values = [1 2 5 10 20 50 100];
range = [1:0.1:100];

prc = prctile(Q_Moran_4420,[1 - 76.91./prc_values/length(Q_Moran_4420)]*100);
prc_full = prctile(Q_Moran_4420,[1 - 76.91./range/length(Q_Moran_4420)]*100);

for id_sim = 1:4
    if id_sim == 4 % MixQ scenario

        Q = zeros(length(Q_Moran_18002020_flood), size(Q_ARR_reach,1));

        boundaries = [ [[2 5 10 20 50 100] - 0.00001] ; [[1 2 5 10 20 50]+ 0.00001]];
        reach_sel = ones(1,length(reach_type));
        reach_sel(reach_type ==4) = 2;

        %Q_Moran_18002020_flood_prc contains the flood class for each flood
        %event for each reach
        Q_Moran_18002020_flood_prc = zeros(length(Q_Moran_18002020_flood),length(reach_type));

        for i=1:size(Q_Moran_18002020_flood,1)

            idx_prcF = find(Q_Moran_18002020_flood(i) < prc_full,1);
            if isempty(idx_prcF); idx_prcF= length(prc_full); end

            Q_Moran_18002020_flood_prc(i,:) = repmat(range(idx_prcF), [1,length(reach_type)]) ;

            idx_prc = find(Q_Moran_18002020_flood_prc(i,1) < prc_values,1)-1;
            if isempty(idx_prc); idx_prc= length(prc_values)-1; end    

            for n= 1: length(reach_type)
                if Q_Moran_18002020_flood_prc(i,n) < boundaries(reach_sel(n), idx_prc)
                    Q_Moran_18002020_flood_prc(i,n) = prc_values(idx_prc);

                    Q(i,n) =  Q_ARR_reach(n,idx_prc);

                else
                    Q_Moran_18002020_flood_prc(i,n) = prc_values(idx_prc+1);
                    Q(i,n) =  Q_ARR_reach(n,idx_prc+1);

                end
            end
        end

        Q_scenarios{id_sim} = Q;
    else

        switch id_sim
            case 1
                boundaries = [1 2 5 10 20 50 100]+ 0.00001; %high simulation
            case 2 
                boundaries = [1.5 3.5 7.5 15 35 75]; %medium simulation
            case 3   
                boundaries = [2 5 10 20 50 100 100] - 0.00001; %lower simulation
        end
        Q = zeros(length(Q_Moran_18002020_flood), size(Q_ARR_reach,1));

        Q_Moran_18002020_flood_prc = zeros(length(Q_Moran_18002020_flood),1);

        for i=1:size(Q_Moran_18002020_flood,1)

            idx_prcF = find(Q_Moran_18002020_flood(i) < prc_full,1);
            if isempty(idx_prcF); idx_prcF= length(prc_full); end

            Q_Moran_18002020_flood_prc(i) = range(idx_prcF);

            idx_prc = find(Q_Moran_18002020_flood_prc(i) < prc_values,1)-1;
             if isempty(idx_prc); idx_prc= length(prc_values)-1; end

            % approximate the percentile of the flow to one of the known percentiles defined in prc_values
            if Q_Moran_18002020_flood_prc(i) < boundaries(idx_prc)
                Q_Moran_18002020_flood_prc(i) = prc_values(idx_prc);
            else
                Q_Moran_18002020_flood_prc(i) = prc_values(idx_prc+1);
            end

            % calculate Q by doing liner interpolation of the flow betweeen the
            % prc_values that bound the Q_Moran_18002020_flood_prc value
            Q(i,:) =  Q_ARR_reach(:,idx_prc) + ( Q_Moran_18002020_flood_prc(i) - prc_values(idx_prc) ) ./ (prc_values(idx_prc+1) - prc_values(idx_prc))  .* ( Q_ARR_reach(:,idx_prc+1) - Q_ARR_reach(:,idx_prc)) ;

        end

        Q_scenarios{id_sim} = Q;

    end
    

end

clear prc_values range Q boundaries reach_sel idx_prcF idx_prc i n  prc prc_full

%% run the D-CASCADE model for a single scenario

id_sim = 4;
%1 - highQ
%2 - medQ
%3 - lowQ
%4 - mixQ

Q = Q_scenarios{id_sim};

% Find river network struct, containing 
Network = graph_preprocessing(ReachData);

%run D-CASCADE
data_plot = DCASCADEmodel_Bega( ReachData , Network , Q , Q_ARR_reach, n_Man , Width_change );

%% run the D-CASCADE model for all the discharge scenarios

for id_sim = 1:4  
    Q = Q_scenarios{id_sim};
    data_plot_multiple{id_sim,1} = DCASCADEmodel_Bega( ReachData , Network , Q , Q_ARR_reach,n_Man , Width_change );
end

%% define validation struct

% the validation parameters are contained in the SimVal struct, which is
% already defined and contains in the fields 'Field Data' and 'preES Data'
% the values of the validation parameters for the present-day and
% before European Settlement respectively

limt_2000 = 212; % timestep used for the extraction of the validation data
reach_type = [ReachData.reach_type]; %extract the reach type from ReachData

total_dep = [ReachData.deposit].*[ReachData.Length];

% SedDelRt
% the matrix SedDelRt_field contains all the available field measurements
% for the sediment delivery ratio  
for i=1:length(SedDelRt_field)
    for j = 1:size(data_plot_multiple,1)
        SimVal(j).([ 'SedDelRt_' num2str(SedDelRt_field(i,1))]) = data_plot_multiple{j,1} {9,2} (limt_2000,SedDelRt_field(i,1));
    end
end
 
% Type 4 - % incised
prc_eroded_sw = num2cell(1 - cell2mat(cellfun( @(x)sum(x{3,2}( limt_2000 , reach_type == 4  )) , data_plot_multiple(:,1),'UniformOutput',0)) ./ sum(total_dep(reach_type == 4)) )     ;

[SimVal(1:size(data_plot_multiple,1)).prc_eroded_sw] = prc_eroded_sw{:};

% Type 3 - Width
Wac_ll =  cellfun( @(x)mean(x{1,2}( limt_2000 , 94:99  )) , data_plot_multiple(:,1),'UniformOutput',0 );
[SimVal(1:size(data_plot_multiple,1)).Wac_ll] = Wac_ll{:} ;

% Type 3 - D50
D50_ll =  cellfun( @(x)mean(x{5,2}( limt_2000-6:limt_2000+6 , 94:99  ),'all') , data_plot_multiple(:,1),'UniformOutput',0 );
[SimVal(1:size(data_plot_multiple,1)).D50_ll] = D50_ll{:} ;

clear i limt_2000 prc_eroded_sw Wac_ll D50_ll total_dep change_class

%% compact the D-CASCADE outputs for the reaches in the lowland sections displayed in the figures

% reaches_compact_group contains the IDs of the reaches to be united into a
% single section

data_plot_multiple_compacted = data_plot_multiple;

for m=1:size(data_plot_multiple,1)
    for i=1:length(reach_compact_group)
    reaches_compact = reach_compact_group{i};

        for j=1:9

            if any(j==[2,3])
                data_plot_multiple_compacted{m,1}{j,2}(:,reaches_compact) = repmat(sum(data_plot_multiple_compacted{m,1}{j,2}(:, reaches_compact),2) , [1, length(reaches_compact)] ) ;
            else
                data_plot_multiple_compacted{m,1}{j,2}(:,reaches_compact) = repmat(mean(data_plot_multiple_compacted{m,1}{j,2}(:, reaches_compact),2) , [1, length(reaches_compact)] ) ;       
            end

        end
    data_plot_multiple_compacted{m,1}{10,2}{1}(:,reaches_compact) = repmat(sum(data_plot_multiple_compacted{m,1}{10,2}{2}(:, reaches_compact),2) , [1, length(reaches_compact)] ) ;
    data_plot_multiple_compacted{m,1}{10,2}{2}(:,reaches_compact) = repmat(sum(data_plot_multiple_compacted{m,1}{10,2}{2}(:, reaches_compact),2) , [1, length(reaches_compact)] ) ;
    end
end

clear m j reaches_compact