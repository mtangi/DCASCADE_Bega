function [data_plot_scenarios] = DCASCADEmodel_Bega_future( ReachData , Network , Q_future_scenarios , Q_ARR_reach , data_plot, extended_output, ID_scenario )
%% D-CASCADE for the Bega river network case study, for simulation of future discharge scenarios

%input:
% ReachData      = n x 1 Struct defining the features of the n network reaches. The network features refers to pre-ES contitions
% Network        = 1x1 struct containing for each node info on upstream and downstream nodes
% Q_future_scenarios = n x 1 cell struct containing the discharge input of each scenario
% Q_ARR_reach    = Bega-specific matrix, reporting discharge for each reach in each of the 7 flood classes considered [1- 2- 5- 10- 20- 50- 100-year return period]
% data_plot , extended_output   = cell structs containing the outputs of the D-CASCADE simulation used as baseline for the future trajectories
% ID_scenario    = defines which future scenario we are simulating, either the present-day vegetation scenario (ID_scenario = 1) or the pre-recovery one

%% Define simulation parameters

%reach properties
Lngt = [ReachData.Length];
NH = Network.NH;

Slope = [ReachData.Slope]; 

reach_type = [ReachData.reach_type]'; %reach change type
dep = [ReachData.deposit];
total_dep = dep.*[ReachData.Length];

% sediment velocity parameters
min_velocity = 0.00001;

global roundpar;
roundpar = 1;
  
%% define future scenarios length and number

timescale = size(data_plot{8,2},1); % timescale indicates the numer of timesteps in the historic simulation
future_time = size(Q_future_scenarios{1},1);
n_scenarios = length(Q_future_scenarios);

%% Sediment classes definition
   
sed_range = [-5.5 , 3.5]; %range of sediment sizes considered in the model
class_size = 1.5; %amplitude of the sediment classes

global psi
psi =  sed_range(1):class_size:sed_range(2);   

%% sources reaches identification

sources = find (cellfun(@isempty,Network.Upstream.Node)==1);
    
%% Routing scheme

wb = waitbar(1/n_scenarios, ['scenario 1/' num2str(n_scenarios)]); %open waitbar
timerout = 1000000;

%loop for all scenarios
for s=1:n_scenarios
    
    %print waitbar
    time1   = clock;
    waitbar(s/n_scenarios, wb, ['scenario ' num2str(s) '/' num2str(n_scenarios) '  -  ' num2str(ceil(timerout * (n_scenarios-s)./60)) ' min left' ]); % update waitbar

    %variabiel initialization
    Q_future = [data_plot{8,2};Q_future_scenarios{s}];
    overbank_values_future = repmat([ReachData.overbank_prc], [timescale + future_time  1 ]);
    
    Qbi_tr_future = [extended_output{1,2};cell(future_time,1)];
    Qbi_mob_future = [extended_output{2,2};cell(future_time,1)];
    Q_out_future = [extended_output{3,2};cell(future_time,1)];
    Qbi_dep_future = [extended_output{4,2}; cell(future_time,length(NH))];
    Fi_r_act_future = [extended_output{5,2};cell(future_time,1)];
    D50_AL_future = [extended_output{6,2}; zeros(future_time,length(NH))];
    tr_cap_sum_future = [data_plot{6,2}; zeros(future_time,length(NH))];
    Wac_future = [data_plot{1,2}; zeros(future_time,length(NH))];
    avg_vel_future = [extended_output{7,2}; zeros(future_time,length(NH))];
    dep_increase_future = [ extended_output{8,2}; zeros(future_time,length(NH))];
    
    %set limit for erosion in 1 timestep, given by the parameter mlim, that
    %is the maximum depth in meter that can be eroded
    mlim = ones(1,length(ReachData))  * 1 ; 
    mlim(reach_type == 4) = 10;
    V_lim_tot = round (mlim .* Wac_future(timescale,:) .* Lngt , roundpar) ;

    if ID_scenario == 1 % present-day vegetation scenario
        n_Man_future = [data_plot{4,2}; repmat(data_plot{4,2}(end,:), [future_time,1])];
    else % pre-recovery vegetation scenario
        n_Man_future = [data_plot{4,2}; repmat(data_plot{4,2}(160,:), [future_time,1])];
    end

    %loop for each future timesteps
    for t = timescale: timescale + future_time
        
        %variables initialization
        Qbi_tr_future{t} = zeros( size(Qbi_tr_future{1}) );
        Qbi_mob_future{t} = zeros( size(Qbi_mob_future{1}) );

        Q_out_future{t} = zeros( size(Q_out_future{1}) );
        
        %calculate new water dept for all reaches
        
        %Manning
        h = (Q_future(t,:).*n_Man_future(t,:)./(Wac_future(t-1).*sqrt( Slope ))).^(3/5);
        v = 1./n_Man_future(t,:).*h.^(2/3).*sqrt( Slope);
        
%       % hydraulic solver
%       hydraulicData = hydraulic_solver(Slope , Q(t,:) , Wac, 0);
%       h = hydraulicData(:, 1);
           
        %loop for all reaches
        for n = NH
            
             % add new sed deposit made available by deforestation 

             V_dep_old = Qbi_dep_future{t-1,n}; % extract the deposit layer of the reach from the relative cell in the previous timestep         
                
            if any(n == sources)
                
                Qbi_dep_future{t,n} = V_dep_old; 
                Fi_r_act_future{t}(:,n) = Fi_r_act_future{t-1}(:,n);
                D50_AL_future(t,n) = D_finder(Fi_r_act_future{t}(:,n));
                
                tr_cap = Ackers_White_tr_cap( Fi_r_act_future{t}(:,n) , Slope(n) , Q_future(t,n) , v(n) , h(n) )' .* 24.*60.*60;
                
                %V_mob = [n min([ tr_cap  ; V_lim_tot(n).*Fi_r_act{t}(:,n)' ]) ] ;
                V_mob = [n 0.*Fi_r_act_future{t}(:,n)' ] ;
                
                Qbi_mob_future{t}(V_mob(:,1),n,:) = round( V_mob(:,2:end), roundpar);  % Qbi_mob contains the volume mobilized in the reach, that is about to be transfer downstream
                
                % i sum the volumes transported in reach n with all the other
                % volumes mobilized by all the other reaches in time t
                tr_cap_sum_future(t,n) = sum(tr_cap);

            else
                %%% Step 1a) extract the deposit layer from the storage matrix and load the incoming cascades

                Qbi_incoming = [(1:length(NH))' squeeze(Qbi_tr_future{t-1}(:, n,:)) ]; 
                Qbi_incoming(sum(Qbi_incoming(:,2:end),2)==0,:) = [];
                
                if isempty(Qbi_incoming); Qbi_incoming = [n zeros(1,length(psi))]; end %put an empty cascade if no incoming volumes are present (for computation)

                a=sum(Qbi_incoming(:,2:end),'all');
                %sort incoming matrix accoring to distance, in this way
                %sediment coming from closer reaches will be deposited first 
                [Qbi_incoming] = sortdistance(Qbi_incoming, Network.Upstream.distancelist{n} );

                if round(a,roundpar) ~= round(sum(Qbi_incoming(:,2:end),'all'),roundpar)
                    warning([' sorting error in reach ' num2str(n) ' at time ' num2str(t) ' = ' num2str( a - sum(Qbi_incoming(:,2:end),'all'))])
                end

                %%% Step 1b) find cascades to be included into the active layer according to the limit V_lim_tot, and use the cumulative GSD to compute tr_cap

                %find the volume of sediment from the incoming load (V_inc2act) and deposit layer (V_dep2act) to be included in the active layer
                [V_inc2act , V_dep2act ,  V_dep , Fi_r_act_future{t}(:,n)] = layer_search (Qbi_incoming, V_dep_old, V_lim_tot(n));
                    
                %if i have almost no sediment in my reach (less then 1/100 of the active layer limit, i put the sed distribution to the values in the previous timestep)
                if sum(Qbi_incoming(:,2:end),'all') + sum(V_dep_old(:,2:end),'all') < V_lim_tot(n)/100
                    Fi_r_act_future{t}(:,n) = Fi_r_act_future{t-1}(:,n);
                end

                if sum(Fi_r_act_future{t}(:,n))==0; Fi_r_act_future{t}(:,n) = Fi_r_act_future{t-1}(:,n); end % in case the active layer is empty, i use the GSD of the previous timesteep

                D50_AL_future(t,n) = D_finder(Fi_r_act_future{t}(:,n));

                %calculate transport capacity using the Fi of the active layer, the resulting tr_cap is in m3/s and is converted in m3/day
                tr_cap = Ackers_White_tr_cap( Fi_r_act_future{t}(:,n) , Slope(n) , Q_future(t,n) , v(n) , h(n) )' .* 24.*60.*60;

                tr_cap_sum_future(t,n) = sum(tr_cap);

                % Overbank flooding add-on
                %if the river goes overbank, the volume that can be eroded
                %is limited by the flow actually in the river, defined by the
                %values in overbank_values 
                if Q_future(t,n) > Q_ARR_reach(n,overbank_values_future(t,n))
                    
                    Q_overbank = Q_ARR_reach(n,overbank_values_future(t,n));
                    h_overbank = (Q_overbank.*n_Man_future(t,n)./(Wac_future(t-1,n).*sqrt( Slope(n) ))).^(3/5);
                    v_overbank = 1./n_Man_future(t,n).*h_overbank.^(2/3).*sqrt( Slope(n));

                    tr_cap_overbank = Ackers_White_tr_cap( Fi_r_act_future{t}(:,n) , Slope(n) ,Q_overbank , v_overbank , h_overbank )' .* 24.*60.*60;
                    
                    [V_dep2act, V_dep ] = tr_cap_deposit_overbank(V_dep2act, V_dep, tr_cap_overbank);
                end
                
                %%% Step 1c) Deposit the cascades in the active layer until the volume mobilized for each class is equal to the tr_cap

                if sum(tr_cap) < ( sum(V_dep2act(:,2:end),'all') + sum(V_inc2act(:,2:end),'all') ) %if the total transport capacity is lower than the active layer volume...
                      %... deposit part of the active layer cascades, 
                      %    proportionally to their volume and the volume of the active layer

                      [V_mob, V_dep ] = tr_cap_deposit( V_inc2act, V_dep2act, V_dep, tr_cap);      
                                        
                else
                    % if not, the mobilized layer is equal to the active
                    % layer
                    V_mob = matrix_compact([V_dep2act ; V_inc2act]);

                end 

                %%% Step 1d) after this passage, V_mob contains only the volumes actually mobilized    
                Qbi_dep_future{t,n} = V_dep;
                
                %remove empty rows
                Qbi_dep_future{t,n}(sum(Qbi_dep_future{t,n}(:,2:end),2)==0,:)=[];
                
                % Qbi_mob contains the volume mobilized in the reach, that is about to be transfer downstream
                Qbi_mob_future{t}(V_mob(:,1),n,:) = V_mob(:,2:end) ; 

                % if removing empty rows leaves only an Qbi_dep{t,n} empty
                % matrix, put an empty layer
                if isempty( Qbi_dep_future{t,n}); Qbi_dep_future{t,n} = [n zeros(1,length(psi))]; end
              
            %(end of the source if)
            end

        %%% Step 2) is skipped as we assume no further channel expansion
        %%% will occurr in the future
        
        Wac_future(t, n) =  Wac_future(t-1,n) ;

        %(end of the reach loop)
        end
        
        %%% Step 3) Move the mobilized volumes to the destination reaches according to the sediment velocity
        
        % compute sediment velocity in the timestep  
        [ v_sed ] = velocity_AW( Fi_r_act_future{t} , Slope , Wac_future(t-1,:) , Q_future(t,:) , v , h , min_velocity ) ;
        
        avg_vel_future(t,:) = mean(v_sed);
        
        %loop for all reaches again, now that i have the Fi_r and thus can compute transfer rates for all reaches
        for n = NH
            V_mob = zeros(length(NH),length(psi)+1);
            V_mob(:,1) = 1:length(NH);
            
            V_mob(:,2:length(psi)+1) = squeeze(Qbi_mob_future{t}(:,n,:));
            V_mob = matrix_compact(V_mob);
            Fi_mob = sum(V_mob(:,2:end),1)'./sum(V_mob(:,2:end),'all');
            
            if isnan(Fi_mob); Fi_mob = Fi_r_act_future{t}(:,n);           end
            
            [ v_sed ] = velocity_AW( repmat( Fi_mob,[1 length(NH)] ) , Slope , Wac_future(t-1,:) , Q_future(t,:) , v , h , min_velocity ) ;
            
            [Qbi_tr_t, Q_out_t  ] = sed_transfer( matrix_compact(V_mob) , n , v_sed.*(60*60*24) , Lngt, Network );

            % i sum the volumes transported in reach n with all the other
            % volumes mobilized by all the other reaches in time t
            Qbi_tr_future{t} = Qbi_tr_future{t} + Qbi_tr_t;
            Q_out_future{t} =  Q_out_future{t} + Q_out_t;

        end
                
        %change the erosion limit accroding to the new width 
        V_lim_tot = round (mlim .* Wac_future(t,:) .* Lngt , roundpar) ;       
            
    %end of the time loop
    end
    
    % output processing for the single future scenario
    QB_sum_future = cell2mat(cellfun(@(x)sum(x,[1,3]), Qbi_tr_future(timescale+1:end)   ,'UniformOutput',0));

    QB_mob_future = cellfun(@(x)sum(x,3),Qbi_mob_future,'UniformOutput',0);
    QB_mob_sum_future = cell2mat(cellfun(@(x)sum(x,[1,3]),Qbi_mob_future(timescale+1:end),'UniformOutput',0));

    V_class_dep_future = cellfun(@(x)sum(x(:,2:end),1),Qbi_dep_future(1:end,:),'UniformOutput',0);

    %Material in a reach in each timestep, classified according to their provenance 
    %tot_sed_reaches{1} contains the material in a reach originated from type 1, 2 and 3 reaches
    %tot_sed_reaches{1} contains the material in a reach originated from type 4 reaches

    tot_sed_future = cell2mat(cellfun(@(x)sum(x(:,2:end),'all'),Qbi_dep_future(timescale+1:end,:),'UniformOutput',0)) + cell2mat(cellfun(@(x)sum(x,1), QB_mob_future(timescale+1:end), 'UniformOutput', false)) ;
   
    reaches = find(reach_type ~= 4);
    tot_sed_reaches_future{1} = cell2mat( cellfun(@(x)sum(ExtractSedReaches(x, reaches, 1),'all') ,Qbi_dep_future(timescale+1:end,:),'UniformOutput',0) ) + cell2mat(cellfun(@(x)sum(x(reach_type ~= 4,:),1), QB_mob_future(timescale+1:end), 'UniformOutput', false)) ;
    tot_sed_reaches_future{2} = cell2mat( cellfun(@(x)sum(ExtractSedReaches(x, reaches, 0),'all') ,Qbi_dep_future(timescale+1:end,:),'UniformOutput',0) ) + cell2mat(cellfun(@(x)sum(x(reach_type == 4,:),1), QB_mob_future(timescale+1:end), 'UniformOutput', false)) ;

    %total material in a reach in each timestep including initial banks deposit
    csum = cumsum( dep_increase_future );
    tot_sed_deposit_future = tot_sed_future + total_dep.*(reach_type~=4)' - csum(timescale+1:end,:);

    % D50 of the material in a reach for each timestep 
    Qbi_mob_class = cellfun(@(x)squeeze(sum(x,1)),Qbi_mob_future,'UniformOutput',0);
    tot_sed_class = Qbi_mob_class;

    D50_tot = zeros(size(Qbi_mob_class,1), size(V_class_dep_future,2));

    for t = timescale:timescale + future_time
        
        tot_sed_class{t} = Qbi_mob_class{t} + reshape(cell2mat(V_class_dep_future(t,:)),[length(psi) , size(V_class_dep_future,2) ])';
        a_fi = tot_sed_class{t} ./ sum(tot_sed_class{t},2);
        a_fi(isnan(a_fi)) = 0;
        for i=1: size(V_class_dep_future,2)
            D50_tot(t,i) = D_finder(a_fi(i,:)');
        end
    end

    % percentage of sediment leaving the subcatchment
    % perc_sed_subcatch conains the percentage of sediment that left the
    % subcatchment of each reach relative to the sediment volume at the
    % beginning of the simulation

    tot_sed_subcatch_future = zeros(size(QB_sum_future));

    for t = 1:size(tot_sed_subcatch_future,1)
        for n = 1:  size(tot_sed_subcatch_future,2)   
            sub_nodes = find((Network.Upstream.Distance{n} ~= Inf) == 1);    % find nodes in the subcatchment 
            sub_nodes(sub_nodes == n) = [];
            tot_sed_subcatch_future(t,n) = sum(tot_sed_deposit_future(t,sub_nodes));
        end
    end

    perc_sed_subcatch =  (1 - tot_sed_subcatch_future ./ extended_output{9,2}(2,:)) .*100;
    perc_sed_subcatch(1,:) = perc_sed_subcatch(2,:);
    perc_sed_subcatch(isnan(perc_sed_subcatch)) = 0;
    perc_sed_subcatch(perc_sed_subcatch<0) = 0;

    % save the most important outputs of the future simulation as data_plot 
    data_plot_future = cell(1,2);

    data_plot_future{1,2} = Wac_future(timescale+1:end,:); 
    data_plot_future{2,2} = QB_mob_sum_future;
    data_plot_future{3,2} = tot_sed_deposit_future;
    data_plot_future{4,2} = n_Man_future(timescale+1:end,:);
    data_plot_future{5,2} = D50_tot(timescale+1:end,:).*1000; 
    data_plot_future{6,2} = tr_cap_sum_future(timescale+1:end,:);
    data_plot_future{7,2} = avg_vel_future(timescale+1:end,:).*60.*60.*24./1000;
    data_plot_future{8,2} = Q_future(timescale+1:end,:);
    data_plot_future{9,2} = perc_sed_subcatch;
    data_plot_future{10,2} = tot_sed_reaches_future;

    data_plot_future{1,1} = 'Wac';
    data_plot_future{2,1} = 'QB mob';
    data_plot_future{3,1} = 'Total sed in the reach';
    data_plot_future{4,1} = 'n_Man';
    data_plot_future{5,1} = 'D50 tot';
    data_plot_future{6,1} = 'trcap sum';
    data_plot_future{7,1} = 'avg km/d vel';
    data_plot_future{8,1} = 'Q';
    data_plot_future{9,1} = 'perc_sed_subcatch';
    data_plot_future{10,1} = 'tot sed - provenance';

    % store the data_plot cell for the scenario into data_plot_scenarios
    data_plot_scenarios{s,1} = data_plot_future;
    data_plot_scenarios{s,2} = ['avgQ: ' num2str(round(mean(Q_future(timescale+1:end,58)),1))];

    %measure time of routing for the waitbar
    time2   = clock;
    timerout = etime(time2, time1);

end

close(wb); clear wb;

end