function [data_output,extended_output] = DCASCADEmodel_Bega( ReachData , Network , Q , Q_ARR_reach, n_Man , Width_change )
%% D-CASCADE for the Bega river network case study

%input:
% ReachData      = n x 1 Struct defining the features of the n network reaches. The network features refers to pre-ES conditions
% Network        = 1x1 struct containing for each node info on upstream and downstream nodes
% Q              = n x t meatrix containing daily discharge [m3/s] for each reach in each timestep t
% Q_ARR_reach    = Bega-specific matrix, reporting discharge for each reach in each of the 7 flood classes considered [1- 2- 5- 10- 20- 50- 100-year return period]
% n_Man          = n x t matrix containing the Manning's roughness parameter for each reach in each timestep t
% Width_change   = 3xn matrix reporting the changes in the pre-ES channel width for eahc reach in phase 1, 2 and 3 of the simulation

%% Static variables extraction from ReachData 

%reach properties
Lngt = [ReachData.Length];
NH = Network.NH;

Slope = [ReachData.Slope]; 

reach_type = [ReachData.reach_type]'; %reach change type

global roundpar;
roundpar = 1;
  
%% Sediment classes definition
   
sed_range = [-5.5 , 3.5]; %range of sediment sizes considered in the model
class_size = 1.5; %amplitude of the sediment classes

global psi
psi =  sed_range(1):class_size:sed_range(2);   

%% outlet and sources reaches identification

outlet = Network.NH(end);

sources = find (cellfun(@isempty,Network.Upstream.Node)==1);

%% time initialization 

time_intervals = [50 50 10 70 37 1 ];

timescale = sum(time_intervals);

%% available deposit definiton
    
dep = [ReachData.deposit];
total_dep = dep.*[ReachData.Length];
dep_increase = zeros(timescale, length(ReachData));

if isnan(roundpar)
    %increase deposit for lowland reaches in phase 1, with a linear increase in deposit
    dep_increase(time_intervals(1)+1:time_intervals(1)+time_intervals(2),reach_type ~= 4) = repmat([total_dep(reach_type ~= 4)./time_intervals(2)] ,time_intervals(2) , 1 );
else
    dep_increase(time_intervals(1)+1:time_intervals(1)+time_intervals(2),reach_type ~= 4) = round(repmat([total_dep(reach_type ~= 4)./time_intervals(2)] ,time_intervals(2) , 1 ) , roundpar );
end
    
%% Wac increase definiton

Width = zeros(timescale, length(ReachData));
Width(1,:) = [ReachData.Wac];

Width_change_phase1 = Width_change (1,:) ; %change in width due to dbank erosion (type 2-3)

Width_change_phase2 = Width_change (2,:) ; %change in width between swamps and channelized valley fills (type 4)

Width_change_phase3 = Width_change (3,:) ; %change in width caused by incision in valley fills (type 4)

%% variables initialization 

%Variables to be calculated during the routing
clear Qbi_tr; Qbi_tr  = cell(timescale,1);     %Qbi_tr reports the cascades mobilized present in the reach AFTER transfer
clear Qbi_mob; Qbi_mob  = cell(timescale,1);   %Qbi_mob reports the cascades mobilized present in the reach BEFORE transfer

clear Qbi_dep; Qbi_dep  = cell(timescale,length(Network.II));

clear Fi_r_act; Fi_r_act = cell(timescale,1); %Fi_r_act contains the sed distribution of the active layer
clear D50_AL; D50_AL = zeros(timescale,length(NH)); %D50_AL contains the sed D50 of the active layer in each timestep
clear Qbi_input; Qbi_input  = cell(timescale,1); [Qbi_input{:}] = deal(zeros(length(Network.II), length(psi)));

clear Q_out ; Q_out = cell(timescale,1);  %Q_out reports the cascades leaving the network trought the outlet 

clear avg_vel; avg_vel = zeros(timescale,length(NH)); %average sed. velocoty in each timestep in each reach
    
%% variable initialization for 1st timestep
    
% sediment velocity parameters
phi = 0.4; 
min_velocity = 0.00001;

% initialize matrixes for first timestep 
Qbi_tr{1} = zeros([size(Network.II),length(psi)]); 
Qbi_mob{1} = zeros([size(Network.II),length(psi)]); 

Q_out{1} = zeros( length(Network.II), length(psi));
D50 = [ReachData.D50];

D50_AL(1,:) = D50;

% deposit lauer initialization
startlayer = 2; % define number of cascades into which to divie the initialized deposit layer

for n = NH 
    
    % initialize the deposiut only in the type 4 reaches, in the other,
    % material will be added as n_Man decreases in the first phase
    if reach_type(n) == 4

        Fi_r_act{1}(:,n) = GSDcurvefit( ReachData(n).D16 , D50(n), ReachData(n).D84 ); %Find the GSD for each sed. class given the D16, D50 and D84 of each reach
        Qbi_dep{1,n} = round(repmat([n, total_dep(n) / startlayer .* Fi_r_act{1}(:,n)' ], startlayer,1), roundpar);
    else  
        Fi_r_act{1}(:,n) = GSDcurvefit( ReachData(n).D16 , D50(n), ReachData(n).D84 );
        Qbi_dep{1,n} = [n, zeros(1,length(psi)) ];
    end
end  

%set limit for erosion in 1 timestep, given by the parameter mlim, that
%is the maximum depth in meter that can be eroded
mlim = ones(1,length(ReachData))  * 1 ; 
mlim(reach_type == 4) = 10;
V_lim_tot = round (mlim .* Width(1,:) .* Lngt , roundpar) ;

% tr_cap_sum contains the total transport capacity in each reach for each
% timestep
tr_cap_sum = zeros(timescale, length(ReachData));

%overbank_values indicates for each reach and for each t which of the return periods defined in prc_values (equal to [1 2 5 10 20 50 100] for the Bega) 
%to reduce the flow to if higher, ONLY for determinig the volume eroded. This indicates that if the flow is higher the water is
%actually overbank and thus not capable of eroding the river bed.

overbank_values = repmat([ReachData.overbank_prc], [timescale 1 ]);
overbank_values(1:sum(time_intervals(1:2)), or(reach_type == 2 , reach_type == 3 ) ) = 7;

%% Routing scheme

wb = waitbar(1/timescale, ['timestep 1/' num2str(timescale)]); %open waitbar
timerout = 1000000;

for t = 2:timescale

    %print waitbar
    time1   = clock;
    waitbar(t/timescale, wb, ['timestep ' num2str(t) '/' num2str(timescale) '  -  ' num2str(ceil(timerout * (timescale-t))) ' sec left' ]); % update waitbar

    %variables initialization
    Qbi_tr{t} = zeros( size(Qbi_tr{1}) );
    Qbi_mob{t} = zeros( size(Qbi_mob{1}) );

    Q_out{t} = zeros( size(Q_out{1}) );

    %calculate new water dept and velocity for all reaches w/ Manning
    h = (Q(t,:).*n_Man(t,:)./(Width(t-1).*sqrt( Slope ))).^(3/5);
    v = 1./n_Man(t,:).*h.^(2/3).*sqrt( Slope );

    %loop for all reaches
    for n = NH

         % add new sed deposit made available by deforestation 
        if dep_increase(t,n) ~= 0
            V_dep_old = [[ n dep_increase(t,n) .* Fi_r_act{1}(:,n)'] ; Qbi_dep{t-1,n} ];
            V_dep_old(:,2:end) = round(V_dep_old(:,2:end),roundpar); 

        else
            V_dep_old = Qbi_dep{t-1,n}; % extract the deposit layer of the reach from the relative cell in the previous timestep
        end

        % this component assures that the deposit from deforestation is
        % always on top of the deposit layer 
        if reach_type(n) ~= 4
            to_be_eroded = sum(V_dep_old(V_dep_old(:,1) == n,2:end),1);
            V_dep_old(V_dep_old(:,1) == n,:) = []; 
            V_dep_old = [ V_dep_old ; [n to_be_eroded]];
        end

        if any(n == sources) 

            Qbi_dep{t,n} = V_dep_old; 
            Fi_r_act{t}(:,n) = Fi_r_act{t-1}(:,n);
            D50_AL(t,n) = D_finder(Fi_r_act{t}(:,n));

            tr_cap = Ackers_White_tr_cap( Fi_r_act{t}(:,n) , Slope(n) , Q(t,n) , v(n) , h(n) )' .* 24.*60.*60;

            V_mob = [n 0.*Fi_r_act{t}(:,n)' ] ;

            Qbi_mob{t}(V_mob(:,1),n,:) = round( V_mob(:,2:end), roundpar);  % Qbi_mob contains the volume mobilized in the reach, that is about to be transfer downstream

            % i sum the volumes transported in reach n with all the other
            % volumes mobilized by all the other reaches in time t
            tr_cap_sum(t,n) = sum(tr_cap);

        else
            %%% Step 1a) extract the deposit layer from the storage matrix and load the incoming cascades

            Qbi_incoming = [(1:length(NH))' squeeze(Qbi_tr{t-1}(:, n,:)) ; n Qbi_input{t}(n,:) ]; 
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
            [V_inc2act , V_dep2act ,  V_dep , Fi_r_act{t}(:,n)] = layer_search (Qbi_incoming, V_dep_old, V_lim_tot(n));

            %if i have almost no sediment in my reach (less then 1/100 of the active layer limit, i put the sed distribution to the values in the previous timestep)
            if sum(Qbi_incoming(:,2:end),'all') + sum(V_dep_old(:,2:end),'all') < V_lim_tot(n)/100
                Fi_r_act{t}(:,n) = Fi_r_act{t-1}(:,n);
            end

            if sum(Fi_r_act{t}(:,n))==0; Fi_r_act{t}(:,n) = Fi_r_act{t-1}(:,n); end % in case the active layer is empty, i use the GSD of the previous timesteep

            D50_AL(t,n) = D_finder(Fi_r_act{t}(:,n));

            %calculate transport capacity using the Fi of the active layer, the resulting tr_cap is in m3/s and is converted in m3/day
            tr_cap = Ackers_White_tr_cap( Fi_r_act{t}(:,n) , Slope(n) , Q(t,n) , v(n) , h(n) )' .* 24.*60.*60;

            tr_cap_sum(t,n) = sum(tr_cap);

            % ADD-ON COMPONENT: Overbank flooding 
            %if the river goes overbank, the volume that can be eroded
            %is limited by the flow actually in the river, defined by the
            %values in overbank_values 
            
            if Q(t,n) > Q_ARR_reach(n,overbank_values(t,n))

                Q_overbank = Q_ARR_reach(n,overbank_values(t,n));
                h_overbank = (Q_overbank.*n_Man(t,n)./(Width(t-1,n).*sqrt( Slope(n) ))).^(3/5);
                v_overbank = 1./n_Man(t,n).*h_overbank.^(2/3).*sqrt( Slope(n));

                tr_cap_overbank = Ackers_White_tr_cap( Fi_r_act{t}(:,n) , Slope(n) ,Q_overbank , v_overbank , h_overbank )' .* 24.*60.*60;

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
            Qbi_dep{t,n} = V_dep;

            %remove empty rows
            Qbi_dep{t,n}(sum(Qbi_dep{t,n}(:,2:end),2)==0,:)=[];

            % Qbi_mob contains the volume mobilized in the reach, that is about to be transfer downstream
            Qbi_mob{t}(V_mob(:,1),n,:) = V_mob(:,2:end) ; 

            % if removing empty rows leaves only an Qbi_dep{t,n} empty
            % matrix, put an empty layer
            if isempty( Qbi_dep{t,n}); Qbi_dep{t,n} = [n zeros(1,length(psi))]; end

            % if we are in the initialization phase of the simulation, do not change deposition and elevation
            if t < time_intervals(1)
                Qbi_dep{t,n} = Qbi_dep{t-1,n};
            end

        %(end of the source if)
        end

    %%% Step 2) Update network width with channel expansion add-on
    
    %ADD_ON COMPONENT: Width change    
    % increase the width according to the material that was removed (if there was material added)
    if sum( dep_increase(1:t,n) ) > 0 && reach_type(n) ~= 4
        % The sum of the material not yet available and the material
        % still in the reach is the total material no yet eroded.
        % The Width increase proportionally to the percentage of
        % material eroded 
        Width_change = (1-(sum( dep_increase(t+1:end,n) )  +  sum( Qbi_dep{t,n}( Qbi_dep{t,n}(:,1) == n ,2:end),'all'))/total_dep(n) )  * Width_change_phase1(n);   
        Width(t, n) =  Width(1,n) + Width_change;     

    elseif  reach_type(n) == 4 && t > sum(time_intervals(1:2)) && t <= sum(time_intervals(1:3)) % channel width reduction from drainage in type 4 reaches in phase 2

        Width_change =  Width_change_phase2(n) / time_intervals(3) *  min((t -  sum(time_intervals(1:2))) , time_intervals(3) ) ;
        Width(t, n) =  Width(1,n) + Width_change;   

    elseif  reach_type(n) == 4 && t > sum(time_intervals(1:3)) && total_dep(n)>0 % channel width expansion via incision in type 4 reaches in phase 3

        % we expand the channel proportionally to the erosion of the material available
        Width_change = (total_dep(n) -  sum( Qbi_dep{t,n}(Qbi_dep{t,n}(:,1) == n,2:end) ,'all'))/ total_dep(n) * Width_change_phase3(n) ; 
        Width(t, n) =  Width(sum(time_intervals(1:3)),n) + Width_change;   

    else
        Width(t, n) = Width(t-1, n);
    end

    
    %(end of the reach loop)
    end

    %%% Step 3) Move the mobilized volumes to the destination reaches according to the sediment velocity

    % compute sediment velocity in the timestep  
    [ v_sed ] = velocity_AW( Fi_r_act{t} , Slope , Width(t-1,:) , Q(t,:) , v , h , min_velocity ) ;

    avg_vel(t,:) = mean(v_sed);

    %loop for all reaches again, now that i have the Fi_r and thus can compute transfer rates of the sediment volumes for all reaches
    for n = NH
        V_mob = zeros(length(NH),length(psi)+1);
        V_mob(:,1) = 1:length(NH);

        V_mob(:,2:length(psi)+1) = squeeze(Qbi_mob{t}(:,n,:));
        V_mob = matrix_compact(V_mob);
        Fi_mob = sum(V_mob(:,2:end),1)'./sum(V_mob(:,2:end),'all');

        if isnan(Fi_mob); Fi_mob = Fi_r_act{t}(:,n);           end

        [ v_sed ] = velocity_AW( repmat( Fi_mob,[1 length(NH)] ) , Slope , Width(t-1,:) , Q(t,:) , v , h , min_velocity ) ;

        [Qbi_tr_t, Q_out_t  ] = sed_transfer( matrix_compact(V_mob) , n , v_sed.*(60*60*24) , Lngt, Network );

        % i sum the volumes transported in reach n with all the other
        % volumes mobilized by all the other reaches in time t
        Qbi_tr{t} = Qbi_tr{t} + Qbi_tr_t;
        Q_out{t} =  Q_out{t} + Q_out_t;
    end

    %change the active layer volume according to the new width 
    V_lim_tot = round (mlim .* Width(t,:) .* Lngt , roundpar) ;

    %measure time of routing
    time2   = clock;

    if mod(round(timescale/10) , t) == 0  %save time only at certain timesteps
         timerout = etime(time2, time1);
    end


%end of the time loop
end

close(wb); clear wb;

%% output identification

% aggregated matrixes
QB_mob = cellfun(@(x)sum(x,3),Qbi_mob,'UniformOutput',0);
QB_mob_sum = cell2mat(cellfun(@(x)sum(x,[1,3]),Qbi_mob,'UniformOutput',0));
V_class_dep = cellfun(@(x)sum(x(:,2:end),1),Qbi_dep,'UniformOutput',0);

%total material in a reach in each timestep
tot_sed = cell2mat(cellfun(@(x)sum(x(:,2:end),'all'),Qbi_dep,'UniformOutput',0)) + cell2mat(cellfun(@(x)sum(x,1), QB_mob, 'UniformOutput', false)) ;

%Material in a reach in each timestep, classified according to their provenance 
%tot_sed_reaches{1} contains the material in a reach originated from type 1, 2 and 3 reaches
%tot_sed_reaches{1} contains the material in a reach originated from type 4 reaches

reaches = find(reach_type ~= 4);
tot_sed_reaches{1} = cell2mat( cellfun(@(x)sum(ExtractSedReaches(x, reaches, 1),'all') ,Qbi_dep,'UniformOutput',0) ) + cell2mat(cellfun(@(x)sum(x(reach_type ~= 4,:),1), QB_mob, 'UniformOutput', false)) + total_dep.*(reach_type~=4)' - cumsum(dep_increase) ;
tot_sed_reaches{2} = cell2mat( cellfun(@(x)sum(ExtractSedReaches(x, reaches, 0),'all') ,Qbi_dep,'UniformOutput',0) ) + cell2mat(cellfun(@(x)sum(x(reach_type == 4,:),1), QB_mob, 'UniformOutput', false)) ;

%total material in a reach in each timestep including initial banks deposit
tot_sed_deposit = tot_sed + total_dep.*(reach_type~=4)' - cumsum(dep_increase);

% D50 of the material in a reach for each timestep 
Qbi_mob_class = cellfun(@(x)squeeze(sum(x,1)),Qbi_mob,'UniformOutput',0);
tot_sed_class = Qbi_mob_class;

D50_tot = zeros(size(Qbi_mob_class,1), size(V_class_dep,2));

for t = 1:length(Qbi_mob_class)
tot_sed_class{t} = Qbi_mob_class{t} + reshape(cell2mat(V_class_dep(t,:)),[length(psi) , size(V_class_dep,2) ])';
a_fi = tot_sed_class{t} ./ sum(tot_sed_class{t},2);
a_fi(isnan(a_fi)) = 0;
for i=1: size(V_class_dep,2)
    D50_tot(t,i) = D_finder(a_fi(i,:)');
end
end

% percentage of sediment leaving the subcatchment
% perc_sed_subcatch contains the percentage of sediment that left the
% subcatchment of each reach relative to the sediment volume at the
% beginning of the simulation

tot_sed_subcatch = zeros(size(Qbi_dep));

for t = 1:size(tot_sed_subcatch,1)
    for n = 1:  size(tot_sed_subcatch,2)   
        sub_nodes = find((Network.Upstream.Distance{n} ~= Inf) == 1);    % find nodes in the subcatchment     
        sub_nodes(sub_nodes == n) = [];
        tot_sed_subcatch(t,n) = sum(tot_sed_deposit(t,sub_nodes));
    end
end

perc_sed_subcatch =  (1 - tot_sed_subcatch ./ tot_sed_subcatch(2,:)) .*100;
perc_sed_subcatch(1,:) = perc_sed_subcatch(2,:);
perc_sed_subcatch(isnan(perc_sed_subcatch)) = 0;
perc_sed_subcatch(perc_sed_subcatch<0) = 0;

%% output struct definition
% data_plot contais the most important D_CASCADE outputs
data_output = cell(1,2);

data_output{1,2} = Width; 
data_output{2,2} = QB_mob_sum;
data_output{3,2} = tot_sed_deposit;
data_output{4,2} = n_Man;
data_output{5,2} = D50_tot.*1000; 
data_output{6,2} = tr_cap_sum;
data_output{7,2} = avg_vel.*60.*60.*24./1000;
data_output{8,2} = Q;
data_output{9,2} = perc_sed_subcatch;
data_output{10,2} = tot_sed_reaches;

data_output{1,1} = 'Wac';
data_output{2,1} = 'QB mob';
data_output{3,1} = 'Total sed in the reach';
data_output{4,1} = 'n_Man';
data_output{5,1} = 'D50 tot';
data_output{6,1} = 'trcap sum';
data_output{7,1} = 'avg km/d vel';
data_output{8,1} = 'Q';
data_output{9,1} = 'SedDelRt';
data_output{10,1} = 'tot sed - provenance';

% all other outputs are included in the extended_output cell variable
extended_output = cell(1,2);

extended_output{1,2} = Qbi_tr; 
extended_output{2,2} = Qbi_mob;
extended_output{3,2} = Q_out;
extended_output{4,2} = Qbi_dep;
extended_output{5,2} = Fi_r_act; 
extended_output{6,2} = D50_AL;
extended_output{7,2} = avg_vel;
extended_output{8,2} = dep_increase;
extended_output{9,2} = tot_sed_subcatch;

extended_output{1,1} = 'Qbi_tr'; 
extended_output{2,1} = 'Qbi_mob';
extended_output{3,1} = 'Q_out';
extended_output{4,1} = 'Qbi_dep';
extended_output{5,1} = 'Fi_r_ac'; 
extended_output{6,1} = 'D50_AL';
extended_output{7,1} = 'avg_vel';
extended_output{8,1} = 'dep_increase';
extended_output{9,1} = 'tot_sed_subcatch';

end