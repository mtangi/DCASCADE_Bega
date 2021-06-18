function [V_inc2act, V_dep2act , V_dep, Fi_r_reach] = layer_search(Qbi_incoming, V_dep_old, V_lim_tot_n)
%LAYER_SEARCH put part of the incoming and deposited sediment volumes into the
%active layer according to the active layer volume and the incoming and
%deposited volumes

global roundpar 

%% V_dep and V_act identification

 % if, considering the incoming volume, I am still under the threshold of the active layer volume...    
if V_lim_tot_n - sum(Qbi_incoming(:,2:end),'all') >0  

    % ... I put sediment from the deposit layer into the active layer
    V_lim_dep = V_lim_tot_n - sum(Qbi_incoming(:,2:end),'all'); %remaining active layer volume after considering incoming sediment cascades
    csum = cumsum(sum(V_dep_old(:,2:end),2),1,'reverse');

    V_inc2act = Qbi_incoming;  % all the incoming volume will end up in the active layer

    %find active layer 
    if isempty(find(csum>V_lim_dep,1,'last'))
         %if the cascades in the deposit have combined
         %volume that is less then the active layer volume
         %(i've reached the bottom)
         
         V_dep2act = V_dep_old;  % i put all the deposit into the active layer
         V_dep = [V_dep_old(1,1) zeros(1, size(Qbi_incoming,2)-1)]; % i put the deposit volumes to 0

    else
        % if i have multiple deposit layers, put the upper layers into the active layer until i reach the threshold. 
        % The layer on the threshold (defined by position index) gets divided according to perc_layer
         index = find(csum>=V_lim_dep,1,'last');

         perc_layer = (V_lim_dep - sum(V_dep_old(csum<V_lim_dep,2:end),'all'))/sum(V_dep_old(index,2:end)); %percentage to be lifted from the layer on the threshold
         
         perc_layer = max(0,perc_layer); % remove small negative values that can arise from the difference being very close to 0
         
         if ~isnan( roundpar )
             V_dep2act = [[V_dep_old(index,1) round(V_dep_old(index,2:end).*perc_layer,roundpar) ]; V_dep_old(csum<V_lim_dep,:)];
             V_dep = [V_dep_old(1:index-1,:); [V_dep_old(index,1) round(V_dep_old(index,2:end).* (1-perc_layer), roundpar)]];
         else
             V_dep2act = [[V_dep_old(index,1) V_dep_old(index,2:end).*perc_layer]; V_dep_old(csum<V_lim_dep,:)];
             V_dep = [V_dep_old(1:index-1,:); [V_dep_old(index,1) V_dep_old(index,2:end).* (1-perc_layer)]];
         end
    end

    % in the end, except for the case find(csum<V_lim_dep,1) == 1, sum([V_dep2act V_inc2act],'all') == V_lim_tot_n 

else %if the incoming sediment volume is enough to completely fill the active layer...

    %... deposit part of the incoming cascades 
    %    proportionally to their volume and the volume of the active layer, 
    %    and put the rest into the active layer

    perc_dep = V_lim_tot_n / sum(Qbi_incoming(:,2:end),'all');  %percentage of the volume to be put in the active layer for all the cascades
    
     if ~isnan( roundpar )
         Qbi_incoming_dep = round(Qbi_incoming(:,2:end).*(1-perc_dep), roundpar);
     else
         Qbi_incoming_dep = Qbi_incoming(:,2:end).*(1-perc_dep); %this contains the fraction of the incoming volume to be deposited
     end
     
    V_inc2act = [Qbi_incoming(:,1) (Qbi_incoming(:,2:end) - Qbi_incoming_dep) ];
    V_dep2act = [V_dep_old(1,1) zeros(1, size(Qbi_incoming,2)-1)];
    
    if any(sum(Qbi_incoming(:,2:end).*(1-perc_dep))) ~= 0 %if, given the round, the deposited volume of the incoming cascades is 0...
        V_dep = [V_dep_old ; [Qbi_incoming(:,1) Qbi_incoming_dep]];
    else
        V_dep = V_dep_old; %... i leave the deposit as it was.
    end
end

%V_inc2act = matrix_compact(V_inc2act);

% remove empty rows (if the matrix is not already empty)
if sum(V_dep2act(:,2:end))~=0
    V_dep2act(sum(V_dep2act(:,2:end),2)==0,:) = [];
end


%% find active layer GSD


Fi_r_reach =  ( sum(V_dep2act(:,2:end),1) + sum(V_inc2act(:,2:end),1) )./  ( sum(V_dep2act(:,2:end),'all') + sum(V_inc2act(:,2:end),'all') ); %i find the GSD of the active layer, for the transport capacity calculation
Fi_r_reach(isnan(Fi_r_reach)) = 0;  %if V_act is empty, i put Fi_r equal to 0 for all classes
    

   % elseif isempty(find(csum<V_lim_dep,1))
        %if the first cascade in the deposit has volume
        %that is more then the remaining active layer volume
%         index = size(V_dep,1) ;
%         perc_layer = V_lim_dep/sum(V_dep(index,2:end)); %percentage to be lifted from the single cascade
% 
%         V_act = [ [V_dep(index,1) V_dep(index,2:end).*perc_layer] ; Qbi_incoming]; % i put part of the first cascade into the active layer, along with the incoming cascades
% 
%         V_dep(index,2:end) = V_dep(index,2:end).*(1-perc_layer);
% 
%     else



end

