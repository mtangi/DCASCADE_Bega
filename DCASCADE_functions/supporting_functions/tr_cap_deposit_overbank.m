function [V_dep2act_new, V_dep ] = tr_cap_deposit_overbank(V_dep2act, V_dep, tr_cap_overbank)
%TR_CAP_DEPOSIT_OVERBANK deposit part of the deposited sediment volumes in
%the active layer according to the transport capacity of the reach
%considering the flow corrected for overbank. The result is a new V_dep2act
%that will be used in the next step to determine the mobilized volume Vmob.

%This function is only used if the volume contained in V_dep2act is higher
%that the tr.cap overbank of the reach, and is used to limite the sediemnt
%that can be lifted form the reach

% the overall structure of the function is similar to tr_cap_deposit

%% V_dep and V_act identification

global roundpar 

class_sup_dep = tr_cap_overbank < sum(V_dep2act(:,2:end),1); % classes for which the tr_cap_ovdrbank is less than the volumes in the active layer from the deposit

%if there are sed classes for which the tr cap is lower than the volume in
%V_dep2act...

if any(class_sup_dep)
    %...  sediment from V_dep2act will have to be deposited back in V_dep,
    % taking into consideration the sediment stratigraphy 
    %(lower layers get deposited first) 
    
    V_to_be_eroded = tr_cap_overbank(class_sup_dep); % total volume per class to be eroded
    V_dep2act_class = V_dep2act(:,[0 class_sup_dep]==1); % i take only the columns with the cascades of the classes class_sup_dep
    
    csum = cumsum(V_dep2act_class,1,'reverse'); %cumulative volume

    % find the indexes of the first cascade above the tr_cap threshold, for each class
    map=bsxfun(@gt,csum,V_to_be_eroded);
    
    map(1,sum(map,1) == 0) = 1;
    
    % find position of the layer to be splitted between deposit and erosion
    [~, firstoverthresh] = min(map,[],1);
    firstoverthresh = firstoverthresh - 1;
    firstoverthresh(firstoverthresh==0) = size((csum),1); 
    
    mapfirst = zeros(size(map));
    mapfirst( sub2ind(size(V_dep2act_class), firstoverthresh, 1:1:sum(class_sup_dep)) ) = 1;

    perc_dep = min(  (V_to_be_eroded - sum(V_dep2act_class.*~map,1)) ./ V_dep2act_class(sub2ind(size(V_dep2act_class), firstoverthresh, 1:1:sum(class_sup_dep))) ,1); %percentage to be lifted from the layer "on the threshold"

    map_perc = mapfirst.*perc_dep + ~map;
    
    % the matrix V_dep2act_new contains the mobilized cascades from the deposit layer, now corrected according to the tr_cap
    V_dep2act_new = zeros(size(V_dep2act));
    V_dep2act_new(:,1) = V_dep2act(:,1);
    V_dep2act_new(:,[0 class_sup_dep]==1) = map_perc.* V_dep2act_class;
    
    if ~isnan( roundpar );  V_dep2act_new(:,2:end)  = round( V_dep2act_new(:,2:end) , roundpar );   end

    % the matrix V_2dep contains the cascades that will be deposited into
    % the deposit layer.
    %  (the new volumes for the classes in class_sup_dep and all the
    %   volumes in the remaining classes)
    V_2dep = zeros(size(V_dep2act));
    V_2dep(:,[1 ~class_sup_dep]==1) = V_dep2act(:, [1 ~class_sup_dep]==1);   
    V_2dep(:,[0 class_sup_dep]==1) = (1 - map_perc) .* V_dep2act_class;
    
    if ~isnan( roundpar );  V_2dep(:,2:end)  = round( V_2dep(:,2:end) , roundpar );   end

    % remove empty rows (if the matrix is not already empty)
    if sum(V_dep2act(:,2:end),'all') ~= 0
        V_dep2act_new(sum(V_dep2act(:,2:end),2)==0,:) = [];
    end
    
else
    
    V_dep2act_new = V_dep2act;
    V_2dep = zeros(1, size(V_dep2act,2));
    V_2dep(1) = 1;

end



%% Put the volume exceeding the overbank transport capacity back in the deposit

%If the upper layer in the deposit and the lower layer in the volume to be
%deposited are from the same reach, i sum them
if (V_dep (end,1) == V_2dep(1,1))
        
    V_dep (end,2:end) = V_dep (end,2:end) + sum(V_2dep(1,2:end),1);
    V_dep = [V_dep; V_2dep(2:end,:)];   
        
else
    
    V_dep = [V_dep; V_2dep];
    
end
    
% remove empty rows
if sum(V_dep(:,2:end),'all')~=0
    V_dep(sum(V_dep(:,2:end),2)==0,:) = [];
end

end