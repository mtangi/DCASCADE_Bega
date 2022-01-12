function [V_mob, V_dep ] = tr_cap_deposit(V_inc2act, V_dep2act, V_dep, tr_cap)
%TR_CAP_DEPOSIT deposit part of the incoming and deposited sediment volumes
%according to the transport capapcity

global roundpar

%% V_dep and V_act identification

class_sup_dep = tr_cap > sum(V_inc2act(:,2:end),1); % classes for which the tr_cap is more than the incoming volumes in the active layer

%if there are sed classes for which the tr cap is more than the volume in
%V_inc2act...

if any(class_sup_dep)
    %...  sediment from V_dep2act will have to be mobilized, taking into consideration
    % the sediment stratigraphy (upper layers get mobilized first)
    
    tr_cap_remaining = tr_cap(class_sup_dep) - sum(V_inc2act(:,[0 class_sup_dep]==1),1); %remaining active layer volume per class after considering V_inc2act
    V_dep2act_class = V_dep2act(:,[0 class_sup_dep]==1); % i take only the columns with the cascades of the classes class_sup_dep
    
    csum = cumsum(V_dep2act_class,1,'reverse'); %

    % find the indexes of the first cascade above the tr_cap threshold, for each class
    map=bsxfun(@gt,csum,tr_cap_remaining);
    
    map(1,sum(map,1) == 0) = 1;
    
    % find position of the layer to be splitted between deposit and erosion
    [~, firstoverthresh] = min(map,[],1);
    firstoverthresh = firstoverthresh - 1;
    firstoverthresh(firstoverthresh==0) = size((csum),1); 
    
    mapfirst = zeros(size(map));
    mapfirst( sub2ind(size(V_dep2act_class), firstoverthresh, 1:1:sum(class_sup_dep)) ) = 1;

    perc_dep = min(  (tr_cap_remaining - sum(V_dep2act_class.*~map,1)) ./ V_dep2act_class(sub2ind(size(V_dep2act_class), firstoverthresh, 1:1:sum(class_sup_dep))) ,1); %percentage to be lifted from the layer "on the threshold"

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

else   
    
    V_dep2act_new = zeros(1, size(V_dep2act,2));
    V_dep2act_new(1) = 1;
    V_2dep = V_dep2act; % I re-deposit all the matrix V_dep2act into the deposit layer
    
end

% for the classes where V_inc2act is enough, I deposit the cascades
% proportionally

perc_inc = tr_cap(~class_sup_dep)./ sum(V_inc2act(:,[0 ~class_sup_dep]==1),1);
perc_inc(isnan(perc_inc)) = 0; %change NaN to 0 (naN appears when both tr_cap and sum(V_inc2act) are 0)
class_perc_inc = zeros(size(class_sup_dep));
class_perc_inc(class_sup_dep==0) = perc_inc;

V_mob = matrix_compact([V_dep2act_new ; V_inc2act.*[1 class_sup_dep ] + V_inc2act.*[0 class_perc_inc] ]); % i put in the active layer the volumes of V_inc2act for the class_sup_dep classes

if ~isnan( roundpar );  V_mob(:,2:end)  = round( V_mob(:,2:end) , roundpar );   end

class_residual = zeros(size(class_sup_dep));
class_residual(class_sup_dep==0) = 1 - perc_inc;

V_2dep = [V_2dep ; V_inc2act.*[1 class_residual]];
   
if ~isnan( roundpar );  V_2dep(:,2:end)  = round( V_2dep(:,2:end) , roundpar );   end

%% Put the volume exceeding the transport capacity back in the deposit

%If the upper layer in the deposit and the lower layer in the volume to be
%deposited are from the same reach, i sum them
if (V_dep (end,1) == V_2dep(1,1))
        
    V_dep (end,2:end) = V_dep (end,2:end) + sum(V_2dep(1,2:end),1);
    V_dep = [V_dep; V_2dep(2:end,:)];   
        
else
    
    V_dep = [V_dep; V_2dep];
    
end
    
% remove empty rows
if sum(V_dep2act(:,2:end))~=0
    V_dep(sum(V_dep(:,2:end),2)==0,:) = [];
end

%% OPTIONAL : if the first n rows of V_2dep and the last row of V_dep are from the same reach, i can sum them.

%V_dep = [V_dep; V_2dep];

%V_dep = layer_compact([V_dep; V_2dep]);

% if (V_dep (end,1) == V_2dep(1,1))
%     
%     index_last = find(diff(V_2dep(:,1) == V_dep (end,1))  == -1 ,1); %find the index of the last consecutive row with same index as the last row of V_dep
%     
%     if isempty(index_last); index_last = size(V_2dep,1); end % if all the rows come from the same reach of the last row of V_dep, the index is the last row
%     
%     V_dep (end,2:end) = V_dep (end,2:end) + sum(V_2dep(1:index_last,2:end),1);
%     V_dep = [V_dep; V_2dep(index_last+1:end,:)];
%        
% else
%     
%     V_dep = [V_dep; V_2dep];
%     
% end

end

