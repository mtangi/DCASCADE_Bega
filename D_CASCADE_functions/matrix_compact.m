function [V_layer_cmpct] = matrix_compact(V_layer)
%MATRIX_COMPACT takes a stratigraphy matrix V_layer and compact it by
%summing all the layers with the same source reach id

% method 1 - accumarray using arrayfun and cell2mat

[ID,~,ord] = unique(V_layer(:,1));
V_layer_cmpct = [ID, cell2mat(arrayfun(@(x) accumarray(ord,V_layer(:,x)),2:size(V_layer,2),'un',0))];

if size(V_layer_cmpct,1)>1
    V_layer_cmpct(sum(V_layer_cmpct(:,2:end),2)==0,:) = [];
end

if isempty(V_layer_cmpct)
    V_layer_cmpct = [ID(1) zeros(1, size(V_layer(:,2:end),2))];
end
% 
%% method 2 - accumarray for each sed class
% [ID,~,ord] = unique(V_layer(:,1));
% 
% V_layer_cmpct = zeros(length(ID), size(Qbi_tr_t,3)+1);
% V_layer_cmpct(:,1) = ID;
% for i=1:size(Qbi_tr_t,3)
%     V_layer_cmpct(:,i+1) = accumarray(ord,V_layer(:,i+1));
% end
% 
% 
%% method 3 - for loop for each unique reach ID in V_act
%  
% ID = unique(V_layer(:,1));
% V_layer_cmpct = zeros(length(ID), size(Qbi_tr_t,3)+1);
% 
% if ~isempty(ID)
%     for i = ID
%        V_layer_cmpct(i,:) = [ID(i) sum(V_layer(V_layer(:,1) == i,2:end), 1)];
%     end
% end

end

