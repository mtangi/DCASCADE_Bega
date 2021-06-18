function Q_cell_prov = ExtractSedReaches(Q_cell,reaches, id)
% EXTRACTSEDREACHES is a support function used to extract only the cascades
% with provenance defined in vector reaches, from the deposit layer Q_cell.
% if id = 1 the function extract the cascade with provenance = reaches,
% otherwirse it extract the all the cascades not originated in reaches.

%%
if id == 1
    id_reaches =  any( Q_cell(:,1) == reaches',2 ) ;
else
    id_reaches =  sum( Q_cell(:,1) ~= reaches',2 ) == length(reaches) ;
end

Q_cell_prov = Q_cell(id_reaches,2:end);

end

