function [Qbi_sort] = sortdistance(Qbi, distancelist )
% SORTDISTACE sort the rows of the Qbi_incoming matrix by 
% increasing distance from the reach

%% code
[index,~] = find(  Qbi(:,1) == distancelist(distancelist ~= inf ) );

if ~isempty(index)
    
    Qbi_sort = Qbi(index,:);
    
else
    
    Qbi_sort = Qbi;
    
end

if round(sum(Qbi_sort(:,2:end),'all'),5) ~= round(sum(Qbi(:,2:end),'all'),5)
    
    warning(['the sortdistance function produced results that are erroneous'] )

end
    
    

end
