function [Qbi_tr_t, Q_out_t , setplace, perc_out] = sed_transfer(V_mob , n , v_sed_day , Lngt, Network )
% SED_TRANSFER takes the matrix of the mobilized layers(V_mob) and the
% vector of the sed velocity for each class(v_sed_id) in a reach (n) and
% returns the 3D matrices containing the sed volumes in V_mob in their new
% position at the end of the timestep. 

% back = point in the sed.volume placed on the FromN
% front = point in the sed.volume placed on the ToN

%% initailize parameters
global roundpar 
outlet = Network.NH(end);

%% find start and end reach of the sed volume after the timestep

    %reach_s is the id of the reach where the back of the sed. volume stops after the timestep 
    %p_start is the position from the from_node of the id reach where the back of the sed. volume stops after the timestep 

    start_pos = max( 1 - min( v_sed_day(:,n) , Lngt(n)) / Lngt(n) , 0) + 0.00000001; % the buffer is to avoid errors if the parcel falls exactly on the n reach ToNode
    [p_start , reach_s] = track_sed_position( n , v_sed_day , Lngt , Network , start_pos);
    
    %reach_e is the id of the reach where the front of the sed. volume stops after the timestep 
    %p_end is the position from the from_node of the id reach where the front of the sed. volume stops after the timestep 

    if n == outlet  
       reach_e = repmat( n , length(reach_s) ,1);
       p_end = Lngt(n) + v_sed_day(:,n);    
    else   
        %to find p_end, i track the position of a sediment parcel starting
        %from the To_Node of the reach n (i.e. the From_node of the downstream reach).
        
        [p_end , reach_e] = track_sed_position( cell2mat(Network.Downstream.Node(n)) , v_sed_day , Lngt , Network ); 
        p_end = p_end + Lngt(n);
    end

    % the length of the volume after the transition Lngt_vol is needed to
    % compute the fractions
    Lngt_vol = max( p_end - p_start , 0.00001);

    %downdist contains the d
    downdist = Network.Downstream.Distance{n,1} ;
    
%% find fractions of the volume to be attributed in each reach

    %perc_s contains the fraction of the volume in reach_s at the end of the timestep
    %perc_e contains the fraction of the volume in reach_e at the end of the timestep

    perc_s = (Lngt(reach_s)' + downdist(reach_s)' - p_start) ./ Lngt_vol;
    perc_e = (p_end - downdist(reach_e)') ./ Lngt_vol;

    %if the entire volume mobilized is contained in the same reach (and thus
    %perc1 and perc2 >1), set one perc values at 1 and the other at 0
    flag_fill = or( and(perc_s >= 1, perc_e >= 1) , Lngt_vol == 0.00001) ;
    perc_s = perc_s.*(~flag_fill) + flag_fill;
    perc_e = perc_e.*(~flag_fill);

    %if some of the sediment volume mobilized leave from the outlet, i find
    %the fraction of leaving volume and save it in perc_out
    perc_out = zeros(size(perc_e));

    if any(p_end - Lngt(reach_e)' - downdist(reach_e)' > 0)

        find_rows = find(p_end - Lngt(reach_e)' - downdist(reach_e)' > 0);
        %perc_out(find_rows) = min( (p_end(find_rows) - Lngt(reach_e)' - downdist(reach_e)') / Lngt(id),1);

        perc_out(find_rows) = (p_end(find_rows) - Lngt(reach_e(find_rows))' - downdist(reach_e(find_rows))') ./ Lngt_vol(find_rows) ;

        % if all the volume leave (perc_out >=1) set all the other par to 0
        perc_e(perc_out>=1) = 0;
        perc_s(perc_out>=1) = 0;
        perc_out(perc_out>=1) = 1;

        % if the volume is spread among outside the network, the outlet reach and others
        % upstream
        index_spread = and(p_end - Lngt(reach_e)' - downdist(reach_e)' > 0 ,reach_s ~= reach_e); %rows where the spread happens
        perc_e( index_spread ) = max(perc_e(index_spread) - perc_out(index_spread),0);

        % if the volume is spread between outside the network and the outlet reach
        index_not_spread = and(p_end - Lngt(reach_e)' - downdist(reach_e)' > 0 ,reach_s == reach_e);
        perc_e(index_not_spread) = 0;

    end

    % matrix setplace contains for each reach and for each class the
    % fraction of sediment left after the timestep has passed
    setplace = zeros([size(v_sed_day,1) length(downdist) ]);

    %if there are whole reaches covered by the sediment volume between reach_s and
    %reach_e, mark them with -1
    setplace( and(downdist>downdist(reach_s)', downdist<downdist(reach_e)')) = -1;

    %i put the values of perc_s and perc_e in their respective places
    setplace(sub2ind(size(setplace) , 1:size(v_sed_day,1) ,  reach_s' )) = perc_s;
    setplace(sub2ind(size(setplace) , 1:size(v_sed_day,1) ,  reach_e' )) =  perc_e + setplace(sub2ind(size(setplace) , 1:size(v_sed_day,1) ,  reach_e' ))'; % i do the sum in case reach_s and reach_e are the same, to not remove the values. In most of the cases the second term will be 0

    %if there are whole reaches covered by the sediment volume between reach_s
    %and reach_e...
    if any(setplace==-1,'all')
        % ... attribute them a fraction of the remaining sed volume according
        % to their length
        [rows_find , ~] = find(setplace==-1);

        for indx_discrepancy=1:length(rows_find)

            setrow = setplace(rows_find(indx_discrepancy),:);

            setrow(setrow==-1) = (1 - perc_s(rows_find(indx_discrepancy)) - perc_e(rows_find(indx_discrepancy)) - perc_out(rows_find(indx_discrepancy)) ) * Lngt(setrow==-1)  ./ sum(Lngt(setrow==-1)) ;

            setplace(rows_find(indx_discrepancy),:) = setrow;
        end

    end

%% change setplace to account for nonuniform distribution of sediment due to changing velocity
    
    res_time = 1./v_sed_day;

    % i do a simple weighted average using the residence time as weigth, i
    % attibute to the leaving sediment the residence time of the outlet reach
    setplace = (setplace .* res_time)./sum( [setplace .* res_time , perc_out.* res_time(:,outlet) ] ,2);
    perc_out = (perc_out .* res_time(:,outlet))./sum( [setplace .* res_time , perc_out.* res_time(:,outlet) ] ,2);
    
    % if a row in setplace is all empty (thus i have nans as a result of
    % the division), all sediment have left the network 
    setplace(isnan(setplace)) = 0; 
  
%% save results

    % i save the results in matrix Qbi_tr_t and Q_out_t
    Qbi_tr_t = zeros (length(Lngt), length(Lngt) , size(setplace,1));
    
    for c = 1:size(setplace,1)
        
        if ~isnan(roundpar)
            Qbi_tr_t(V_mob(:,1),:,c) = floor( V_mob(:,c+1) .* setplace(c,:)*roundpar*10)/(roundpar*10);          
            %Qbi_tr_t(V_mob(:,1),:,c) = round( V_mob(:,c+1) .* setplace(c,:), roundpar);            
        else
            Qbi_tr_t(V_mob(:,1),:,c) = V_mob(:,c+1) .* setplace(c,:);
        end
 
    end
    
    Q_out_t = zeros (length(Lngt), size(setplace,1));
    
    if ~isnan(roundpar)
        Q_out_t(V_mob(:,1),:) = floor( V_mob(:,2:end) .* perc_out'*roundpar*10)/(roundpar*10);          
        %Q_out_t(V_mob(:,1),:) = round( V_mob(:,2:end) .* perc_out' , roundpar);
    else
        Q_out_t(V_mob(:,1),:) = V_mob(:,2:end) .* perc_out';
    end

    if ~isnan(roundpar)   
        for c = 1:size(setplace,1)    
            %if the rounding lead to removed or added volume from V_mob to Qbi_tr_t...
            if any (round(sum(V_mob(:,c+1) ,2) - sum(Q_out_t(V_mob(:,1),c),2) - sum(Qbi_tr_t(V_mob(:,1),:,c),2),roundpar*3)~=0)
                
                diff = round( sum(V_mob(:,c+1) ,2) -  sum(Q_out_t(V_mob(:,1),c),2) - sum(Qbi_tr_t(V_mob(:,1),:,c),2),roundpar*3);
                indx_discrepancy = find(diff~=0);

                for i = 1:length(indx_discrepancy) %for each reach in V_mob where I have displacement errors
                    
                   % find the position of the reaches with the highest delivered volume in the network
%                    max_reach = find(Qbi_tr_t(V_mob(indx_discrepancy(i),1),:,c) == max(Qbi_tr_t(V_mob(indx_discrepancy(i),1),:,c),[],2),1 );
%                    idx = sub2ind(size(Qbi_tr_t),V_mob(indx_discrepancy(i),1),max_reach,c) ; 
                   
                   % find the position of the reaches with the highest delivered volume in the network
                   [~,max_reach] = max([setplace(c,:) perc_out(c)]); % set the outlet as a new reach at the tail of the vector
                   if max_reach == length(Lngt)+1 % if the highest position is the outlet (set as the last reach) i add the approx value to the outvolume
                       Q_out_t(V_mob(indx_discrepancy(i),1),c) = Q_out_t(V_mob(indx_discrepancy(i),1),c) + diff(indx_discrepancy(i));
                   else
                       idx = sub2ind(size(Qbi_tr_t),V_mob(indx_discrepancy(i),1),max_reach,c) ; 
                       %add the volume in diff to the reach max_reach with the highest volume delivered.
                       Qbi_tr_t(idx) = Qbi_tr_t(idx) + diff(indx_discrepancy(i)); %add the volume in diff to the reach max_reach with the highest volume delivered.
                       % remove negative values arising from errors in the difference
                       if Qbi_tr_t(idx)<0
                            Qbi_tr_t(idx) = 0;
                       end
                   end
                   

%                    if Qbi_tr_t(idx)<0
%                        
%                        diff = diff-Qbi_tr_t(idx);
%                        Qbi_tr_t(idx) = 0;
%                        max_reach = find(Qbi_tr_t(V_mob(indx_discrepancy(i),1),:,c) == max(Qbi_tr_t(V_mob(indx_discrepancy(i),1),:,c),[],2),1 );
%                        idx = sub2ind(size(Qbi_tr_t),V_mob(indx_discrepancy(i),1),max_reach,c) ; % find the position of the reaches with the highest delivered volume in the network
%                    end
                end
            end
        end
    end
    
end

