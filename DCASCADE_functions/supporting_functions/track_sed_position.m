function [sed_pos , end_reach_ID, outind] = track_sed_position( n , v_sed_day , Lngt , Network , start_pos)
%TRACK_SED_POSITION_TRCAP finds the position of a sediment parcel starting
%from reach n after the timestep has passed, defined as the reach ID and
%the position from the From_node of the starting reach. 
%
%To satisfy the transport capacity in the ToNode section, the starting
%position of the volume is positioned in a location that guarantees that
%all of it passes through the ToNode and leaves the reach

%% define starting position
% the starting position is the position on the reach n from which the
% parcel start, defined as fraction of reach length
% if start_pos = 0, the parcel starts form the From_Node
% if start_pos = 1, the parcel starts form the To_Node

global psi

if nargin < 5 
    start_pos = 0;
end


%% find path downstream

% start_pos (between 0 and 1) defines the position in the reach n where the 
% sediment parcel start, if 1 start form the From_node, if 0 starts from
% the To_Node
% if nargin < 5
%     start_pos = 1;
% end

timestep = 1;

outlet = Network.NH(end);

%path2out contains the path from reach n to the outlet, defined as the IDs of the
%reach downstream ordered.
%downdist_path contains the distance from the reaches in path2out to the
%reach n From_node
path2out = cell2mat(Network.Downstream.Path{n,1}(1,outlet));
downdist_path = Network.Downstream.Distance{n,1}(path2out) ;

%
%% find position and destination reach ID

%isolate the velocity of the downstream reaches
v_sed_path = v_sed_day(:,path2out);

%change the length of the starting reach according to the starting
%position, different for each tr.cap
Lngt_pathout = repmat(Lngt(path2out),length(psi),1);
Lngt_pathout(:,1)  = Lngt_pathout(:,1) .* (1 - start_pos); 

%calculate the time (in days) it takes to completely cross a reach
transit_time = Lngt_pathout./v_sed_path;

% the cumulate of transit_time defines how long it takes to reach each
% downstream To_Node comnsidering the whole path to the reach
cum_tr_time = cumsum(transit_time,2);

%given cum_tr_time, i can find the reach where the parcel is
%after the timestep  
find_point = cum_tr_time - timestep;
find_point(find_point<0) = nan;
[~,indx_pos] = min(find_point,[],2); %(note, this is not the reach ID, it is the position of the reach in the downstream path)

indx_pos = indx_pos';
indx_pos(isnan(find_point(:,end))) = length(path2out); %if the whole row in find_point is nan, it means that the parcel left the network

%I can find the time remaining for the parcel after if enters reach
%indx_pos, needed to find the position of the parcel
find_time = timestep - cum_tr_time;
find_time(find_time<0) = nan;
[~,indx_t] = min(find_time,[],2); %indx_t is the reach before indx_pos 

time_left = find_time( sub2ind( size(find_time) , 1:size(find_time,1) , indx_t'  ) ) ;
time_left(isnan(time_left)) = timestep ; %if time_left is nan, it means that the parcel remained in the starting reach m

if n == outlet

    %the if is necessary since the velocities are extracted as a column
    %vector, and need to be transposed
    sed_pos = time_left' .* v_sed_path( sub2ind(size(transit_time) , 1:length(psi) ,indx_pos ) )' + downdist_path(indx_pos);
    sed_pos(isnan(find_point(:,end))) = downdist_path(length(path2out)) + Lngt(outlet) + v_sed_path(isnan(find_point(:,end)), length(path2out)) .* time_left(isnan(find_point(:,end))) ;
      
else
  
    %find the position of the parcel, given the length of the reaches
    %already passed and the remaining time and velocity in the last
    %reach(indx_pos).

    sed_pos = time_left .* v_sed_path( sub2ind(size(transit_time) , 1:length(psi) ,indx_pos ) ) + downdist_path(indx_pos);

    % If the whole row in find_point is nan (the parcel left the network),
    % use the velocity of the outlet to determine the final position
    % (that will be outside the network)

    sed_pos(isnan(find_point(:,end))) = downdist_path(length(path2out)) + Lngt(outlet) + ...
         v_sed_path(isnan(find_point(:,end)), length(path2out))' .* time_left(isnan(find_point(:,end))) ;
     
end

%outind tells for which sed. size the point fell outside the
%network (1 - outside, 0 - inside)
outind = isnan(find_point(:,end));

sed_pos = sed_pos';

%sed_pos = sed_pos + Lngt(n) * (1 - start_pos);
end_reach_ID = path2out(indx_pos)'; % i find the ID of the destination reach from indx_pos, given the order defined by path2out 

end

