function [ v_sed, Qtr_cap ] = velocity_AW( Fi_r_reach ,  Slope , Wac , Q , v , h , min_velocity)
%VELOCITY_AW returns the velocity of the sediment (in m/s) for each sediment
%class for each reach using the Ackers and White equations
%
%OUTPUTS:
% 
% v_sed: [cxn] matrix reporting the veloctivy for each sediment class c fro
%        each reach n
%% references
% Ackers P., White W.R. Sediment transport: New approach and analysis (1973)

%% optional input selection

if nargin < 5
    min_velocity = 0;
end

h_act = 0.1; %characteristic vertical length scale for transport.

%% variables initialization

[ Qtr_cap, pci] = Ackers_White_tr_cap( Fi_r_reach,  Slope , Q, v, h) ; 

%% calculate velocity

v_sed = max( Qtr_cap./( Wac .* h_act .* pci ) , min_velocity);

end