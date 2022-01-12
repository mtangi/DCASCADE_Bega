function [ Pci ] = Molinas( Fi_r_reach, h, v, Slope, D_changes)
% MOLINAS returns the Molinas coefficient of fractional transport rates Pci, to be multiplied
% by the total sediment load to split it into different classes (function adapted for dynamic mode) 

%% references

% Molinas, A., & Wu, B. (2000). Comparison of fractional bed-material load computation methods in sand?bed channels. Earth Surface Processes and Landforms: The Journal of the British Geomorphological Research Group

%% Molinas and wu coefficients 

%sediment classes diameter (mm)
global psi 
dmi = 2.^(-psi)./1000;
dmi_finer= 2.^(-psi); 

%convert sediment classes diameter (m to mm)
D16 = D_changes(1,:).*1000;
D50 = D_changes(2,:).*1000;
D84 = D_changes(3,:).*1000;

% Molinas requires D50 and dmi in mm
g = 9.81; 
   
%geometric standard deviation
std = sqrt(D84 ./ D16);

% Hydraulic parameters in each flow percentile for the current reach
Dn = ( 1 + (std-1).^1.5 ) .* D50; %scaling size of bed material

tau = 1000*9.81*h.*Slope;
vstar = sqrt(tau / 1000);
FR = v ./ sqrt( g * h );     %Froude number

% alpha, beta, and Zeta parameter for each flow percentile (columns), and each grain size (rows)
% % EQ 24 , 25 , 26 , Molinas and Wu (2000)  
alpha = - 2.9 * exp(-1000*(v./vstar).^2 .* (h./D50).^(-2));
beta = 0.2 .* std;
Zeta = 2.8 .*FR.^(-1.2) .* std.^(-3); 
Zeta(isinf(Zeta)) = 0; % Zeta gets inf when there is only a single grain size. 

% alpha, beta, and Zeta parameter for each flow percentile (columns), and each grain size (rows)
% EQ 17 , 18 , 19 , Molinas and Wu (2003)  
% alpha = - 2.85* exp(-1000*(v./vstar).^2.*(h./D50).^(-2));
% beta = 0.2.* GSD_std(Fi_r_reach,dmi);
% Zeta = 2.16.*FR.^(-1);
% Zeta(isinf(Zeta)) = 0; 

% fractioning factor for each reach (columns), and each grain size (rows) 
frac = Fi_r_reach .* ( (dmi_finer'./Dn).^(alpha) + Zeta.*(dmi_finer'./Dn).^(beta) ) ; % Nominator in EQ 23, Molinas and Wu (2000) 
Pci = frac./sum(frac,1);

end

