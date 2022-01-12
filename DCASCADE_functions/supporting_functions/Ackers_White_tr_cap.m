function [ Qtr_cap , Pci ] = Ackers_White_tr_cap( Fi_r_reach,  Slope , Q , v, h)

%ACKERS_WHITE_TR_CAP returns the value of the transport capacity (in Kg/s) for each
%sediment class in the reach measured using the Ackers-White equations

%% references

% Ackers P., White W.R. Sediment transport: New approach and analysis (1973)

%% variables initialization

rho_w = 1000; %water density
rho_s = 2650; %sediment density
g = 9.81; %gravity acceleration

R = rho_s / rho_w - 1; %submerged specific gravity of sediment 
%FR = v/sqrt(g*h);     %Froude number

nu = 1.02*1E-6; % kinematic viscosity @ 20�C: http://onlinelibrary.wiley.com/doi/10.1002/9781118131473.app3/pdf
%nu = 0.000011337;  % kinematic viscosity (ft2/s)

%alpha = 10; %coefficient in the rough turbulent equation with a value of 10;
alpha = 7.18697;

%% compute changes in D16, D50, D85 along the river (without par)
    
D_values = [35 16 50 84];

[D_changes] = D_finder(Fi_r_reach, D_values );

% Ackers - White suggest to use the D35 instead of the D50, in m
D_AW = D_changes(1,:);

%% Ackers_White parameters

D_gr = D_AW * ( g * R / nu^2 )^(1/3); %dimensionless grain size

C = (D_gr < 60) .* ( 10 .^ ( 2.79 .* log10(D_gr) - 0.98 .* log10(D_gr).^2 - 3.46 ))  +  (D_gr >= 60).* 0.025 ;
m = (D_gr < 60) .* ( 6.83 ./ D_gr + 1.67 )  +  (D_gr >= 60) .* 1.78 ;        %alterative formula: m = 9.66 / D_gr + 1.34;
A = (D_gr < 60) .* ( 0.23./ sqrt(D_gr) + 0.14 )  +  (D_gr >= 60).* 0.17 ;
n = (D_gr < 60) .* (1 - 0.56 .* log10(D_gr)) +  (D_gr >= 60) .* 0 ;

%shear velocity
u_ast = sqrt(g .* h .* Slope);

%% Transport capacity 

% mobility factor
F_gr = u_ast.^n ./ sqrt(g .* D_AW .* R) .* ( v ./ (sqrt(32) * log10(alpha .* h ./ D_AW ) ) ) .^(1-n);
% alterative formula: F_gr = u_ast / sqrt(g * D_AW * R) * ( v / (sqrt(32) * log10( h /D_AW ) ) ) ;

% dimensionless transport
G_gr = C .* ( max(F_gr ./ A - 1 , 0) ) .^m ;

% weight concentration of bed material (Kg_sed / Kg_water)
QS_conc = G_gr .* (R + 1) .* D_AW ./ ( (u_ast ./ v ).^n .* h);

% transport capacity (Kg_sed / s)
QS_AW = rho_w .* Q .* QS_conc;
%QS_AW_td = QS_AW/1000*60*60*24 ;

%Pci contains the fraction rates for each sediment class
Pci = Molinas(Fi_r_reach, h, v, Slope , D_changes(2:end,:));

% transport capacity for each class (m3 / s)
Qtr_cap = Pci.*QS_AW ./ rho_s;

end

