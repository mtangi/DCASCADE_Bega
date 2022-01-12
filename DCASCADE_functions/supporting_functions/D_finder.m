function [D_changes] = D_finder(Fi_r, D_values )
%D_FINDER find the value of granulometry for the specified D_values, for
%the sediment distribution Fi_r

%% Define inputs
%if D values not given, find only the D50
if nargin == 1
    D_values =  50 ;
end

%%

global psi
dmi= 2.^(-psi)./1000;

D_changes = zeros(length(D_values),size(Fi_r,2));
Perc_finer=zeros(length(dmi),size(Fi_r,2));
Perc_finer(1,:)=100;

for i=2:size(Perc_finer,1)
    Perc_finer(i,:) = Perc_finer(i-1,:) - (Fi_r(i-1,:)*100);  
end

for k=1:size(Perc_finer,2)
   for j=1:length(D_values)
       a = min( find( Perc_finer(:,k) >  D_values(j), 1, 'last' ),length (psi)-1);
       D_changes(j,k) = (D_values(j) - Perc_finer(a+1,k))/(Perc_finer(a,k) - Perc_finer(a+1,k))*(-psi(a)+psi(a+1))-psi(a+1);
       D_changes(j,k) =  2^(D_changes(j,k))./1000;
       D_changes(j,k) = D_changes(j,k)*(D_changes(j,k)>0) + dmi(end)*(D_changes(j,k)<0);
   end
end

%% only for D50 

%D_changes = 2.^sum(-psi.* Fi_r_act(:,n)') /1000;

end

