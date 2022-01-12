%% Figure 5b) - plot sediment volumes and channel width for the lowland reaches, for one simulation 

id_sim = 2;
%1 - highQ
%2 - medQ
%3 - lowQ
%4 - mixQ

reaches =  [86 90 92 95 98 52 55 58]; %define ID of the plotted reaches
titles_subplot = {'A';'B';'C';'D';'E';'F';'G';'H'};  %define name of the plotted reaches

timeinterval = [50,218]; %define time boundaries of the plot
colorplot = [0.6 0.9 0.6;0.9 0.6 0.6; 1 0 0];

fontsize = 16;
ID2ndData = 1; % the second data plotted is the Width (based on the order in data_plot)

color_light = [1 1 1 ];
color_dark  = [0.4 0.4 1 ];
cMap_bcg = [linspace(color_light(1),color_dark(1),size(Q_ARR_reach,2))', linspace(color_light(2),color_dark(2),size(Q_ARR_reach,2))', linspace(color_light(3),color_dark(3),size(Q_ARR_reach,2))'];      
    
plot_time_changes_provenance(data_plot_multiple_compacted{id_sim,1}, reaches,'equal_axis',1,'ID_prov',2,'ID2nddata',ID2ndData,'background',8,'titles', titles_subplot,'reachID',[ReachData(reaches).ID_Bega],'Q_ARR_reach',Q_ARR_reach,'showlegend','off','colorplot',colorplot,'cMap_bcg',cMap_bcg,'timeinterval',timeinterval,'Fontsize', fontsize)

for i=1:length(reaches)
    sb = subplot(3,3,i);   
    set(sb, 'XTickLabels', { '1850','1900','1950','2000'})
    ax = gca;
    ax.YAxis(2).Color = [0.5 0 0];

    ylim([0 180])

    if any(i == [1,4,7])
        yyaxis left
        ylabel('Sediment volume [m^3]')
    end
    if any(i == [3,6,8])
        yyaxis right
        ylabel('Width [m]')

    end

end

sb = subplot(3,3,9);     
cMap = hsv (length(data_plot_multiple{1,1}(:,1)));  

a1=area(1, 1 , 'FaceColor', colorplot(1,:) , 'LineWidth' , 0.001);
hold on
a2=area([0], 0  , 'FaceColor', colorplot(2 ,:) , 'LineWidth' , 0.001);
a3=plot([0] , 'Color', colorplot(3,:) , 'LineWidth' , 2);
ylim([100,200])
xlim([10 100])
legend( { 'Sediment volume [m^3] from type 2-3' ,'Sediment volume [m^3] from type 4' , 'Channel width [m]' } ,'FontSize',15,'location', 'northwest')

colormap(cMap_bcg);
clb = colorbar('south','YTick',[1/7/2:1/7:1],'YTickLabel', { 'Q1','Q2','Q5','Q10','Q20','Q50','Q100'},'FontSize',11);
set(sb, 'Color', 'None')
set(sb,'visible','off')

title(clb,'Discharge return period (3 timestep average)','FontSize',15)

%% Figure 5c) - plot sediment volumes and channel width for the type 4 reaches, for one simulation 

id_sim = 2;
%1 - highQ
%2 - medQ
%3 - lowQ
%4 - mixQ

ID2ndData = 1; % the second data plotted is the Width (based on the order in data_plot)
reaches =  sort([187;188;211;215;233;236;248;262]);
[~,ord] = sort([ReachData(reaches).ID_Bega]);

reaches= reaches(ord);
%
colorplot = [0.9 0.6 0.6;1 0 0];
    
plot_time_changes_2data(data_plot_multiple_compacted{id_sim,1}, reaches, 'plotID', [3 ID2ndData] ,'equal_axis',1,'background',8,'reachID',[ReachData(reaches).ID_Bega],'Q_ARR_reach',Q_ARR_reach,'showlegend','off','colorplot',colorplot,'title', {ReachData(reaches).river_Name},'timeinterval',timeinterval,'cMap_bcg',cMap_bcg,'Fontsize', fontsize)

for i=1:length(reaches)
    sb = subplot(3,3,i);   
    set(sb, 'XTickLabels', { '1850','1900','1950','2000'})
    ax = gca;
    ax.YAxis(2).Color = [0.5 0 0];
    ylim([0 180])

     if any(i == [1,4,7])
        yyaxis left
        ylabel('Sediment volume [m^3]')
    end
    if any(i == [3,6,8])
        yyaxis right
        ylabel('Width [m]')

    end
    
end

sb = subplot(3,3,9);     
cMap = hsv (length(data_plot_multiple{1,1}(:,1)));  

a1=area(1, 1 , 'FaceColor', colorplot(1  ,:) , 'LineWidth' , 0.001);
hold on
a3=plot([0] , 'Color', colorplot(2,:) , 'LineWidth' , 2);
ylim([100,200])
xlim([10 100])
legend( { 'Sediment volume [m^3]' , 'Channel width [m]' } ,'FontSize',15,'location', 'northwest')

colormap(cMap_bcg);
clb = colorbar('south','YTick',[1/7/2:1/7:1],'YTickLabel', { 'Q1','Q2','Q5','Q10','Q20','Q50','Q100'},'FontSize',11);
set(sb, 'Color', 'None')
set(sb,'visible','off')

title(clb,'Discharge return period (3 timestep average)','FontSize',15)

%% Figure 6a) - plot a comparision of the sediment volumes and channel width in the lowland sections for different discharge scenarios

reaches =  [86 90 92 95 98 52 55 58];

plot_id = [ 3 1 ];
    
%plot results
plot_time_changes_multiple(data_plot_multiple_compacted, reaches , 'plotID', plot_id ,'equal_axis',1,'background',8,'titles', titles_subplot,'reachID',[ReachData(reaches).ID_Bega],'showlegend','off','Q_ARR_reach',Q_ARR_reach,'cMap_bcg',cMap_bcg,'timeinterval',timeinterval,'Fontsize', fontsize)

%add labels, captions and legend
for i=1:length(reaches)
    sb = subplot(3,3,i); 
    set(sb, 'XTickLabels', { '1850','1900','1950','2000'})
    ax = gca;
    ax.YAxis(2).Color = [0.1,0.4,0];
    ylim([0 150])

    if any(i == [1,4,7])
        yyaxis left
        ylabel('Sediment volume [m^3]')
    end
    if any(i == [3,6,8])
        yyaxis right
        ylabel('Width [m]')

    end
    
end

sb = subplot(3,3,9);  
name_cMap1 = 'copper'; name_cMap2 = 'summer';
cMapName1 = [name_cMap1 '(' num2str(size(data_plot_multiple,1)) ')'];
cMap1 = eval(cMapName1);  

cMapName2 = [name_cMap2 '(' num2str(size(data_plot_multiple,1)) ')'];
cMap2 = eval(cMapName2);
cMap = [cMap1;cMap2];

for i=1:length(data_plot_multiple_compacted)*2 %plot fake lines for legend

    plot(1 , 'Color', cMap(i,:) , 'LineWidth' , 2) 
    hold on
end

leg = legend([data_plot_multiple(:,2); data_plot_multiple(:,2)] ,'FontSize',15,'location', 'southeast','NumColumns',2);
str = {['   Width [m]         Sed.volume [m^3]']};
title(leg,str)

colormap(cMap_bcg);
clb = colorbar('south','YTick',[1/7/2:1/7:1],'YTickLabel', { 'Q1','Q2','Q5','Q10','Q20','Q50','Q100'},'FontSize',11);
set(sb, 'Color', 'None')
set(sb,'visible','off')

title(clb,'Discharge return period','FontSize',15)

%% Figure 6b) - plot scenarios perfomances in a parallel plot

plot_attributes = [16:21];
matres = [];
fns = fieldnames(SimVal);
for i = 1:length(plot_attributes)
    matres = [matres; SimVal.(fns{plot_attributes(i)})];
end
matres = matres';

labels =  {'SedDelRt - C';'SedDelRt - E';'SedDelRt - H';'Type 4 - % incised';'Type 3 - Width';'Type 3 - D50'};

GrpData = categorical({' SCENARIOS ';'HighQ';'MedQ';'LowQ';'MixQ';'Present-day data';'pre-ES conditions'}); 

matres = [[ 100 100 100 1 132 1.4 ]; matres];

figure
p = parallelplot(matres,'GroupData',GrpData','Jitter',0); %,'Labels',labels,'Group',{SimVal.simname},'Standardize','on')
p.LineWidth = 4;
p.Color = [ 1 1 1 ; hsv(length(GrpData)-3); 0 0 0 ; 0.7 0.7 0.7 ];

p.CoordinateTickLabels = [labels];
set(gcf,'color','w');
set(gca,'FontSize',16)


%% Figure 7a) - plot future trajectories in the Lowland sections - with 2 scenarios (no restoration and restoration)

%note: to obtain this plot, we need to run the DCASCADE_Begafuturesim.m
%script two times, to obtain outputs for both scenarios.
%Here we assume the data_plot_scenarios_compacted variable contains the
%outputs for the restoration scenario with present-day vegetation,
%while data_plot_scenarios_compacted_nores contains the outputs for the
%pre-restoration vegetation scenario.  

reaches = [86 90 92 95 98 52 55 58];
color_plot = [0 0 1;1 0 0];

titles_subplot = {'A';'B';'C';'D';'E';'F';'G';'H'};  %define name of the plotted reaches

fontsize = 16;
ID2ndData = 1; % the second data plotted is the Width (based on the order in data_plot)

id_sim = 4;

prc_plot = [0:10:100];

data_plot_scenario_mix = data_plot_scenarios_compacted;

plot_id =  [7 3];

for i=1:length(data_plot_scenario_mix)
    
    data_plot_scenario_mix{i,1}{7,2} = data_plot_scenarios_compacted_nores{i,1}{3,2};
    data_plot_scenario_mix{i,1}{7,1} = 'Total sed in the reach - no_restoration scenario';
end

data_plot_multiple_mix = data_plot_multiple_compacted{id_sim,1};
data_plot_multiple_mix{7,2} = data_plot_multiple_mix{3,2}; 

%plot future scenario results
plot_time_changes_future(data_plot_multiple_mix, data_plot_scenario_mix,reaches, 'plotID', plot_id ,'futureT',100,'equal_axis',1,'startT',50,'colors', color_plot,'titles', titles_subplot,'reachID',[ReachData(reaches).ID_Bega],'showlegend','off','prc_plot',prc_plot,'Fontsize', fontsize)

%add labels, captions and legend
color_default = [0.9 0.9 0.9];
cMap_future = shadesOfColor(color_default, floor(length(prc_plot)/2));
cMap_future(:,1) = cMap_future(:,2);

for i=1:8
    sb = subplot(3,3,i); 
    set(sb, 'XTickLabels', { '1850','1900','1950','2000','2050','2100'})
    ax = gca;
    ax.YAxis(2).Color = [ 0 0 0 ];
    ax.YAxis(2).Limits = ax.YAxis(1).Limits;
    set(gca,'ytick',[])

    if any(i == [1,4,7])
        yyaxis left
        ylabel('Sediment volume [m^3]')
    end

end

sb = subplot(3,3,9);  
for i=1:size(color_plot,1) %plot fake lines for legend
    plot(1 , 'Color', color_plot(i,:).*0.5 , 'LineWidth' , 2) 
    hold on
end

colormap(flip(cMap_future));
clb = colorbar('south','YTick',[0:1/(size(cMap_future,1)):1],'YTickLabel', num2cell(flip(prc_plot(1:2:end))) ,'FontSize',11);
title(clb,'Percentile range','FontSize',fontsize)

legend( {'Sediment volume [m^3] - pre-recovery roughness','Sediment volume [m^3] - present-day roughness'} ,'FontSize',fontsize,'location', 'northwest')
set(sb, 'Color', 'None')
set(sb,'visible','off')

clear data_plot_scenario_mix clb data_plot_multiple_mix sb ax color_default cMap_future prc_plot data_plot_scenario_mix
 
%% Figure 7b) - plot future scenarios perfomances in a parallel plot
SimValopt1 = SimVal_future;
SimValopt2 = SimVal_future_nores;

reach_type = [ReachData.reach_type];
total_dep = [ReachData.deposit].*[ReachData.Length];

id_sim = 4;
lim_final = 218; %timestep to which the data on the reference scenario are defined
plot_attributes = [16:21];
matres = [];
fns = fieldnames(SimVal);
for i = 1:length(plot_attributes)
    matres = [matres; SimValopt1.(fns{plot_attributes(i)})];
end
matres = matres';

avgQ_future = cell2mat(cellfun(@(x)str2double(x(7:end)),{SimValopt1.simname}, 'UniformOutput', false))';

n=10;
prc = (max(avgQ_future) - min(avgQ_future))/n * [0:1:n] + min(avgQ_future) ;

GrpData = sum(avgQ_future >= prc,2);
GrpData = prc(GrpData);

[GrpData,b] = sort(GrpData);
matres = matres(b,:);

% find values of val parameters for the norestoration scenario
matres_nores = [];
fns = fieldnames(SimVal);
for i = 1:length(plot_attributes)
    matres_nores = [matres_nores; SimValopt2.(fns{plot_attributes(i)})];
end
matres_nores = matres_nores';

avgQ_future = cell2mat(cellfun(@(x)str2double(x(7:end)),{SimValopt2.simname}, 'UniformOutput', false))';

n=10;
prc = (max(avgQ_future) - min(avgQ_future))/n * [0:1:n] + min(avgQ_future) ;

GrpData = sum(avgQ_future >= prc,2);
GrpData = prc(GrpData);

[GrpData,b] = sort(GrpData);
matres_nores = matres_nores(b,:);

%data of the scenario used for comparison
data_id_SedDelRt = data_plot_multiple{id_sim,1}{9,2}(lim_final,SedDelRt_field(15:17,1));
data_id_swamp = 1 - sum(data_plot_multiple{id_sim,1}{3,2}(lim_final, reach_type == 4))/sum(total_dep(reach_type == 4));
data_id_D50 = mean(data_plot_multiple{id_sim,1}{5,2}(lim_final-6:lim_final , 94:99  ),'all');
data_id_width = mean(data_plot_multiple{id_sim,1}{1,2}(lim_final, 94:99  ),'all');

data_id_sim = [data_id_SedDelRt , data_id_swamp , data_id_width, data_id_D50];

%data to define maximum and minimum in the graph 
matres(end+1:end+3,:) = [[0,0,0, 0,42.5000000000000,0.455452025285418];[100,100,100,1,132.800000000000,1.41000000000000]; data_id_sim];
matres = [matres; matres_nores];
GrpData = categorical([round(GrpData,-2)+50, 0 , 0 , 0 ,round(GrpData,-2)]);
GrpData(1,103) = 'MixQ (2020)';
GrpData(1,101:102) = '. ';

matres = matres([104:end,1:100, 101:103],:);
GrpData = GrpData([104:end,1:100, 101:103]);

figure
p = parallelplot(matres,'GroupData',GrpData','Jitter',0,'LineWidth',2); %,'Labels',labels,'Group',{SimVal.simname},'Standardize','on')

labels=fns(plot_attributes);
labels =  {'SedDelRt - C';'SedDelRt - E';'SedDelRt - H';'Type 4 - % incised';'Type 3 - Width';'Type 3 - D50'};
p.CoordinateTickLabels = [labels];
p.LineWidth = [ones(23,1)*1.5; 2.5];
set(gcf,'color','w');
set(gca,'FontSize',16)

p.Color = [winter(length(prc)); copper(length(prc)); 1 1 1 ; 1 0 0; ];

clear data_id_SedDelRt data_id_swamp data_id_D50 data_id_width data_id_slope data_id_sim SimValopt1 SimValopt2 p GrpData matres

%% Figure S1a) bar plot of historic flood

% the matrix depth_hist contains the historic flood depth record in Bega Township, with the associated discharge if available 
figure
depth_hist = [1851,552.228056500000,9.78000000000000;1861,157.154454300000,5.90000000000000;1870,443.061730500000,9.10000000000000;1873,386.708472000000,8.68000000000000;1874,247.321499300000,7.30000000000000;1878,178.893222700000,6.30000000000000;1891,110.050473100000,4.80000000000000;1892,119.332546800000,5.05000000000000;1893,341.924210900000,8.30000000000000;1898,415.270081300000,8.90000000000000;1900,121.280871400000,5.10000000000000;1913,224.419957500000,7;1914,243.348386700000,7.25000000000000;1919,387.963051400000,8.69000000000000;1922,106.543047700000,4.70000000000000;1926,119.332546800000,5.05000000000000;1928,247.321499300000,7.30000000000000;1934,402.034982800000,8.80000000000000;1935,231.807922100000,7.10000000000000;1942,119.332546800000,5.05000000000000;1943,0,6.30000000000000;1949,0,6.05000000000000;1950,250.767500000000,6.15000000000000;1952,465.520910000000,7.20000000000000;1953,217.385230000000,6.10000000000000;1955,0,6.72000000000000;1956,345.679480000000,6.75000000000000;1958,0,6.52000000000000;1960,177.817650000000,5.95000000000000;1961,0,6.50000000000000;1962,212.772790000000,8.10000000000000;1963,0,5.20000000000000;1965,181.125280000000,5.95000000000000;1968,0,5.55000000000000;1969,87.5826500000000,5.50000000000000;1971,622.771720000000,9.75000000000000;1973,0,5.25000000000000;1974,0,6.40000000000000;1975,316.778380000000,7.30000000000000;1976,0,5;1977,0,6.24000000000000;1978,0,6.25000000000000;1979,0,6.28000000000000;1981,236.728820000000,7.10000000000000;1983,0,6.44000000000000;1984,0,6.49000000000000;1985,218.512360000000,7;1988,157.027520000000,6.60000000000000;1990,0,6.65000000000000;1991,0,5.26000000000000;1992,0,7.30000000000000;1993,0,7.10000000000000;1997,0,6.90000000000000;2010,143.252780000000,6.70000000000000;2011,264.346190000000,8.40000000000000;2012,95.5839600000000,6.50000000000000];

color = [0 0 0.5];
bar(depth_hist(:,1),depth_hist(:,3),0.2,'FaceColor',color,'EdgeColor',color);
ylim([4,10])
xlim([1850,2013])
set(gcf,'color','w');
xticks(1850:10:2020)
xtickangle(75)
set(gca, 'YGrid', 'on', 'XGrid', 'off','FontSize',16)
set(gca,'TickDir','out');
box off
ylabel('Gauge Heigth [m]')
xlabel('Year')
title('Flood heigth record at North Bega gauge','FontSize',20)
clear color

%% clear all variables
clear a1 a2 a3 ax p sb str timeinterval titles_subplot ord fns cMap colorplot labels matres id_sim ID2ndData fontsize GrpData i cMapName2 cMapName1 cMap1 cMap2 leg clb plot_id color_light cMap_bcg color_dark color_plot reaches ans name_cMap1 name_cMap2 plot_attributes

