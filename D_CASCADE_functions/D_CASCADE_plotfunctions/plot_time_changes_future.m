function plot_time_changes_future(data_plot , data_plot_scenarios, reaches , varargin)
%PLOT_SNAP_NETWORK plots different subplots of the network showing the
%changes in [Txn] matrix data with time

%% default settings 

def_linewidth = 2;
def_cMap = 'hsv';
def_startT = 2;
def_futureT = size(data_plot_scenarios{1,1}{1,2},1);

%% read additional inputs 

p = inputParser;
addOptional(p,'cMap',def_cMap);
addOptional(p,'colors', []);
addOptional(p,'startT',def_startT);
addOptional(p,'futureT',def_futureT);
addOptional(p,'LineWidth',def_linewidth);
addOptional(p,'FontSize',12);
addOptional(p,'titles',[]);
addOptional(p,'reachID',[]);
addOptional(p,'plotID',[ 1 2 ]);
addOptional(p,'equal_axis', 0);
addOptional(p,'showlegend', 'on');
addOptional(p,'prc_plot', [0:10:100]);

parse(p,varargin{:})

linewidth = p.Results.LineWidth;
startT = p.Results.startT;
futureT = p.Results.futureT;
name_colormap = p.Results.cMap ;
colors = p.Results.colors ;
Fsize = p.Results.FontSize;
titles_subplot = p.Results.titles;
reach_id = p.Results.reachID;
plot_id = p.Results.plotID;
equal_axis = p.Results.equal_axis;
showlegend = p.Results.showlegend;
prc_plot = p.Results.prc_plot;

%% find constant colorbar for plot

cMapName = [name_colormap '(' num2str(length(data_plot(:,1))) ')'];
cMap = eval(cMapName);  

color1 = cMap(plot_id(1),:);
color2 = cMap(plot_id(2),:);

if ~isempty(colors)
    
    color1 = colors(1,:);
    color2 = colors(2,:);
  
end

%% find subplot design
if length(reaches) > 12
    n_subplot = [4 4];
elseif length(reaches) > 9
    n_subplot = [3 4];
elseif length(reaches) > 6
    n_subplot = [3 3];
elseif length(reaches) > 4
    n_subplot = [3 2];
elseif length(reaches) > 2
    n_subplot = [2 2];
elseif length(reaches) > 1
    n_subplot = [1 2];
else
    n_subplot = [1 1];
end

%% calculate time length of past and future sims and the number of future scenarios

t_hist = min(cellfun(@(x)size(x,1),data_plot(1:end-1,2)));
n_scenarios = size(data_plot_scenarios,1);

%% for each reach, for each future timestep extract the values in data_plot_scenarios for the 2data considered
future_scenarios_data = cell(length(reaches),1);

for n = 1:length(reaches)
    
    r = reaches(n);
    future_scenarios_data{n} = cell(2,1);
    
    for d = 1:length(plot_id)
        future_scenarios_data{n}{d} = cell2mat(cellfun( @(x)x{plot_id(d),2}(:,r) , data_plot_scenarios(:,1),'UniformOutput',0 )' );
    end
end

%% set y axis limits, if Equal axis option is decided

if equal_axis 
    yinterval1 = zeros(1,2);
    yinterval1(1) = min( [ data_plot{plot_id(1),2}(startT:end,reaches ); cell2mat( cellfun(@(x)min(x{1},[],2), future_scenarios_data,'UniformOutput',0 )') ]  ,[], 'all');
    yinterval1(2) = max( [ data_plot{plot_id(1),2}(startT:end,reaches ); cell2mat( cellfun(@(x)max(x{1},[],2), future_scenarios_data,'UniformOutput',0 )') ]  ,[], 'all');

    yinterval2 = zeros(1,2);
    yinterval2(1) = min( [ data_plot{plot_id(2),2}(startT:end,reaches ); cell2mat( cellfun(@(x)min(x{2},[],2), future_scenarios_data,'UniformOutput',0 )') ]  ,[], 'all');
    yinterval2(2) = max( [ data_plot{plot_id(2),2}(startT:end,reaches ); cell2mat( cellfun(@(x)max(x{2},[],2), future_scenarios_data,'UniformOutput',0 )') ]  ,[], 'all');
end

%% plot data in frame
fig = figure;
set(gcf,'color','w');

set (fig, 'Units', 'normalized', 'Position', [0.1,0.1,0.8,0.8]);
for i=1:length(reaches)
    subplot(n_subplot(1),n_subplot(2),i);
    
    if i==length(reaches) && strcmp(showlegend,'on') %only for last subplot...
        %plot "fake" lines for legend
        plot(yinterval1(1),'Color',color1, 'LineWidth', 4)
        hold on
        plot(yinterval1(1),'Color',color2 ,'LineWidth', 4)
        hold on
    end
    
    data_allscenario = [ repmat( data_plot{plot_id(1),2}(:,reaches(i)),[1 n_scenarios]) ; future_scenarios_data{i}{1} ];
    fanChart(startT :t_hist + futureT , data_allscenario (startT :t_hist + futureT,:),[], prc_plot ,'cCenter',color1.*0.5,'colormap',{'shadesOfColor', color1 } );
    if equal_axis; ylim( yinterval1 ); end
    
    yyaxis right
    hold on
    data_allscenario = [ repmat( data_plot{plot_id(2),2}(:,reaches(i)),[1 n_scenarios]) ; future_scenarios_data{i}{2} ];
    fanChart(startT :t_hist + futureT , data_allscenario (startT :t_hist + futureT,:) ,[],prc_plot,'cCenter', color2.*0.5 ,'colormap',{'shadesOfColor', color2 } );
    
    if equal_axis; ylim( yinterval2 ); end
           
    %show title on the subplot
    if isempty(reach_id) %Set different title form reach_id if specified
        reach_name = num2str(reaches(i));
    else
        reach_name =  num2str(reach_id(i));
    end
    
    if isempty(titles_subplot) %add detail to the title if specified
        title( [ 'Reach ' reach_name] );
    else      
        title( [ titles_subplot{i}] );
    end

    hold off
    
    xlim([startT t_hist + futureT])

    grid minor
    set(gca,'XGrid','on','FontSize',Fsize)

end

%% set legend
if strcmp(showlegend,'on')
    legend( data_plot(plot_id,1),'FontSize',15,'location', 'southeast')
end

end
