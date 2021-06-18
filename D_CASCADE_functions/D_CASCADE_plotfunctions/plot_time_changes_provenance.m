function plot_time_changes_provenance(data_plot , reaches , varargin)
%PLOT_TIME_CHANGES_PROVENANCE plots different subplots of the network showing the
%changes in two data in the data_plot matrix with time.
%The first data visualized is always the total sediment volume in the
%reach, divided according to the provenance

%% default settings 

def_linewidth = 2;
def_cMap = 'hsv';
def_timeinterval = [2 min(cellfun('size',data_plot(1:8,2),1))];

%% read additional inputs 

p = inputParser;
addOptional(p,'cMap',def_cMap);
addOptional(p,'timeinterval',def_timeinterval);
addOptional(p,'LineWidth',def_linewidth);
addOptional(p,'FontSize',12);
addOptional(p,'titles',[]);
addOptional(p,'reachID',[]);
addOptional(p,'ID2ndData',[ 1 ]);
addOptional(p,'equal_axis', 0);
addOptional(p,'ID_prov', 1);
addOptional(p,'background', []);
addOptional(p,'Q_ARR_reach', []);
addOptional(p,'showlegend', 'on');
addOptional(p,'colorplot', [] );
addOptional(p,'cMap_bcg', [] );

parse(p,varargin{:})

linewidth = p.Results.LineWidth;
timeinterval = p.Results.timeinterval;
name_colormap = p.Results.cMap ;
Fsize = p.Results.FontSize;
titles_subplot = p.Results.titles;
reach_id = p.Results.reachID;
ID2ndData = p.Results.ID2ndData;
equal_axis = p.Results.equal_axis;
ID_prov = p.Results.ID_prov;
bcg_id = p.Results.background;
Q_ARR_reach = p.Results.Q_ARR_reach;
showlegend = p.Results.showlegend;
colorplot = p.Results.colorplot;
cMap_bcg = p.Results.cMap_bcg;
ID_tot_sed = 3;
ID_tot_sed_prov = 10;

%% find constant colorbar for plot

cMapName = [name_colormap '(' num2str(length(data_plot(:,1))) ')'];
cMap = eval(cMapName);  

if isempty(colorplot)
    colorplot = cMap([ID_tot_sed,ID_tot_sed_prov,ID2ndData],:);
end

%% find subplot design
if length(reaches) > 16
    n_subplot = [4 4];
    reaches = reaches(1:16);
    warning('max number of displayed reaches surpassed')
elseif length(reaches) > 12
    n_subplot = [4 4];    
elseif length(reaches) > 9
    n_subplot = [4 3];
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

%% define background colors (optional)
% if I have info on the RP of the floods events, i can plot bars in the
% backgroud showing flood intensity
if ~isempty(bcg_id)
    bcg_window = 2;
    if isempty(cMap_bcg)
        color_light = [0.9 1 1 ];
        color_dark  = [0.4 0.4 1 ];
     
        if isempty(Q_ARR_reach)
           cMap_bcg = [linspace(color_light(1),color_dark(1),101)', linspace(color_light(2),color_dark(2),101)', linspace(color_light(3),color_dark(3),101)'];
        else
           cMap_bcg = [linspace(color_light(1),color_dark(1),size(Q_ARR_reach,2))', linspace(color_light(2),color_dark(2),size(Q_ARR_reach,2))', linspace(color_light(3),color_dark(3),size(Q_ARR_reach,2))'];      
        end
    end
end

%% set y axis limits, if Equal axis option is decided

if equal_axis 
    yinterval1 = [ min( data_plot{ID_tot_sed,2}(:,reaches ) ,[], 'all') ; max(data_plot{ID_tot_sed,2}(:,reaches ) ,[], 'all')] ; 
    yinterval2 = [ min( data_plot{ID2ndData,2}(:,reaches ) ,[], 'all') ; max(data_plot{ID2ndData,2}(:,reaches ) ,[], 'all')] ; 
end

%% plot data in frame
fig = figure;
set(gcf,'color','w');

set (fig, 'Units', 'normalized', 'Position', [0.1,0.1,0.8,0.8]);

for i=1:length(reaches)
    subplot(n_subplot(1),n_subplot(2),i);
    
    % if i want the background coloration
    if ~isempty(bcg_id) && equal_axis
        lm = max(yinterval1(2),yinterval2(2)); %define the heigth of the bars, give by the limits of the plot
        M = movmean(data_plot{bcg_id,2}(:,i),bcg_window);  %I use if i want the moving average
        datacolor = round((M-min(M))/(max(M)-min(M)) *100 ) + 1; % attribute to each value in the timestep a class from 1 to 100 given their value
        
        if ~isempty(Q_ARR_reach)
            M = matlab.tall.movingWindow(@(x) max(x),bcg_window,data_plot{bcg_id,2}(:,i)); % this function find for each element the maximum in a window defined by bcg_window
            datacolor = sum(M >= Q_ARR_reach(i,:),2); 
        end %if the background is Q, use the color of the prcentiles in Q_ARR_reach
        
        if max(M) == min(M) %if there is no variation in the backgroud plot variable
            mapcolor = repmat(cMap_bcg(1,:),size(datacolor)); % i pur all colors as the ligth color
        else
            mapcolor = cMap_bcg(datacolor,:); % Attribute to each class a color in the colorscale cMap_bcg previously defined 
        end
        
        b = bar(repmat(lm ,[1 max(timeinterval)]),1,'EdgeColor', 'None'); %plot n bars, allof the same heigth given by lm, one for each timestep considered in timeinterval
        b.FaceColor = 'flat';
        b.CData = mapcolor; %color each bar according to the color given by mapcolor
        hold on
    end
      
    area(data_plot{ID_tot_sed,2}(:,reaches(i) ) , 'FaceColor', colorplot(1,:) , 'LineWidth' , linewidth)
    hold on
    area(data_plot{ID_tot_sed_prov,2}{ID_prov}(:,reaches(i) ) , 'FaceColor', colorplot(2,:) , 'LineWidth' , linewidth)
    
    if equal_axis; ylim( yinterval1 ); end
    
    yyaxis right
    hold on
    
    plot(data_plot{ID2ndData,2}(:,reaches(i)) , 'Color', colorplot(3,:) , 'LineWidth' , linewidth)
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
        title( [ titles_subplot{i}]);
    end
    hold off
    
    xlim(timeinterval)
    
    grid minor
    set(gca,'XGrid','on','FontSize',Fsize)

end

%% set legend
if strcmp(showlegend,'on')
    if ~isempty(bcg_id)
        legend( { ['bckground: ' data_plot{[bcg_id] , 1}] ,'sed - lowland' ,'sed - swamps' , data_plot{[ID2ndData] , 1}} ,'FontSize',15,'location', 'southeast')
    else
        legend( { 'sed - lowland' ,'sed - swamps' , data_plot{[ID2ndData] , 1}} ,'FontSize',15,'location', 'southeast')
    end
end

end