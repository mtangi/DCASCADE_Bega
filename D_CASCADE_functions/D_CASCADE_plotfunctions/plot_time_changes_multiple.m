function plot_time_changes_multiple(data_plot_multiple , reaches , varargin)
%PLOT_TIME_CHANGES_PROVENANCE plots different subplots of the network showing the
%changes in two data in the data_plot matrix with time, for different
%simulation scenarios defined in data_plot_multiple.

%% default settings 

def_linewidth = 2.5;
name_cMap1 = 'summer';
name_cMap2 = 'copper';
def_timeinterval = [];

%% read additional inputs 

p = inputParser;
addOptional(p,'timeinterval',def_timeinterval);
addOptional(p,'LineWidth',def_linewidth);
addOptional(p,'FontSize',12);
addOptional(p,'titles',[]);
addOptional(p,'reachID',[]);
addOptional(p,'plotID',[ 1 2 ]);
addOptional(p,'equal_axis', 0);  %equal axis guarantee all subplots have the same scale in the y axis
addOptional(p,'background', []);
addOptional(p,'Q_ARR_reach', []);
addOptional(p,'showlegend', 'on');
addOptional(p,'cMap_bcg', [] );


parse(p,varargin{:})

linewidth = p.Results.LineWidth;
timeinterval = p.Results.timeinterval;
Fsize = p.Results.FontSize;
titles_subplot = p.Results.titles;
reach_id = p.Results.reachID;
plot_id = p.Results.plotID;
equal_axis = p.Results.equal_axis;
bcg_id = p.Results.background;
Q_ARR_reach = p.Results.Q_ARR_reach;
showlegend = p.Results.showlegend;
cMap_bcg = p.Results.cMap_bcg;

if isempty(timeinterval )
    timeinterval = [2 min(cellfun('size',data_plot_multiple{1,1}(plot_id,2),1))];
end
    
%% find constant colorbar for plot

cMapName1 = [name_cMap1 '(' num2str(size(data_plot_multiple,1)) ')'];
cMap1 = eval(cMapName1);

cMapName2 = [name_cMap2 '(' num2str(size(data_plot_multiple,1)+1) ')'];
cMap2 = eval(cMapName2);  

%% find subplot design
if length(reaches) > 12
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

if ~isempty(bcg_id)
    bcg_window = 3;
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
    yinterval1 = [ min(cell2mat(cellfun(@(x)min(x{plot_id(1),2}(:, reaches),[],'all'),data_plot_multiple(:,1),'UniformOutput',0)))  ; max(cell2mat(cellfun(@(x)max(x{plot_id(1),2}(:, reaches),[],'all'),data_plot_multiple(:,1),'UniformOutput',0))) ] ; 
    yinterval2 = [ min(cell2mat(cellfun(@(x)min(x{plot_id(2),2}(:, reaches),[],'all'),data_plot_multiple(:,1),'UniformOutput',0)))  ; max(cell2mat(cellfun(@(x)max(x{plot_id(2),2}(:, reaches),[],'all'),data_plot_multiple(:,1),'UniformOutput',0))) ] ; 
end

% plot data in frame
fig = figure;
set(gcf,'color','w');

set (fig, 'Units', 'normalized', 'Position', [0.1,0.1,0.8,0.8]);
for i=1:length(reaches)
    subplot(n_subplot(1),n_subplot(2),i);
      
        % if i want the background coloration
    if ~isempty(bcg_id) && equal_axis
        lm = max(yinterval1(2),yinterval2(2)); %define the heigth of the bars, give by the limits of the plot
        M = movmean(data_plot_multiple{1,1}{bcg_id,2}(:,i),bcg_window);  %I use if i want the moving average
        datacolor = round((M-min(M))/(max(M)-min(M)) *100 ) + 1; % attribute to each value in the timestep a class form 1 to 100 given their value
        
        if ~isempty(Q_ARR_reach) %if the background is Q, use the color of the prcentiles in Q_ARR_reach
            M = matlab.tall.movingWindow(@(x) max(x),bcg_window,data_plot_multiple{2,1}{bcg_id,2}(:,i)); % this function find for each element the maximum in a window defined by bcg_window
            datacolor = sum(M >= Q_ARR_reach(i,:),2); 
        end 
        
        if max(M) == min(M) %if there is no variation in the backgroud plot variable
            mapcolor = repmat(cMap_bcg(1,:),size(datacolor)); % i pur all colors as the ligth color
        else
            mapcolor = cMap_bcg(datacolor,:); % Attribute to each class a color in the colorscale cMap_bcg previously defined 
        end
        b = bar(repmat(lm ,[1 max(timeinterval)]),1,'EdgeColor', 'None');
        b.FaceColor = 'flat';
        b.CData = mapcolor;
        hold on
    end
    
    %plot n1
    for p = 1:size(data_plot_multiple,1)
        plot(data_plot_multiple{p,1}{plot_id(1),2}(:,reaches(i) ) , 'Color', cMap1(p,:) , 'LineWidth' , linewidth)
        hold on
    end
    
    if equal_axis; ylim( yinterval1 ); end

    yyaxis right

    %plot n2
    for p = 1:length(data_plot_multiple)    
        plot(data_plot_multiple{p,1}{plot_id(2),2}(:,reaches(i)) , 'Color', cMap2(p,:) , 'LineWidth' , linewidth, 'LineStyle','-')
        hold on
    end
    
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
        title( [titles_subplot{i}]);
    end
    
    hold off
    
    xlim(timeinterval)

    grid minor
    set(gca,'XGrid','on','FontSize',Fsize)

end

%% set legend

if strcmp(showlegend,'on')
    if ~isempty(bcg_id)
        legend([['bckground : ' data_plot_multiple{1,1}{bcg_id,1}]; data_plot_multiple(:,2); data_plot_multiple(:,2)] ,'FontSize',15,'location', 'southeast')
    else
        legend([data_plot_multiple(:,2); data_plot_multiple(:,2)] ,'FontSize',15,'location', 'southeast')
    end
end


%% put annotation with legend

% str = {[name_cMap1 ' : ' data_plot_multiple{1, 1}{plot_id(1),1}  ]  [name_cMap2 ' : ' data_plot_multiple{1, 1}{plot_id(2),1}] };
% a = annotation('textbox',[0.70 0 0.3 0.3],'String',str,'FitBoxToText','on',...
%  'BackgroundColor', 'w','FontSize',12);
% 
end
