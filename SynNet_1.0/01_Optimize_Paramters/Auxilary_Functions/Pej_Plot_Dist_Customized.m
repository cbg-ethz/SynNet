% This function makes a violin distribution plots of the columns in "data"
% place at position "x", with "color", in the "direction", and full or
% "Dashed".
% Example: 
% Pej_Plot_Dist([4 4], repmat(X,1,2),  [1 0 0],[-1; 1]);
% this plots X in a symmetric violin at point 4, in red.
% NOTE: if you want to plot some log-normal data, just set the figure YScal
% to log and it goes automatically:
% figure; 
% set(gca, 'YScale', 'log');
% Pej_Plot_Dist([4 4], repmat(X,1,2),  [1 0 0],[-1; 1]);


% Written by pejman.m@gmail.com
% 2013 Sept.
% -----------------

function Pej_Plot_Dist_Customized(x,data, Color, Direction, Dashed, Scale)
if nargin<2; data = x; x = 1:size(data,2); end

MidPoints =  geomean(data);

if strcmpi(get(gca, 'YScale'), 'log')
    %% estimate densities from the log(data), and plot it back on a semilogY scale.
    data = log(data);
    isLog = true;
else
    isLog = false;
end

N = size(data,2);
if nargin<3 || isempty(Color);       Color      = repmat([.6,0,.05], N, 1); end
if nargin<4 || isempty(Direction);   Direction  = ones(N,1); end
if nargin<5 || isempty(Dashed);      Dashed     = false; end
if nargin<6 || isempty(Scale);       Scale      = nan; end

if size(Color,1)    == 1; Color     = repmat(Color    ,N,1); end
if size(Direction,1)== 1; Direction = repmat(Direction,N,1); end
if length(Scale)    == 1; Scale     = repmat(Scale    ,N,1); end

for i = 1:N
    [f,xi] = ksdensity(data(:,i),linspace(min(data(:,i)), max(data(:,i)), 300));
    if isnan(Scale(i))
    f = f/max(f) * .4;
    else
    f = f * Scale(i);
    end
    if Dashed
        
        tmpF = [f; f];
        tmpF(1,1:2:end)=0; tmpF(2,2:2:end)=0;
        tmpXi= [xi;xi];
        
        f = tmpF(:)';
        xi = tmpXi(:)';
    end
    
%     X975 = quantile(data(:,i), .975);
%     X025 = quantile(data(:,i), .025);
%     f(xi>X975)  = []; f(xi<X025)  = [];
%     xi(xi>X975) = []; xi(xi<X025) = [];
    if isLog
        xi = exp(xi);
    end
    hold on
    fill(x(i) + Direction(i) * [f 0 0],[xi xi(end) xi(1)], Color(i,:), 'edgecolor', Color(i,:));
%     plot(x(i)+ [0 max(f)* Direction(i)], [1 1 ] * MidPoints(i), '-', 'color', 1*[1 1 1], 'linewidth', 1.2  , 'markersize', 5 )
    %plot(x(i), MidPoints(i), 's', 'color', Color(i,:), 'linewidth', .5   , 'markersize', 3 , 'markerfacecolor', Color(i,:))
    plot(x(i), MidPoints(i), 'o', 'color', [1 1 1], 'linewidth', .5   , 'markersize', 6, 'markerfacecolor', 0*[1 1 1])
    plot(x(i), MidPoints(i), '+', 'color', [1 1 1], 'linewidth', 1  , 'markersize', 4 )
end
end 
