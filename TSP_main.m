%% Traveling Salesman Problem - Parallel Search

% This code searches for an optimized route in a list of cities by
% simultaneously running (in a parallel pool) 6 settings variations of the
% Metropolis algorithm.

% A maximum number of iterations and maximum timeout time can be set to
% stop the search, both work inside each call of the function so in order
% for the timeout to work properly the parallel pool requires at least 6
% workers (or the number of calls, if changed).

% Both here and in the function TSP_F the whole route distance calculation
% is optimized by indexing and being a single line statement (intermediate
% variables result in longer times).

%% Settings
% Data to load
year = 2019; % 2018 or 2019
dateNum = 2; % Date number, from 0 to 5

% General settings
maxIter = 1e8 * ones(1,6); % 100 M, might be overridden by timeout
minsToTimeout = 2; % Number of minutes

% Temperature settings (one per parallel process)
initialTemp = [1 10 1 1 10 100];
finalTemp = [1e-12 1e-12 1e-12 1e-12 1e-12 1e-12];
linearTemp = [1 1 1 0 0 0]; % 1 for linear, 0 for exponential

% Calculated timeout setting
timeout = minsToTimeout*60 - 10; % Substract 10 seconds for overhead

%% Load data
dataName = sprintf('%i_fecha%i',year,dateNum);
load(dataName)
cities = vec_all(:,2:3);

%% Calculation of the best route
% Start overall timing
tic

% Allocate for results
bestRoutes = cell(1,6);

% Simultaneously call the 6 variations of the settings
parfor c = 1:6
    bestRoutes{c} = TSP_F(maxIter(c),initialTemp(c),finalTemp(c),linearTemp(c),cities,timeout,sprintf('TravelingSalesman_%s%i',dataName,c));
end

% Evaluate distances
N = length(cities);
idxStart = 1:N; % Precalculated indices for the route (starting point)
idxEnd = idxStart + 1; % Precalculated indices for the route (ending point)
allDistances = zeros(1,6); % Allocate space
for c = 1:6
    allDistances(c) = sum(sqrt(sum((cities(bestRoutes{c}(idxEnd),:) - cities(bestRoutes{c}(idxStart),:)).^2,2)));
end

% Choose best distance
bestIdx = find(allDistances == min(allDistances),1);

%% Results
% Print results
bestRoute = bestRoutes{bestIdx};
fprintf('The best route was the following:\n')
for n = 1:N
    fprintf('\t%02i.-\t%02i => %02i\n',n,bestRoute(n),bestRoute(n+1))
end

% Plot results
f = figure(1);
hold off
scatter(cities(:,1),cities(:,2))
hold on
plot([cities(bestRoute(idxStart),1); cities(bestRoute(end),1)],[cities(bestRoute(idxStart),2); cities(bestRoute(end),2)])
xRange = max(cities(:,1)) - min(cities(:,1));
xMin = min(cities(:,1)) - 0.1 * xRange;
xMax = max(cities(:,1)) + 0.1 * xRange;
yRange = max(cities(:,2)) - min(cities(:,2));
yMin = min(cities(:,2)) - 0.1 * yRange;
yMax = max(cities(:,2)) + 0.1 * yRange;
axis([xMin xMax yMin yMax])

% Print overall time
elapsedTime = toc;
fprintf('Total elapsed time is %0.6f seconds.\n',elapsedTime);

% Print shortest distance
fprintf('The shortest distance found was %0.6f\n',allDistances(bestIdx))
