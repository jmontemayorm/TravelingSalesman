function bestRoute = TSP_F(maxIter,initialTemp,finalTemp,linearTemp,cities,timeout,partialFileName)
    %% Traveling Salesman Problem as Function

    % Input will be a list of cities and their respective coordinates
    % The shortest route that goes through each city only once and returns to
    % the origin city needs to be calculated

    % Requirements:
    % - The number of cities will vary between 5 and 50
    % - The code must accept two vectors, the combination of those vectors
    %   results in the nth city coordinate
    % - The number of iterations needs to be flexible to be defined just before
    %   a new set of cities
    % - The output must include the route and the calculated distance (save)

    %% Settings
    alpha = (finalTemp / initialTemp)^(1 / maxIter); % Exponential factor
    eta = (initialTemp - finalTemp) / maxIter; % Linear factor

    % Print and save settings
    printElapsedTime = 1; % Set to 1 to print execution time
    saveMetadata = 1; % Settings and additional information about the run
    saveResults = 1;

    %% Main
    tic

    N = length(cities);

    % Initial arbitrary order (shortest distances)
    route = zeros(1,N + 1); % Allocate
    route(1) = 1; % Start city
    remainingCities = true(1,N); % Remaining cities
    remainingCities(1) = false; % Remove first city

    % Add the closest city going from the start point
    for n = 2:(N - 1)
        % Get all neighbor distances
        neighborDistances = sqrt(sum((repmat(cities(route(n-1),:),N-n+1,1) - cities(remainingCities,:)).^2,2));

        % Find shortest distance
        shortestIdx = find(neighborDistances == min(neighborDistances),1);

        % Find the corresponding index (last of mappedIdx) and remove city
        mappedIdx = find(remainingCities,shortestIdx);
        remainingCities(mappedIdx(end)) = false;

        % Update route
        route(n) = mappedIdx(end);
    end

    % Add the remaining city and go back to the origin
    route(N) = find(remainingCities,1);
    route(N + 1) = route(1);

    % Initial distance
    idxStart = 1:N; % Precalculated indices for the route (starting point)
    idxEnd = idxStart + 1; % Precalculated indices for the route (ending point)
    distance = sum(sqrt(sum((cities(route(idxEnd),:) - cities(route(idxStart),:)).^2,2)));

    % Results storage
    bestRoute = route;
    bestDistance = distance;

    % Optimiaztion loop
    iter = 0;
    T = initialTemp;
    elapsedTime = toc;
    while iter < maxIter && elapsedTime < timeout
        % Get the indices to swap
        swappers = randi(N,[1,2]);

        % Swap the cities
        newRoute = route;
        newRoute(swappers(1)) = route(swappers(2));
        newRoute(swappers(2)) = route(swappers(1));
        newRoute(end) = newRoute(1);

        % Calculate the distance of the new route
        newDistance = sum(sqrt(sum((cities(newRoute(idxEnd),:) - cities(newRoute(idxStart),:)).^2,2)));

        % Accept or reject
        if newDistance <= distance
            route = newRoute;
            bestRoute = newRoute;
            distance = newDistance;
            bestDistance = newDistance;
        elseif exp((distance - newDistance)/T) >= rand % Extra probability of acceptance
            route = newRoute;
            distance = newDistance;
        end

        % Update temperature
        if linearTemp == 1
            T = T - eta;
        else
            T = T * alpha;
        end

        % Check for timeOut
        if mod(iter,10000)
            elapsedTime = toc;
        end

        % Increase the counter!
        iter = iter + 1;
    end

    elapsedTime = toc;

    %% Results
    % Print elapsed time
    if printElapsedTime == 1
        fprintf('Elapsed time is %0.6f seconds.\n',elapsedTime)
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

    %% Save reults into a file
    if saveResults == 1
        % Output folder
        if exist('getOutputFolder','file') == 2
            outF = getOutputFolder(mfilename('fullpath'));
        else
            outF = pwd;
        end

        % Create a new file
        fileName = fullfile(outF,sprintf('%s.txt',partialFileName));

        fId = fopen(fileName,'w','n','UTF-8');
        if fId >= 3 % Success opening the file
            % Metadata
            if saveMetadata == 1
                % Iterations
                fprintf(fId,'Number of iterations: %i\n',iter);

                % Temperature
                fprintf(fId,'Initial system ''temperature'': %0.2f\n',initialTemp);
                if linearTemp == 1
                    fprintf(fId,'Used a linear temperature model:\n');
                    fprintf(fId,'\tT(t) = T(0) - eta * t\n');
                    fprintf(fId,'\teta = %0.6f (aprox)\n',eta);
                else
                    fprintf(fId,'Used an exponential temperature model:\n');
                    fprintf(fId,'\tT(t) = T(0) * alpha ^ t\n');
                    fprintf(fId,'\talpha = %0.6f (aprox)\n',alpha);
                end
                fprintf(fId,'\tT(%i) = %0.6f\n',iter,T);

                % Save elapsed time
                fprintf(fId,'Elapsed time was %0.6f seconds\n',elapsedTime);

                % Time 
                dt = datetime('now');
                dt.Format = 'uuuu-MM-dd HH:mm';
                fprintf(fId,'Run on: %s\n\n',char(dt));
            end

            % Write reults
            fprintf(fId,'The shortest distance found was %0.6f\n',bestDistance);
            fprintf(fId,'The best route was the following:\n');
            for n = 1:N
                fprintf(fId,'\t%02i.-\t%02i => %02i\n',n,bestRoute(n),bestRoute(n+1));
            end

            % Close the file
            fclose(fId);
        end

        % Save figure and other data
        saveas(f,fullfile(outF,sprintf('%s.eps',partialFileName)),'epsc')
        save(fullfile(outF,partialFileName),'cities','bestRoute','bestDistance')
    end
end