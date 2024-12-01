function [timeValues, stockPrices, optionSurface] = blackScholesExplicitCovariance(numTimeSteps, numStockSteps, minStockPrice, maxStockPrice, maturity, strikePrice, volatility, riskFreeRate, dividendYield, isCallOption)

    % Create empty option price surface
    % This surface will be filled in using the explicit method
    optionSurface = zeros(1 + numTimeSteps, 1 + numStockSteps);
    
    % Change variables to transformed time and transformed stock price
    maxTau = volatility^2 * maturity / 2;
    deltaTau = maxTau / numTimeSteps;
    tauValues = 0:deltaTau:maxTau;
    
    % Define transformed stock price grid
    % Note: The transformation is based on the Black-Scholes PDE
    maxX = log(maxStockPrice / strikePrice);
    minX = log(minStockPrice / strikePrice);
    deltaX = (maxX - minX) / numStockSteps;
    xValues = minX:deltaX:maxX;
    
    % We define the transformation constants for the Black-Scholes PDE
    % These are used to simplify the PDE and make it easier to solve
    transformConst1 = riskFreeRate / (volatility^2 / 2);
    transformConst2 = (riskFreeRate - dividendYield) / (volatility^2 / 2);
    alpha = -0.5 * (transformConst2 - 1);
    beta = -0.25 * (transformConst2 - 1)^2 - transformConst1;
    
    % Define the boundary condition function for the transformed PDE
    boundaryConditionFunction = @(x, tau) exp((0.25 * (transformConst2 - 1)^2 + transformConst1) * tau) * ...
        max(exp(0.5 * (transformConst2 + 1) * x) - exp(0.5 * (transformConst2 - 1) * x), 0);
    
     % Reverse transformation to get back to original variables
    timeValues = maturity - 2 * tauValues / volatility^2; %Should we use dot division?
    stockPrices = strikePrice * exp(xValues);
    
    % Calculate stability factor and display it
    stabilityFactor = deltaTau / (deltaX^2);
    fprintf('Stability factor is %f (must be <=0.5)\n', stabilityFactor);
    
    % We set the initial and boundary conditions for the transformed PDE
    if isCallOption
        optionSurface(:,1) = 0;
        optionSurface(:,end) = boundaryConditionFunction(maxX, tauValues);
        optionSurface(1,:) = boundaryConditionFunction(xValues, 0);
    else
        optionSurface(:,1) = 0;
        optionSurface(:,end) = strikePrice;
        optionSurface(1,:) = payoff(stockPrices, strikePrice, isCallOption);
    end
    
    % Fill in the rest of the surface by iterating through time and stock price
    % We use the explicit method to solve the transformed PDE
    for timeIndex = 2:numTimeSteps + 1
        for stockIndex = 2:numStockSteps
            optionSurface(timeIndex, stockIndex) = stabilityFactor * optionSurface(timeIndex - 1, stockIndex - 1) + ...
                (1 - 2 * stabilityFactor) * optionSurface(timeIndex - 1, stockIndex) + ...
                stabilityFactor * optionSurface(timeIndex - 1, stockIndex + 1);
        end
        
        % APply boundary condition at stock price = 0
        optionSurface(timeIndex, 2:numStockSteps) = max(optionSurface(timeIndex, 2:numStockSteps), boundaryConditionFunction(xValues(2:numStockSteps), tauValues(timeIndex)));
    end
    
    % Transform the option surface back to the original variables
    [meshX, meshTau] = meshgrid(xValues, tauValues);
    optionSurface = strikePrice * exp(alpha * meshX + beta * meshTau) .* optionSurface;
    
    end
    