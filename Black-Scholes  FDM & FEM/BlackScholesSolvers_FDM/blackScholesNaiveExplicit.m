function [timeValues, stockPrices, optionSurface] = blackScholesNaiveExplicit(numTimeSteps, numStockSteps, minStockPrice, maxStockPrice, maturity, strikePrice, volatility, riskFreeRate, dividendYield, isCallOption)

    % Create empty option price surface
    optionSurface = zeros(1 + numTimeSteps, 1 + numStockSteps);
    
    % Compute time and stock price increments
    deltaTime = maturity / numTimeSteps;
    deltaStock = (maxStockPrice - minStockPrice) / numStockSteps;
    
    % Create grid for time and stock prices
    timeValues = 0:deltaTime:maturity;
    stockPrices = minStockPrice:deltaStock:maxStockPrice;
    
    % Set initial boundary conditions
    if isCallOption
        optionSurface(:,1) = 0;
        optionSurface(:,end) = maxStockPrice - strikePrice;
        optionSurface(end,:) = payoff(stockPrices, strikePrice, isCallOption);
    else
        optionSurface(:,1) = strikePrice;
        optionSurface(:,end) = 0;
        optionSurface(end,:) = payoff(stockPrices, strikePrice, isCallOption);
    end
    
    % Calculate option values using finite differences
    priceIncrementCoeff = @(index) 1 / (1 + riskFreeRate * deltaTime) * ...
        (-0.5 * (riskFreeRate - dividendYield) * index * deltaTime + 0.5 * volatility^2 * index^2 * deltaTime);
    valueDecrementCoeff = @(index) 1 / (1 + riskFreeRate * deltaTime) * ...
        (1 - volatility^2 * index^2 * deltaTime);
    valueIncrementCoeff = @(index) 1 / (1 + riskFreeRate * deltaTime) * ...
        (0.5 * (riskFreeRate - dividendYield) * index * deltaTime + 0.5 * volatility^2 * index^2 * deltaTime);
    
    % Iterate through time and stock prices to calculate option values
    % We use the explicit method to solve the Black-Scholes PDE
    for timeIndex = numTimeSteps:-1:1
        coefficientMatrix = diag(priceIncrementCoeff(2:numStockSteps-1), -1) + ...
                            diag(valueDecrementCoeff(2:numStockSteps)) + ...
                            diag(valueIncrementCoeff(1:numStockSteps-2), 1);
        for stockIndex = 2:numStockSteps
            optionSurface(timeIndex, stockIndex) = priceIncrementCoeff(stockIndex - 1) * optionSurface(timeIndex + 1, stockIndex - 1) + ...
                                                   valueDecrementCoeff(stockIndex) * optionSurface(timeIndex + 1, stockIndex) + ...
                                                   valueIncrementCoeff(stockIndex) * optionSurface(timeIndex + 1, stockIndex + 1);
        end
    
        % We apply the boundary condition at stock price = 0
        optionSurface(timeIndex, 2:numStockSteps) = max(optionSurface(timeIndex, 2:numStockSteps), payoff(stockPrices(2:numStockSteps), strikePrice, isCallOption));
    end
    
    end
    