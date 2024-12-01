function [timeValues, stockPrices, optionSurface] = blackScholesNaiveImplicit(numTimeSteps, numStockSteps, minStockPrice, maxStockPrice, maturity, strikePrice, volatility, riskFreeRate, dividendYield, isCallOption)

    % Create empty option price surface
    optionSurface = zeros(1 + numTimeSteps, 1 + numStockSteps);
    
    % Determine time and stock price increments
    deltaTime = maturity / numTimeSteps;
    deltaStock = (maxStockPrice - minStockPrice) / numStockSteps;
    
    % Create grid for time and stock prices
    timeValues = 0:deltaTime:maturity;
    stockPrices = minStockPrice:deltaStock:maxStockPrice;
    
    % Set boundary conditions
    if isCallOption
        optionSurface(:,1) = 0;
        optionSurface(:,end) = maxStockPrice - strikePrice;
        optionSurface(end,:) = payoff(stockPrices, strikePrice, isCallOption);
    else
        optionSurface(:,1) = strikePrice;
        optionSurface(:,end) = 0;
        optionSurface(end,:) = payoff(stockPrices, strikePrice, isCallOption);
    end
    
    % Coefficients for implicit matrix setup
    coeffA = @(index) 0.5 * (riskFreeRate - dividendYield) * index * deltaTime - 0.5 * volatility^2 * index^2 * deltaTime;
    coeffB = @(index) 1 + volatility^2 * index^2 * deltaTime + riskFreeRate * deltaTime;
    coeffC = @(index) -0.5 * (riskFreeRate - dividendYield) * index * deltaTime - 0.5 * volatility^2 * index^2 * deltaTime;
    
    % Fill in the option values using backward induction
    for timeIndex = numTimeSteps:-1:1   
        matrixA = diag(coeffA(2:numStockSteps-1), -1) + diag(coeffB(2:numStockSteps)) + diag(coeffC(1:numStockSteps-2), 1);
        valuesVector = optionSurface(timeIndex + 1, 2:numStockSteps)';
        valuesVector(1) = valuesVector(1) - coeffA(1) * optionSurface(timeIndex, 1);
        valuesVector(end) = valuesVector(end) - coeffC(numStockSteps + 1) * optionSurface(timeIndex, numStockSteps + 1);
    
        optionSurface(timeIndex, 2:numStockSteps) = matrixA \ valuesVector;
    
        % Enforce free boundary condition
        optionSurface(timeIndex, 2:numStockSteps) = max(optionSurface(timeIndex, 2:numStockSteps), payoff(stockPrices(2:numStockSteps), strikePrice, isCallOption));
    end
    
    end