% Define initial constants for simulation
% We will do Apple Inc. (AAPL) stock option pricing

% Define initial constants for Apple Inc. option pricing simulation
numTimeSteps = 8000;        % Number of time steps in the grid
numStockSteps = 2000;       % Number of stock price steps in the grid
minStockPrice = 0.4;        % Minimum stock price to consider (Smin must be greater than 0)
maxStockPrice = 250;        % Maximum stock price to consider, above current market price

% Option parameters based on current market conditions for Apple Inc.
maturity = 0.5;             % Time to maturity in years (6 months)
strikePrice = 180;          % Strike price of the option
volatility = 0.30;          % Annual volatility (30%)
riskFreeRate = 0.035;       % Annual risk-free interest rate (3.5%)
dividendYield = 0.6;      % Annual dividend yield (0.6%)
isCallOption = true;        % Specifies the option type (true for call, false for put)


% Generate option price surface using the explicit method
[timeValues, stockPrices, optionSurface] = blackScholesExplicitCovariance(numTimeSteps, numStockSteps, minStockPrice, maxStockPrice, maturity, strikePrice, volatility, riskFreeRate, dividendYield, isCallOption);

% Define check points for interpolation
stockCheckPoints = linspace(0, 1.5 * strikePrice, 20);
timeCheckPoints = linspace(0, maturity, 20);

% Initialize matrix for interpolated values
checkSurface = zeros(length(stockCheckPoints), length(timeCheckPoints));
for i = 1:length(timeCheckPoints)
    for j = 1:length(stockCheckPoints)
        checkSurface(i, j) = interp2(stockPrices, timeValues, optionSurface, stockCheckPoints(j), timeCheckPoints(i));
    end
end

% Plotting the interpolated surface
figure; % Ensures a new figure window is created
surfHandle = surf(stockCheckPoints, fliplr(timeCheckPoints), checkSurface);
xlabel('Stock Price');
ylabel('Time Until Maturity');
title('Interpolated Option Price Surface');
colormap(jet); % Applies the 'jet' colormap which is vibrant and visually appealing

% Interpolate a single point and display its value
selectedStockPrice = 30;
selectedTime = 0.25;
interpolatedPoint = interp2(stockPrices, timeValues, optionSurface, selectedStockPrice, selectedTime);
fprintf('At stock price = %f and time = %f, the option price is %f\n', selectedStockPrice, selectedTime, interpolatedPoint);



