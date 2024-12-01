% Define initial constants for computation

% Grid configuration for simulation
numTimeSteps = 4000;    % Number of time steps in the grid
numStockSteps = 1000;   % Number of stock price steps in the grid

% Stock price range
minStockPrice = 0.4;    % Minimum stock price to consider => Must be greater than 0!
maxStockPrice = 1000;   % Maximum stock price to consider, above current market price

% Option characteristics
strikePrice = 10;       % Strike price of the option
maturity = 1;           % Time to maturity in years
isCallOption = true;    % Option type (true for call, false for put)

% Market and volatility parameters
volatility = 0.4;       % Annual volatility
riskFreeRate = 0.02;    % Annual risk-free interest rate
dividendYield = 0;      % Annual dividend yield


% Generate the option price surface for comparison
[timeValues, stockPrices, optionSurface] = blackScholesExplicitCovariance(numTimeSteps, numStockSteps, minStockPrice, maxStockPrice, maturity, strikePrice, volatility, riskFreeRate, dividendYield, isCallOption);

% Compare American and European option prices
stockCheckPoints = linspace(0, 15, 20);
timeCheckPoints = linspace(0, maturity, 20);

% Initialize matrices for comparison
europeanPrices = zeros(length(stockCheckPoints), length(timeCheckPoints));
americanPrices = zeros(length(stockCheckPoints), length(timeCheckPoints));

for i = 1:length(timeCheckPoints)
    for j = 1:length(stockCheckPoints)
        [europeanPrices(i, j), ~] = blsprice(stockCheckPoints(j), strikePrice, riskFreeRate, maturity - timeCheckPoints(i), volatility);
        americanPrices(i, j) = interp2(stockPrices, timeValues, optionSurface, stockCheckPoints(j), timeCheckPoints(i));
    end
end

% Calculate and plot differences
priceDifferences = europeanPrices - americanPrices;
figure; % Creates a new figure window
surfHandle = surf(stockCheckPoints, fliplr(timeCheckPoints), abs(priceDifferences), 'EdgeColor', 'none'); % EdgeColor set to 'none' for smoother surface
xlabel('Stock Price', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Time Until Maturity', 'FontSize', 12, 'FontWeight', 'bold');
zlabel('Price Difference', 'FontSize', 12, 'FontWeight', 'bold');
title('Price Differences Between European and American Options', 'FontSize', 14, 'FontWeight', 'bold');

colorbar;
shading faceted; 
colormap(jet); 
axis tight; % Fits the axes box tightly around the data
grid on; % Enable grid to improve readability of the plot's spatial context
set(gca, 'FontSize', 10, 'FontWeight', 'bold'); % Set the font size and weight for axes

% Find and highlight the maximum difference
[maxDiff, maxIndex] = max(abs(priceDifferences(:)));
[maxRow, maxCol] = ind2sub(size(priceDifferences), maxIndex);
hold on;
plot3(stockCheckPoints(maxCol), timeCheckPoints(length(timeCheckPoints) - maxRow + 1), maxDiff, '+', 'MarkerSize', 15, 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'black');
text(stockCheckPoints(maxCol), timeCheckPoints(length(timeCheckPoints) - maxRow + 1), maxDiff, sprintf(' Max Diff: %.1e', maxDiff), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10, 'Color', 'black');
hold off;


