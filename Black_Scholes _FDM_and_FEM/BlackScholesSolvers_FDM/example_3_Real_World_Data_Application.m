% Use with real-world data
maxStockPrice = 10000;      % Maximum stock price to consider, above current market price
riskFreeRate = 0.0186;      % Annual risk-free interest rate (1.86%)
dividendYield = 0.0291;     % Annual dividend yield (2.91%)
volatility = 0.1325;        % Annual volatility (13.25%)
strikePrice = 140;          % Strike price of the option
maturity = 135 / 365;       % About 135 days to expiration

[timeValues, stockPrices, optionSurface] = blackScholesExplicitCovariance(numTimeSteps, numStockSteps, minStockPrice, maxStockPrice, maturity, strikePrice, volatility, riskFreeRate, dividendYield, isCallOption);

stockInspect = 144.16;
value = interp2(stockPrices, timeValues, optionSurface, stockInspect, 0);
fprintf('Calculated option price: %f\n', value);

%% Plot entire surface (for debugging)
figure; % Creates a new figure window
surfHandle = surf(stockPrices, fliplr(timeValues), optionSurface, 'EdgeColor', 'none'); % EdgeColor set to 'none' for smoother surface
xlabel('Stock Price', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Time Until Maturity', 'FontSize', 12, 'FontWeight', 'bold');
zlabel('Option Price', 'FontSize', 12, 'FontWeight', 'bold');
title('Generated Surface for Real-World Data', 'FontSize', 14, 'FontWeight', 'bold');

colorbar; % Adds a color bar to indicate the scale of option prices
shading interp; % Uses interpolated shading to smooth colors
colormap(jet); % Applies the 'jet' colormap which is vibrant and visually appealing
axis tight; % Fits the axes box tightly around the data
grid on; % Enable grid to improve readability of the plot's spatial context
set(gca, 'FontSize', 10, 'FontWeight', 'bold'); % Set the font size and weight for axes

% Optional: Highlight the maximum option price point
flippedTimeValues = fliplr(timeValues);  % Properly flip timeValues for indexing
[maxOptionPrice, maxIndex] = max(optionSurface(:));
[maxRow, maxCol] = ind2sub(size(optionSurface), maxIndex);
hold on;
plot3(stockPrices(maxCol), flippedTimeValues(maxRow), maxOptionPrice, '+', 'MarkerSize', 15, 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black');
text(stockPrices(maxCol), flippedTimeValues(maxRow), maxOptionPrice, sprintf(' Max Price: %.2f', maxOptionPrice), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'black');
hold off;

