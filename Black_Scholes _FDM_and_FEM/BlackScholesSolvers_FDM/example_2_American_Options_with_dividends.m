% American options with non-zero dividend yield
dividendYield = 0.08;  % Non-zero dividend yield
[timeValues, stockPrices, optionSurface] = blackScholesExplicitCovariance(numTimeSteps, numStockSteps, minStockPrice, maxStockPrice, maturity, strikePrice, volatility, riskFreeRate, dividendYield, isCallOption);

% Check results with non-zero dividend for American calls
americanPricesBinomial = zeros(length(stockCheckPoints), length(timeCheckPoints));
americanPricesFDM = zeros(length(stockCheckPoints), length(timeCheckPoints));
for i = 1:length(timeCheckPoints)
    for j = 1:length(stockCheckPoints)
        
        % Note: for simplicity, the usage of the generated binomial trees is done inefficiently
        % binprice fails at T=0 and also for large values of Smax
        if maturity - timeCheckPoints(i) ~= 0
            [~, option] = binprice(stockCheckPoints(j), strikePrice, riskFreeRate, maturity - timeCheckPoints(i), 0.001, volatility, isCallOption, dividendYield);
            americanPricesBinomial(i,j) = option(1,1);
        else
            americanPricesBinomial(i,j) = payoff(stockCheckPoints(j), strikePrice, isCallOption);
        end
        americanPricesFDM(i,j) = interp2(stockPrices, timeValues, optionSurface, stockCheckPoints(j), timeCheckPoints(i));
    end
end
differences = americanPricesBinomial - americanPricesFDM;

% Plot differences
figure;
surfHandle = surf(stockCheckPoints, fliplr(timeCheckPoints), abs(differences), 'EdgeColor', 'none'); % EdgeColor set to 'none' for smoother surface
xlabel('Stock Price', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Time Until Maturity', 'FontSize', 12, 'FontWeight', 'bold');
title('Differences With binprice (no free boundary check)', 'FontSize', 14, 'FontWeight', 'bold');

colorbar;
shading faceted;
colormap(jet);
axis tight;
grid on; 
set(gca, 'FontSize', 10, 'FontWeight', 'bold');


