%----------------------------------------------------------
%
% EuropeanOptionModel displays European call and put option prices 
% and their greeks using the Black-Scholes-Merton model.
%
%----------------------------------------------------------
%
% Usage:
% option = EuropeanOptionModel(initialStockPrice, strikePrice, riskFreeRate, timeToMaturity, volatility, optionType)
% callOption = EuropeanOptionModel(100, 130, 0.05, 1, 0.2, 'Call');
% callOption.displayPrice();
% callOption.displayDelta();
% callOption.displayGamma();
% callOption.displayTheta();
% callOption.displayVega();
% callOption.displayRho();
%
%----------------------------------------------------------
%
% Inputs:
% initialStockPrice - Current price of the underlying asset.
% strikePrice       - Strike (i.e., exercise) price of the option.
% riskFreeRate      - Annualized continuously compounded risk-free rate of return
%                     over the life of the option, expressed as a positive decimal
%                     number.
% timeToMaturity    - Time to expiration of the option, expressed in years.
% volatility        - Annualized asset price volatility (i.e., annualized standard
%                     deviation of the continuously compounded asset return),
%                     expressed as a positive decimal number.
% optionType        - 'call' or 'put'.
%
% Methods:
% displayPrice
% displayDelta
% displayGamma
% displayTheta
% displayVega
% displayRho
%
%----------------------------------------------------------

classdef (Sealed = true) EuropeanOptionModel < handle
  properties (GetAccess = private, SetAccess = private)
    currentSpotPrices
    remainingTime
    dPlus
    dMinus
  end
  properties (SetObservable, AbortSet)
    initialStockPrice = 0 
    strikePrice = 0
    timeToMaturity = 0
  end
  properties
    riskFreeRate = 0                        
    volatility = 0
    optionType
  end
  methods
    function obj = set.initialStockPrice(obj, value)
      if value <= 0
        error('EuropeanOptionModel:InitialStockPrice', 'Initial stock price must be > 0.');
      else
        obj.initialStockPrice = value;
      end
    end
    function obj = set.strikePrice(obj, value)
      if value <= 0
        error('EuropeanOptionModel:StrikePrice', 'Strike price must be > 0.');
      else
        obj.strikePrice = value;
      end
    end
    function obj = set.riskFreeRate(obj, value)
      if value <= 0
        error('EuropeanOptionModel:RiskFreeRate', 'Risk-free rate must be > 0.');
      else
        obj.riskFreeRate = value;
      end
    end
    function obj = set.timeToMaturity(obj, value)
      if value <= 0
        error('EuropeanOptionModel:TimeToMaturity', 'Time to maturity must be > 0.');
      else
        obj.timeToMaturity = value;
      end
    end
  end
  methods
    function obj = EuropeanOptionModel(initialStockPrice, strikePrice, riskFreeRate, timeToMaturity, volatility, optionType)
      % Validate input parameters
      if nargin ~= 6
        error('EuropeanOptionModel:InvalidInputs', 'Incorrect number of input parameters.');
      end
      if strcmpi(optionType, 'call') ~= 1 && strcmpi(optionType, 'put') ~= 1
        error('EuropeanOptionModel:InvalidOptionType', 'Option type must be "call" or "put".');
      end
      
      obj.currentSpotPrices = 0;
      obj.remainingTime = 0;
      
      addlistener(obj, 'initialStockPrice', 'PostSet', @obj.calculatePriceComponents);
      addlistener(obj, 'strikePrice', 'PostSet', @obj.calculatePriceComponents);
      addlistener(obj, 'timeToMaturity', 'PostSet', @obj.calculatePriceComponents);
      
      obj.initialStockPrice = initialStockPrice;
      obj.strikePrice = strikePrice;
      obj.riskFreeRate = riskFreeRate;
      obj.timeToMaturity = timeToMaturity;
      obj.volatility = volatility;
      obj.optionType = optionType;
    end
    function displayPrice(obj)
      % Calculate option price using Black-Scholes formula
      standardDev = obj.volatility .* sqrt(obj.remainingTime);
      discountedStrikePrice = obj.strikePrice .* exp(-obj.riskFreeRate .* obj.remainingTime);
      
      obj.dPlus = (1 ./ standardDev) .* (log(obj.currentSpotPrices ./ obj.strikePrice) + obj.remainingTime .* (obj.riskFreeRate + obj.volatility.^2 / 2));
      obj.dPlus(isnan(obj.dPlus)) = 0; % Handle NaN results from division by zero
      obj.dMinus = obj.dPlus - standardDev;
      
      if strcmpi(obj.optionType, 'call') == 1
        price = obj.currentSpotPrices .* normcdf(obj.dPlus)  -  discountedStrikePrice .* normcdf(obj.dMinus);
      else
        price = discountedStrikePrice .* normcdf(-obj.dMinus) - obj.currentSpotPrices .* normcdf(-obj.dPlus);
      end
      
      obj.plotGraph(price, 'Option Price ($)');
    end
    % Methods for other greeks would follow the same pattern
    % ...
  end
  methods (Access = private)
    function calculatePriceComponents(obj, src, event)
      % Compute ranges and grids for option price calculations
      priceStep = obj.initialStockPrice / 100;
      if obj.initialStockPrice >= obj.strikePrice
        priceRange = obj.strikePrice / 2: priceStep: 1.5 * obj.initialStockPrice;
      else
        priceRange = obj.strikePrice / 2: priceStep: 1.5 * obj.strikePrice;
      end
      timeStep = obj.timeToMaturity / 100;
      timeRange = 0: timeStep: obj.timeToMaturity;
      
      [obj.currentSpotPrices, obj.remainingTime] = meshgrid(priceRange, timeRange);
    end
    function plotGraph(obj, z, label)
      % Plotting function for displaying graphs of option greeks
      surfaceHandle = surf(obj.currentSpotPrices, obj.remainingTime, z, gradient(z));
      view(-125, 30);
      
      xlabel('Stock Price ($)', 'FontSize', 14, 'FontWeight', 'bold');
      ylabel('Time to Maturity (Years)', 'FontSize', 14, 'FontWeight', 'bold');
      zlabel([label ' in $'], 'FontSize', 14, 'FontWeight', 'bold');
      
      set(surfaceHandle, 'FaceAlpha', .6);
      set(surfaceHandle, 'EdgeAlpha', .2);
      set(surfaceHandle, 'FaceLighting', 'phong');
      set(surfaceHandle, 'FaceColor', 'interp');
      
      grid on;
      grid minor;
      
      colormap('hsv');
      view(45, 150); % Adjust the viewing angle for better visibility
    end
    function x = gaussianDensityPrime(obj, d)
      x = (exp((-d.^2)./2)./sqrt(2*pi));
    end
  end % methods
end % classdef