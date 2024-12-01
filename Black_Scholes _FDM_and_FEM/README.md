# Black-Scholes Option Pricing Implementation

## Overview
This repository contains implementations of the Black-Scholes option pricing model using various numerical methods including Finite Element Method (FEM) and Finite Difference Method (FDM). The implementation focuses on both European and American options, with capabilities for analyzing call and put options.

![Figure1](./Poster.jpg)

## Features
- Multiple implementation methods:
  - Analytical Black-Scholes solution
  - Finite Element Method (FEM)
  - Explicit Finite Difference Method
  - Implicit Finite Difference Method
- Support for both Call and Put options
- Real-world data application examples
- Comparative analysis between numerical and analytical solutions
- Interactive visualization tools
- Greeks calculation (Delta, Gamma, Theta, Vega, Rho)

## Mathematical Foundation
The implementation is based on the Black-Scholes partial differential equation:

```
∂V/∂t + 1/2 σ²S² ∂²V/∂S² + rS ∂V/∂S - rV = 0
```

Where:
- V: Option value
- S: Stock price
- t: Time
- σ: Volatility
- r: Risk-free interest rate

## Usage Examples

### Basic Option Pricing
```matlab
% Example for European Call Option
initialStockPrice = 100;
strikePrice = 100;
riskFreeRate = 0.1;
timeToMaturity = 1;
volatility = 0.6;

callOption = EuropeanOptionModel(initialStockPrice, strikePrice, ...
    riskFreeRate, timeToMaturity, volatility, 'Call');
callOption.displayPrice();
```

### Real-World Application
```matlab
% Tesla (TSLA) option example
strike_price = 900;
current_stock_price = 1000;
implied_volatility = 0.80;
interest_rate = 0.025;
time_to_expiration_years = 255/365.25;
```

## Model Assumptions and Limitations
1. Early Exercise Limitation
2. Constant Parameters
3. Continuous Trading Assumption
4. No Transaction Costs
5. Ignoring Bid/Ask Spreads
6. Lognormal Price Distribution
7. Dividend Omission
8. Ideal Market Conditions

## Numerical Methods Implementation

### Finite Element Method (FEM)
- Uses linear shape functions
- Implements Galerkin weighted residual method
- Handles boundary conditions for both call and put options

### Finite Difference Methods
- Explicit scheme for direct computation
- Implicit scheme for enhanced stability
- Crank-Nicolson method for improved accuracy

## Visualization Capabilities
The implementation includes various visualization tools:
1. 3D surface plots of option prices
2. Comparative plots of analytical vs numerical solutions
3. Greeks visualization
4. Time evolution of option prices

## Requirements
- MATLAB R2019b or later
- Financial Toolbox (optional, for comparative analysis)
- Symbolic Math Toolbox

## References
1. Black, F. and Scholes, M. (1973). The Pricing of Options and Corporate Liabilities. Journal of Political Economy, 81(3), pp.637-654.
2. Capinski, M. and Kopp, E. (2012). The Black-Scholes Model. Cambridge University Press.
3. Hull, J. Fundamentals of Futures and Options Markets, Chapter 13.