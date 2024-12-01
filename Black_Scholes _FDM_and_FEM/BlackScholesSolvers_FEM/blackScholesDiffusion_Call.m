clc, clear, close all
% Given Data => TSLA
strike_price = 900; % Strike Price
current_stock_price = 1000; % Current Stock Price
implied_volatility = 80/100; % Implied Volatility
interest_rate = 0.025; % Assume 2 Percent. Close to US 10 Year Treasury Yield

time_to_expiration_years = 255/365.25; % Years Until Expiration
time_step_days = 25.5/365.25; % 1 Day Time Steps

maximum_stock_price = 2*current_stock_price; % Maximum Stock Price
num_blocks = 25; % Number of Elements
price_step = maximum_stock_price/num_blocks; % Price Step

num_nodes = num_blocks + 1; % Linear Approximation
fixed_nodes = [1, num_nodes]; % Node numbers of fixed DOFs
free_nodes = setxor(1:num_nodes,fixed_nodes); % Node numbers of free DOFs

%Mesh and Time Discretization:
stock_price_vector = (0:price_step:maximum_stock_price)'; % Stock Price Vector from $0 to Maximum Stock Price
time_vector = (0:time_step_days:time_to_expiration_years)'; % Time Vector, where tao = T - t
num_time_nodes = length(time_vector); % Number of Time Nodes

%Boundary and Initial Conditions
initial_option_price_at_zero = 0; % Option Price at Stock Price = $0 - Boundary Condition ($)
option_price_at_max_stock_price = maximum_stock_price - strike_price*exp(-interest_rate*time_vector); % Option Price at Maximum Stock Price - Boundary Condition ($)
initial_option_prices = max(stock_price_vector - strike_price, 0); % Free Global Option Price Vector at tao = 0 - Initial Condition ($)

total_option_price_matrix = zeros(num_nodes, num_time_nodes); % Initialize Option Price Storage Matrix for All Time Steps
total_option_price_matrix(:, 1) = initial_option_prices; % Initial Condition at tao = 0
total_option_price_matrix(end, :) = option_price_at_max_stock_price; % Boundary Condition at Manually Set Maximum Stock Price

syms stock_price
dplus_T = (log(stock_price/strike_price)+(interest_rate+implied_volatility^2/2)*(time_to_expiration_years))/(implied_volatility*(time_to_expiration_years)^(1/2));
dminus_T = (log(stock_price/strike_price)+(interest_rate-implied_volatility^2/2)*(time_to_expiration_years))/(implied_volatility*(time_to_expiration_years)^(1/2));
Exact_Solution = stock_price*normcdf(dplus_T) - strike_price*exp(-interest_rate*time_to_expiration_years)*normcdf(dminus_T);

%% Define Connectivity of Global Free DOFs
connectivity_matrix = zeros(num_blocks, 2); % Initialize element to node connectivity
for i = 1:num_blocks
    connectivity_matrix(i, 1) = i; % Local node 1
    connectivity_matrix(i, 2) = i+1; % Local node 2
end

%% Define Shape Functions
syms n 
N1n = 1/2*(1-n); N2n = 1/2*(1+n); % Linear shape functions in local coordinates
dndS = 2/price_step; % Relationship between dS and dn in the form of dn/dS
dN1dn = diff(N1n); dN2dn = diff(N2n); % Derivative of shape functions wrt n
dN1dS = dN1dn*dndS; dN2dS = dN2dn*dndS; % Derivative of shape functions wrt S

%% Derive Parts of 'Stiffness' Matrix using Method of Weighted Residuals
k11eB1 = int((1-n)*dN1dS*N1n, [-1 1]); k12eB1 = int((1-n)*dN1dS*N2n, [-1 1]);
k21eB1 = int((1-n)*dN2dS*N1n, [-1 1]); k22eB1 = int((1-n)*dN2dS*N2n, [-1 1]);
kB1 = [k11eB1, k12eB1; k21eB1, k22eB1];

k11eB2 = int((1+n)*dN1dS*N1n, [-1 1]); k12eB2 = int((1+n)*dN1dS*N2n, [-1 1]);
k21eB2 = int((1+n)*dN2dS*N1n, [-1 1]); k22eB2 = int((1+n)*dN2dS*N2n, [-1 1]);
kB2 = [k11eB2, k12eB2; k21eB2, k22eB2];

k11eD = int(N1n^2, [-1 1]); k12eD = int(N1n*N2n, [-1 1]);
k21eD = int(N2n*N1n, [-1 1]); k22eD = int(N2n^2, [-1 1]);
kD = [k11eD, k12eD; k21eD, k22eD];

k11eE1 = int((1-n)^2*dN1dS^2, [-1 1]); k12eE1 = int((1-n)^2*dN1dS*dN2dS, [-1 1]);
k21eE1 = int((1-n)^2*dN2dS*dN1dS, [-1 1]); k22eE1 = int((1-n)^2*dN2dS^2, [-1 1]);
kE1 = [k11eE1, k12eE1; k21eE1, k22eE1];

k11eE2 = int((1-n)*(1+n)*dN1dS^2, [-1 1]); k12eE2 = int((1-n)*(1+n)*dN1dS*dN2dS, [-1 1]);
k21eE2 = int((1-n)*(1+n)*dN2dS*dN1dS, [-1 1]); k22eE2 = int((1-n)*(1+n)*dN2dS^2, [-1 1]);
kE2 = [k11eE2, k12eE2; k21eE2, k22eE2];

k11eE3 = int((1+n)^2*dN1dS^2, [-1 1]); k12eE3 = int((1+n)^2*dN1dS*dN2dS, [-1 1]);
k21eE3 = int((1+n)^2*dN2dS*dN1dS, [-1 1]); k22eE3 = int((1+n)^2*dN2dS^2, [-1 1]);
kE3 = [k11eE3, k12eE3; k21eE3, k22eE3];

alpha_coefficient = implied_volatility^2/price_step*[2, -1; -1, 2];

%% Global 'Stiffness' Matrix Assembly (Apply Theta Differencing)
A = zeros(num_nodes, num_nodes); % Initialize Global 'Stiffness' Matrix for the next time step
B = zeros(num_nodes, num_nodes); % Initialize Global 'Stiffness' Matrix for the current time step
theta = 1/2; % Central difference method
for i = 1:num_blocks
    S1 = stock_price_vector(i); S2 = stock_price_vector(i+1); % Local node stock price coordinates
    enodes = connectivity_matrix(i, :); % Element to node connectivity for local to global mapping
    beta = kD\((implied_volatility^2-interest_rate)*(1/2*S1*kB1 + 1/2*S2*kB2) + 1/8*...
        implied_volatility^2*(S1^2*kE1 + 2*S1*S2*kE2 + S2^2*kE3) + interest_rate*kD); % [beta] in the hand calculation
    ket1 = alpha_coefficient\(1/time_step_days - beta*(1-theta)); % Local 'Stiffness' Matrix for the current time step
    ket2 = alpha_coefficient\(1/time_step_days + beta*theta); % Local 'Stiffness' Matrix for the next time step
    A(enodes, enodes) = A(enodes, enodes) + ket2; % Assembly of the global 'stiffness' matrix
    B(enodes, enodes) = B(enodes, enodes) + ket1; % Assembly of the global forcing vector
end

%% Global 'Stiffness' Matrix Partitioning
A_fixed = A(fixed_nodes, fixed_nodes);
A_fixed_free = A(fixed_nodes, free_nodes);
A_free_fixed = A(free_nodes, fixed_nodes);
A_free = A(free_nodes, free_nodes);

%% Solving for Option Pricing for all Stock Price Increments and Time Steps
for t = 2:num_time_nodes
    V_fixed_next = total_option_price_matrix(fixed_nodes, t);
    RS = B*initial_option_prices;
    V_free_next = A_free\(RS(2:end-1) - A_free_fixed*V_fixed_next); % Free DOF Option pricing at next time step
    total_option_price_matrix(2:end-1, t) = V_free_next;
    initial_option_prices = total_option_price_matrix(:, t);
end
total_option_price_matrix

%% Post Processing
% Find the FEM Approximated Option Price at Current Stock Price and tao = T
current_index = find(stock_price_vector == max(stock_price_vector(stock_price_vector <= current_stock_price)));
% Linear Interpolation to Find Option Price In-Between Nodes
option_price_at_current_stock_price = (total_option_price_matrix(current_index+1,end)-total_option_price_matrix(current_index,end))...
    /price_step*(current_stock_price-stock_price_vector(current_index)) + total_option_price_matrix(current_index,end);

%% Plotting
figure(1)
time_in_days = time_vector*365.25; % Change time scale to days for graphing
surf(time_in_days, stock_price_vector, total_option_price_matrix)
colormap jet  % Modern and vibrant color map
lighting phong  % Enhances lighting effects
camlight left  % Adds a light to enhance 3D effect
xlabel('Time Until Expiration (Days)', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Stock Price ($)', 'FontSize', 12, 'FontWeight', 'bold')
zlabel('Option Price ($)', 'FontSize', 12, 'FontWeight', 'bold')
title('Black-Scholes Model - Call Option Pricing', 'FontSize', 14, 'FontWeight', 'bold')
colorbar  % Adds a color bar to indicate the scale of option prices
axis tight  % Removes excess space around the data
grid on  % Adds a grid for better readability



figure(2)
fplot(stock_price, Exact_Solution, [0, maximum_stock_price], 'LineWidth', 2, 'Color', 'blue')
hold on
plot(stock_price_vector, total_option_price_matrix(:, end), 'r--', 'LineWidth', 2, 'Marker', 'o', 'MarkerIndices', 1:10:length(stock_price_vector), 'MarkerFaceColor', 'red')
hold off
xlabel('Stock Price ($)', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Option Price ($)', 'FontSize', 12, 'FontWeight', 'bold')
title('Black-Scholes Equation - Call Option Pricing at tao = T', 'FontSize', 14, 'FontWeight', 'bold')
legend('Black-Scholes Equation', 'FE Approximation', 'Location', 'best')
grid on  % This line adds grid lines to the plot
axis tight  % This line ensures the axes are snugly fit around the data
set(gca, 'FontSize', 10)  % Sets a consistent font size for axis ticks
ylim([0, max(max(total_option_price_matrix(:, end)), double(subs(Exact_Solution, stock_price, maximum_stock_price)))*1.1]);  % This ensures all data is visible and the plot does not cut off significant figures


%% Print Results
fprintf('The current option price using Finite Element Method is $%.4f\n', option_price_at_current_stock_price)
fprintf('The current option price using Black-Scholes Equation is $%.4f\n', subs(Exact_Solution, current_stock_price))
