k = 0;

for i = -10:10

    
my_color = [(k * 1/21),0,1-(k*1/21)];

time_pts = -10 : 0.1 : 10;
y_pts = logistic(time_pts, i, 1, 1);
plot(time_pts, y_pts, 'Color', my_color, 'LineWidth', 3) % The standard curve.
grid on
set(gca, 'FontSize', 20)
title('Logistic Curves')
xlabel('time'); ylabel('y(t)')
hold on

k = k+1;
end

%% 

clear, clc
syms P(t)
r=2; K = 10;
DE = diff(P,t) == r*P * (1 - P/K)
sol = simplify( dsolve(DE, P(0)==1) )
p = matlabFunction(sol) % Use little p for the solution. 

%% 
figure

sol = []
for i = 0:.5:10
    sol =[sol, p(i)];
end

    plot(0:.5:10, sol, 'blue', 'LineWidth', 3) % The standard curve.
    grid on
    set(gca, 'FontSize', 20)
    title('Logistic Equation with r=2, K=10')
    xlabel('time'); ylabel('P(t)')
    axis([0,10,0,12])
    hold on
    
    plot(0,p(0),'Marker','o','MarkerSize',12,'Color','b','MarkerFaceColor','y')
    


%%

%Euler's Method Example and exact solution using dsolve
clear, clc, close all

% First declare the function f(t,y) which gives the slope dy/dt.
r = 2; K = 10;
f = @(t,P) r*P * (1 - P/K)

% Initialize the Variables
dt = 1; % Step size. Also called h in the notes.
tStart = 0; yStart = 1; % Initial point.
tEnd = 10; % Stopping time.

% Define time points and solution vector
t_points = tStart: dt: tEnd;
y_points = zeros(size(t_points)); % Use zeros as place holders for now.

% Initialize the solution at the initial condition.
y_points(1) = yStart;

% Implement Euler's method using a for loop.
for i=2:length(t_points)
 yprime = f(t_points(i-1),y_points(i-1));
 y_points(i) = y_points(i-1) + dt*yprime;
end

% Plot Solutions
figure(1)
plot(t_points, y_points, 'red*:', 'LineWidth', 3)
grid on
xlabel('Time')
ylabel('P(t)')
set(gca, 'FontSize', 20)
title('Logistic Equation with r=2, K=10')
axis([0 10 0 12])
hold on

% Now let's find the exact solution.
syms P(t)
r=2; K = 10;
DE = diff(P,t) == r*P * (1 - P/K)
sol = simplify( dsolve(DE, P(0)==1) )
p = matlabFunction(sol) % Use little p for the solution. 

t_points = tStart: 0.01: tEnd;
plot(t_points, p(t_points), 'blue', 'LineWidth', 3)
plot(tStart, yStart, 'blueo', 'MarkerSize', 16, 'MarkerFaceColor', 'yellow')





% NEXT EULERS

% Initialize the Variables
dt = 1/4; % Step size. Also called h in the notes.
tStart = 0; yStart = 1; % Initial point.
tEnd = 10; % Stopping time.

% Define time points and solution vector
t_points = tStart: dt: tEnd;
y_points = zeros(size(t_points)); % Use zeros as place holders for now.

% Initialize the solution at the initial condition.
y_points(1) = yStart;

% Implement Euler's method using a for loop.
for i=2:length(t_points)
 yprime = f(t_points(i-1),y_points(i-1));
 y_points(i) = y_points(i-1) + dt*yprime;
end

% Plot Solutions
figure(1)
plot(t_points, y_points, 'black*:', 'LineWidth', 3)
grid on
xlabel('Time')
ylabel('P(t)')
set(gca, 'FontSize', 20)
title('Logistic Equation with r=2, K=10')
axis([0 10 0 12])
hold on


% PLOT ODE45
[t_out, y_out] = ode45(f, [tStart tEnd], yStart);

% Plot Solutions


figure(1)
plot(t_out, y_out, 'redd','MarkerFaceColor','green', 'LineWidth', 1)
grid on
xlabel('Time')
ylabel('P(t)')
set(gca, 'FontSize', 20)
title('Logistic Equation with r=2, K=10')
axis([0 10 0 12])
hold on

legend('Euler''s Method with h = 1', 'Exact','Initial Point','Euler''s method with h = 1/4','ode45', 'Location', 'Best')