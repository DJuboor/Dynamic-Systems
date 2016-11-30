%% Beginning of Questions 8-10


%% Custom Direction Field for second-order differential equation:

% ay'' + by' + cy = 0
clear, clc, close all
% BLOCK 1: Convert 2nd-order DE to a system of first order DEs.
% The new system will have the form dx/dt = F(t, x) where x is the state
% vector. So, x = [x1; x2]
% Define the constants here.
% Define F(t, x) as an anonymous function.
% Below, x is a column vector with two components
F = @(t, x) [x(2); -x(1)-((1-((x(1))^2)*x(2))/2)] ;

% BLOCK 2: Set the plot window dimensions here.
% Customize the plot settings.
x1min = -5; x1max= 5; x2min = -5; x2max =5;
figure % Pop up a new figure
axis([x1min x1max x2min x2max ])
axis square
grid on; hold on
set(gca, 'FontSize', 20)
title('Van der Pol Oscillator')
xlabel('x1 = y'); ylabel('x2 = dy/dt')

radius = 0.1 % Control the size of the tick marks
spacing_horizontal = 0.5; % Control their horizontal spacing.
spacing_vertical = 0.5; % Control their vertical spacing.
my_color = [0.25, 0.25, 0.25]; % Control their color.
optional_dots = 1; % Control whether a dot is placed. Use 0 or 1
time = -5:5;

for x1p =x1min: spacing_horizontal: x1max
 for x2p = x2min: spacing_vertical: x2max
 x=[x1p;x2p]; % Update the state vector x.
 dx = F(time,x); % Find dx/dt
 % Add a line here to replace dx by a unit vector.
 mag=(sqrt((dx(1))^2+(dx(2)^2)));
 dx(1)=dx(1)/mag;
 dx(2)=dx(2)/mag;
 dx1 = dx(1) *radius; dx2 = dx(2)*radius; % Scale tick marks.
 % Plot the tick marks!
 plot([x1p - dx1 , x1p + dx1 ], [x2p - dx2 , x2p + dx2] , 'Color', my_color);
 end
end




%plots one solution
tStart = [0]; tEnd = [30];
both_directions = 1; % Change to 0 to draw in the forward direction only.
show_initial_point = 1; % Change to 0 to suppress drawing the initial points.
x0 = [1;0];
[~, x_out] = ode45(F, [tStart, tEnd], x0);
comet(x_out(:,1), x_out(:,2))
% Show initial point using a green dot.
if show_initial_point
 plot(x0(1), x0(2), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'yellow')
end

for n = -5:.1:5
x0 = [n; 0]
% Block of code here
[~, x_out] = ode45(F, [tStart, tEnd], x0);
comet(x_out(:,1), x_out(:,2))
% Show initial point using a green dot.
if show_initial_point
 plot(x0(1), x0(2), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'yellow')
end
end

