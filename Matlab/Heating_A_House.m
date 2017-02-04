

%% Heating for a Three-Room House - Newton's Law of Cooling.
%USING DSOLVE METHOD:

clear, clc, close all
% Heating/cooling constants.
k1 = 0.15; k2 = 0.15; k3 = 0.15; % Insulation is very good for the main room.
k0 = 0.5; k4 =1; % Poor insulation in the attic and basement. Attic is drafty. 


% Temperatures and the Heater
TE = 40; % Temperature of the earth, is constant and cold.
TS = 50 ; % The outside air is a cool 50 degrees Fahrenheit.
H = 10; % The heater can heat the main room at 10 degrees per hour if there are no losses. 

% The differential equations for the three rooms.
syms t TB(t) TM(t) TA(t)
% Temperature TB in the basement
EQ1 = diff(TB, t) == -k0*(TB-TE) -k1*(TB-TM);
% Temperature TM on the main floor
EQ2 = diff(TM,t) == -k1*(TM-TB)-k2*(TM-TS)-k3*(TM-TA)+H;
% Temperature TA in the attic
EQ3 = diff(TA,t) == -k3*(TA-TM)-k4*(TA-TS);

% Tell MATLAB it's OK to use variable-precision arithmetic.
% Use variableprecision arithmetic ; use at least 4 sig figs
% This saves Matlab a lot of effort trying to find exact solns @ 32 digits.
EQNS = vpa([EQ1, EQ2, EQ3], 4) ;


% Assume these initial conditions. TB(0)=50, TM(0)=50, TA(0)=50 & Solve.
sol = dsolve( EQNS, [TB(0)==50, TM(0)==50, TA(0)==50] );

% Let's make these solutions matlab functions:
TB = matlabFunction(sol.TB);
TM = matlabFunction(sol.TM);
TA = matlabFunction(sol.TA);


% Now to Plot/Visualize:
% Plot the temperature in all three rooms:
figure
hold on; grid on
t_pts = 0: 0.01: 24; % Model temperatures in the home for one day.
% Plot the temperature in the basement.
plot(t_pts, TB(t_pts), 'LineWidth', 3 )
plot(t_pts, TM(t_pts), 'LineWidth', 3 )
plot(t_pts, TA(t_pts), 'LineWidth', 3 )
set(gca, 'FontSize', 20)
xlabel('Time in hours')
ylabel('Temperature in ^{\circ}F')
title('Heating a Three-Room House')
legend('Basement','Main Floor', 'Attic', 'Location', 'Best')

%% Solving the same system with Laplace Transforms
close all, clc

% Heating/cooling constants.
k1 = 0.15; k2 = 0.15; k3 = 0.15; % Insulation is very good for the main room.
k0 = 0.5; k4 =1; % Poor insulation in the attic and basement. Attic is drafty.
% Temperatures and the Heater
TE = 40; % Temperature of the earth, is constant and cold.
TS = 50 ; % The outside air is a cool 50 degrees Fahrenheit.
H = 10; % The heater can heat the main room at 10 degrees per hour if there are no losses.
T0 = [50 50 50]'; % All three rooms start off at 50 degrees F.

% Define the matrix:
A = [-(k0+k1), k1, 0; k1, -(k1+k2+k3), k3; 0, k3, -(k3+k4)];
% Define the forcing terms:
F = [k0*TE; k2*TS+H; k4*TS];

% Taking the Laplace transform given:
%
% L{T} = (sI-A)^-1 * [T(0)+(1/s)*F]
%
syms s
I = eye(3);
L_Solns = inv(s*I - A)*(T0+(1/s)*F);

% Now find the solutions in the T domain:
T_Solns = vpa(ilaplace( L_Solns ), 4); % The solution to the time domain.
T_On = matlabFunction(T_Solns); % Make this into a function to plot easily.

% Let's find the final temperature in each room as T -> infinity using the
% final value theorem:
Tinf = vpa(limit(s*L_Solns,0), 4)


% Now to Plot/Visualize:
% Plot the temperature in all three rooms:
figure
hold on; grid on
t_pts = 0: 0.01: 24; % Model temperatures in the home for one day.

% Plot the temperature in the basement.
plot(t_pts, T_On(t_pts), 'LineWidth', 3 )
set(gca, 'FontSize', 20)
xlabel('Time in hours')
ylabel('Temperature in ^{\circ}F')
title('Heating a Three-Room House')
legend('Basement','Main Floor', 'Attic', 'Location', 'Best')
plot(24, Tinf(3), 'redo', 'MarkerFaceColor', 'Yellow') %Final Condition for Attic
plot(24, Tinf(2), 'blueo', 'MarkerFaceColor', 'Red') %Final Condition for Main Floor
plot(24, Tinf(1), 'blueo', 'MarkerFaceColor', 'Blue') %Final Condition for Basement

clearvars -except T_On T_Off
%% Let's prove the power of the Laplace transform domain on a discontinuous system"
close all, clc

% Redefine:
% Heating/cooling constants.
k1 = 0.15; k2 = 0.15; k3 = 0.15; % Insulation is very good for the main room.
k0 = 0.5 ; k4 =1; % Poor insulation in the attic and basement. Attic is drafty.
% Temperatures and the Heater
TE = 40; % Temperature of the earth, is constant and cold.
TS = 50 ; % The outside air is a cool 50 degrees Fahrenheit.
T0 = [50 50 50]'; % All three rooms start off at 50 degrees F.


syms t
H = @(t) 10*(1-heaviside(t-10)); % Turn off the heater at time t = 10.
% Forcing terms again:
F = [k0*TE k2*TS+H(t) k4*TS ]';
Fs = laplace(F);


% Define the matrix:
A = [-(k0+k1), k1, 0; k1, -(k1+k2+k3), k3; 0, k3, -(k3+k4)];


% Taking the Laplace transform given:
%
% L{T} = (sI-A)^-1 * [T(0)+(1/s)*F]
%
syms s
I = eye(3);
L_Solns = inv(s*I - A)*(T0+(1/s)*F);

% Now find the solutions in the T domain:
T_Solns = vpa(ilaplace( L_Solns ), 4); % The solution to the time domain.
T_Off = matlabFunction(T_Solns); % Make this into a function to plot easily.


% Now to Plot/Visualize:
% Plot the temperature in all three rooms:
figure
hold on; grid on
t_pts = 0: 0.01: 24; % Model temperatures in the home for one day.

% Plot the temperature in the basement.
plot(t_pts, T_Off(t_pts), 'LineWidth', 3 )
set(gca, 'FontSize', 20)
xlabel('Time in hours')
ylabel('Temperature in ^{\circ}F')
title('Heating a Three-Room House')
legend('Basement','Main Floor', 'Attic', 'Location', 'Best')

% Notice how easily the Laplacian transform could handle the discontinuity
% in the forcing term.

clearvars -except T_On T_Off
%% A Final Observation on the system: A lag when heater turns on / off
% This SHOULD be giving sharktooth, why is it not?
Range = 0: .8: 24;  % Model temperatures in the home for one day.

A = [];
MF = [];
B = [];
j = 0;

for i = Range; % Model temperatures in the home for one day.

    if j ~= 0
        if MF(j) >= 65
            All_Temp = T_Off(i);
            B = [B, All_Temp(1)];
            MF = [MF,All_Temp(2)];
            A = [A, All_Temp(3)];
        else
            All_Temp = T_On(i);
            B = [B, All_Temp(1)];
            MF = [MF,All_Temp(2)];
            A = [A, All_Temp(3)];
        end
    
    else
    All_Temp = T_On(i);
    B = [B, All_Temp(1)];
    MF = [MF,All_Temp(2)];
    A = [A, All_Temp(3)];
    end
    
    j = j + 1;
end

figure
hold on; grid on
plot(Range, A, 'LineWidth', 3 )
plot(Range, MF, 'LineWidth', 3 )
plot(Range,B, 'LineWidth', 3 )
set(gca, 'FontSize', 20)
xlabel('Time in hours')
ylabel('Temperature in ^{\circ}F')
title('Heating a Three-Room House')
legend('Basement','Main Floor', 'Attic', 'Location', 'Best')


