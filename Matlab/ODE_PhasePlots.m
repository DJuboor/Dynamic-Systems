%%  k == 0
Y = @(t)[exp(-3*t) + 5.*t.*exp(-3.*t)];
Dy = @(t)[2.*exp(-3.*t) - 15.*t.*exp(-3.*t)];

% LINE 1 %
figure
subplot(2,2,1)
t = 0: 0.01 : 10;
plot(t, Y(t), 'b', 'LineWidth',2)
title('Solution Curve y(t)')
xlabel('Time in Seconds')
ylabel('Y(t)')
grid on
hold on 

subplot(2,2,3)
plot(t, Dy(t), 'r', 'LineWidth',2)
title('Derivative of y')
xlabel('Time in Seconds')
ylabel('Dy/Dt')
hold on
grid on

subplot(1,2,2)
plot(Y(t), Dy(t), 'k', 'LineWidth', 2)
title('Phase Plot')
xlabel('Y(t)')
ylabel('Dy/Dt')
hold on
grid on


set(gca, 'FontSize', 20)
subplot(2,2,[2,4]), plot(1, 2, 'redo', 'MarkerFaceColor', 'red') %IC
subplot(2,2,[2,4]), plot(0, 0, 'cyano', 'MarkerFaceColor', 'cyan') %Equilibrium point
legend('phase', 'IC', 'Equilibrium', 'Location', 'Best')

% LINE 2 % k == 10
Y = @(t)[exp(-3.*t) + 5.*t.*exp(-3.*t) + 5.*t.^2.*exp(-3.*t)];
Dy = @(t)[2.*exp(-3.*t) - 5.*t.*exp(-3.*t) - 15.*t.^2.*exp(-3.*t)];

subplot(2,2,1)
t = 0: 0.01 : 10;
plot(t, Y(t), 'b-.', 'LineWidth',2)


subplot(2,2,3)
plot(t, Dy(t), 'r-.', 'LineWidth',2)


subplot(2,2,[2,4])
plot(Y(t), Dy(t),'k-.', 'LineWidth', 2)

% LINE 3 %  k == 20
Y = @(t)[exp(-3.*t) + 5.*t.*exp(-3.*t) + 10.*t.^2.*exp(-3.*t)];
Dy = @(t)[2.*exp(-3.*t) + 5.*t.*exp(-3.*t) - 30.*t.^2.*exp(-3.*t)];

k = 20;
subplot(2,2,1)
t = 0: 0.01 : 10;
plot(t, Y(t),'b:', 'LineWidth',2)


subplot(2,2,3)
plot(t, Dy(t), 'r:', 'LineWidth',2)


subplot(1,2,2)
plot(Y(t), Dy(t),'k:', 'LineWidth', 2)

% LINE 4 % k == 30
Y = @(t)[exp(-3.*t) + 5.*t.*exp(-3.*t) + 15.*t.^2.*exp(-3.*t)];
Dy = @(t)[2.*exp(-3.*t) + 15.*t.*exp(-3.*t) - 45.*t.^2.*exp(-3.*t)];

k = 30;
subplot(2,2,1)
t = 0: 0.01 : 10;
plot(t, Y(t),'b--', 'LineWidth',2)
legend('k=0','k=10','k=20','k=30', 'Location', 'Best')

subplot(2,2,3)
plot(t, Dy(t), 'r--', 'LineWidth',2)
legend('k=0','k=10','k=20','k=30', 'Location', 'Best')

subplot(1,2,2)
plot(Y(t), Dy(t),'k--', 'LineWidth', 2)