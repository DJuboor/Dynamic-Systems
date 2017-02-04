clc, clear
x0 = [2; 0]; 
syms s x1(t) x2(t) 
x = [x1; x2];
Dx = diff(x,t);
A= [0 1;-26 -2];
DE = Dx ==A*x;
Soln = dsolve(DE, x(0)==x0);
 
x1 = matlabFunction(Soln.x1);  % make x1 a function
 
x2 = matlabFunction(Soln.x2);  % make x2 a function
 
% View the answers in nice form. 
x1(t)
x2(t)

%% 

A = [ 0 1 ; -26 -2];
b = [0;20];
Xeq = inv(A)*b;

%% 
clear, clc

syms F(t) F(s)
F(t) = 20*(heaviside(t-5)-heaviside(t-10));
F(t) = matlabFunction(F(t));

F(s) = laplace(F(t));

%% 

x0 = [2; 0];
syms s X_hom(s)
A = [ 0 1 ; -26 -2];
sI = s* eye(2)

X_hom(s) = inv(sI-A)*x0

%% 
syms F(t) F(s) X_forced(s)
F(t) = 20*(heaviside(t-5)-heaviside(t-10));
F(t) = matlabFunction(F(t));
A = [ 0 1 ; -26 -2];
sI = s* eye(2);
B = [0;20];
F(s) = laplace(F(t));

X_forced(s) = inv(sI-A)*B*F(s);

%% 
clearvars -except X_forced X_hom
clc

syms s t
X(s) = X_hom(s) + X_forced(s);

X = matlabFunction(ilaplace(X(s)));

figure
hold on, grid on
plot(2, 0, 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'y')
plot(0, 0, 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'g')
plot(15.2847, -0.0712, 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'r')


t = 0:.01:5;
    temp = X(t);
    plot(temp(1,:),temp(2,:),'g', 'LineWidth',2)

t = 5:.01:10;
    temp = X(t);
    plot(temp(1,:),temp(2,:),'b', 'LineWidth',2)


t = 10:.01:15;
    temp = X(t);
    plot(temp(1,:),temp(2,:),'r', 'LineWidth',2)


