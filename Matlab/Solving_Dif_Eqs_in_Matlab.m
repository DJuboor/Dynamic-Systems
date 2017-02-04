
B = 2;
K = 100;

syms x

f = (2*x^5)/(1+x^2);
int(f,x)

%% 

syms y(t)
DE = diff(y,t) == y - t; 
%dsolve(DE)

Dy = diff(y);
dsolve(DE,y(0)==1)

%% 

B = 2;
K = 100;

syms V(t);
DE = diff(V,t) == B*V*log(K/V);

dsolve(DE, V(0)==1)

%%

syms V(t) B K

eqn1 = diff(V,t) == B*V*log(K/V)
eqn2 = diff(eqn1,t);


