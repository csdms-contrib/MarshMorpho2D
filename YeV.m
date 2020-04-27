function [Hs,Tp]=YeV(fetch,wind,h)
g=9.8;
delta=h*g./wind.^2;
chi=fetch*g./wind.^2;
epsilon=3.64*10^-3*(tanh(0.493*delta.^0.75).*tanh(3.13*10^-3*chi.^0.57./tanh(0.493*delta.^0.75))).^1.74;
ni=0.133*(tanh(0.331*delta.^1.01).*tanh(5.215*10^-4*chi.^0.73./tanh(0.331*delta.^1.01))).^-0.37;
Hs=4*sqrt(wind.^4.*epsilon/g^2);
Tp=wind./ni/g;