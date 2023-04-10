function dx = ODEoralGlucoseMinimalModel(t,x,p,ti,datIns,Gb)
% ODEoralGlucoseMinimalModel: Differential equations of the oral glucose minimal model
% Version 2 - changed Rate of appearance (Ra) function
%   dx = ODEoralGlucoseMinimalModel(t,x,p,datIns,Gb)
%
% input: 
% - t:  time
% - x:  vector of state variables (Q and X)
% - p:  parameter vector (SG, p2, p3, k, sigma, V, D)
% - datIns:  measured insulin concentrations determining the insulin input signal
% - Gb:  basal glucose level
% output:
% - dx:   column vector of derivatives

% parameters:
SG = p(1);
p2 = p(2);
p3 = p(3);
k = p(4);
sigma = p(5);
V = p(6);
D = p(7);

% Qb (plasma glucose mass) can be calculated from Gb (plasma glucose
% concentration) and V (glucose distribution volume)
Qb = Gb*V;

% state variables:
Q = x(1);
X = x(2);

% determine Ra(t) and I(t)
Ra = sigma.*k.^sigma.*t.^(sigma-1).*exp(-(k.*t).^sigma)*D;
I = interp1(ti,datIns,t,'linear');

% calculate derivatives based on ODEs:
dx = zeros(2,1);
dx(1) = -(SG + X)*Q+SG*Qb+Ra;
dx(2) = -p2*X+p3*(I-datIns(1));
end