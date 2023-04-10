function dx = ODEoralGlucoseMinimalModel(t,x,p,ti,datIns,Gb)
% ODEoralGlucoseMinimalModel: Differential equations of the oral glucose minimal model
%
%   dx = ODEoralGlucoseMinimalModel(t,x,p,datIns,Gb)
%
% input: 
% - t:  time
% - x:  vector of state variables (Q and X)
% - p:  parameter vector (SG, p2, p3, V, alpha1,..., alpha8)
% - datIns:  measured insulin concentrations determining the insulin input signal
% - Gb:  basal glucose level
% output:
% - dx:   column vector of derivatives

% parameters:
SG = p(1);
p2 = p(2);
p3 = p(3);
V = p(4);
alpha = [0,p(5:12)];

% Qb (plasma glucose mass) can be calculated from Gb (plasma glucose
% concentration) and V (glucose distribution volume)
Qb = Gb*V;

% state variables:
Q = x(1);
X = x(2);

% determine Ra(t) and I(t)
index = find(ti<=t,1,'last');
if (index < length(alpha))
    Ra = alpha(index)+(alpha(index+1)-alpha(index))/(ti(index+1)-ti(index))*(t-ti(index));
else
    Ra = 0;
end
I = interp1([-20,-10,ti],[datIns(1),datIns(1),datIns],t,'spline');

% calculate derivatives based on ODEs:
dx = zeros(2,1);
dx(1) = -(SG + X)*Q + SG*Qb + Ra;
dx(2) = -p2*X + p3*(I - datIns(1));
end