% mp/mpscripts/bp_neuron01.m
% 8 nov 2013

% Calculation of a space clamped axon at 6.3 degC
% Hobbie p195

% CALLS alpha.m beta.m

%close all
clear all
clc

global Vr

% FIXED PARAMETERS =======================================================
dt = 1e-6;         % time increment
VR = -65e-3;       % resting voltage (V)
Vr = -65;          % resting voltage (mV)
VNa = 50e-3;       % reversal voltage for Na+ (V); default = 50e-3
VK = -77e-3;       % reversal voltage for K+ (V); default = -77e-3
Cm = 1e-6;         % membrane capacitance/area  (F.m^-2)

tmin = 0;          % starting time
tmax = 20e-3;       % finishing time (s)  default  5e-3

gKmax = 36e-3;     % K+ conductance (ohm^-1.cm^-2); default = 36e-3
gNamax = 120e-3;   % Na+ conductance (ohm^-1.cm.-2); default = 120e-3
gLmax = 0.3e-3;    % max leakage conductance (ohm-1.cm-2); default = 0.3e-3

Jext_max = 2e-4;   % max current density for ext stimulus (A.cm^-2; default 1e-4

ts = 0.5e-3;       % stimulus ON
tf = 0.6e-3;       % stimulus OFF
sf = 1e3;          % scale factor for consersion  v to mV and s to ms
T = 18.5;          % temperature (deg C) default 18.5

fs = 14;

% SETUP ==================================================================
num = 5000;

t = linspace(tmin,tmax,num);
dt = t(2) - t(1);

num1 = min(find(t > ts));       % index for stimulus ON
num2 = min(find(t > tf));       % index for stimulus OFF

Jext = zeros(num,1);       % external current density (A.cm^-2)
JNa  = zeros(num,1);       % Na+ current density (A.cm^-2)
JK   = zeros(num,1);       % K+  current density (A.cm^-2)
JL   = zeros(num,1);       % leakage current density (A.cm^-2)
Jm   = zeros(num,1);       % membrane current (A.cm^-2)
V    = zeros(num,1);       % membrane potential (V)
gNa  = zeros(num,1);       % Na+ conductance
gK   = zeros(num,1);       % K+ conductance
gL   = ones(num,1);        % gL conductance
n    = zeros(num,1);       % K+ gate parameter
m    = zeros(num,1);       % Na+ gate parameter
h    = zeros(num,1);       % Na+ gate parameter

V(1) = VR;                   % initial value for membrane potential

%for assignment 2: one pulse
% extext(num1:num2) = xtxt_max;  % external stimulus current

% for assignment 3: Setup Jext as a series of pulses
num_start = 0.05 * num;     
num_end   = 0.05 * num;
num_d     = 0.10 * num;      %controls the time interval between two pulses
num_width = round(0.01 * num); %controls the duration of a pulse

cn_max = 4;       %controls the number of pulses that are generated
num_index = zeros(cn_max,1);

for cn = 1 : cn_max
    num_index(cn) = round(num_start + cn * num_d);
    Jext(num_index(cn) : num_index(cn) + num_width) = Jext_max;
end


% Initial Values
[ An Am Ah ] = alpha(V(1)*1000, T);    % voltage in mV
[ Bn Bm Bh ] = beta(V(1)*1000, T);

n(1) = An / (An + Bn);
m(1) = Am / (Am + Bm);
h(1) = Ah / (Ah + Bh);

gK(1)  = gKmax * n(1)^4;
gNa(1) = gNamax * m(1)^3 * h(1);
gL = gLmax .* gL;

JK(1)  = gK(1)  * (V(1) - VK);
JNa(1) = gNa(1) * (V(1) - VNa);
JL(1)  = gL(1) * (V(1) - VR - 10.6e-3);
Jm(1)  = JNa(1) + JK(1) + JL(1);

V(1) = VR + (dt/Cm) * (-JK(1) - JNa(1) - JL(1) + Jext(1));

for cc = 1 : num-1
    
[ An Am Ah ] = alpha(V(cc)*1000, T);
[ Bn Bm Bh ] = beta(V(cc)*1000, T);
An = sf * An;   Am = sf * Am;   Ah = sf * Ah;  
Bn = sf * Bn;   Bm = sf * Bm;   Bh = sf * Bh; 

n(cc+1) = n(cc) + dt * (An *(1-n(cc)) - Bn * n(cc)); 
m(cc+1) = m(cc) + dt * (Am *(1-m(cc)) - Bm * m(cc)); 
h(cc+1) = h(cc) + dt * (Ah *(1-h(cc)) - Bh * h(cc)); 

gK(cc+1) = n(cc+1)^4 * gKmax;
gNa(cc+1) = m(cc+1)^3 * h(cc+1) * gNamax;

JK(cc+1)  = gK(cc+1)  * (V(cc) - VK);
JNa(cc+1) = gNa(cc+1) * (V(cc) - VNa);
JL(cc+1)  = gL(cc+1) * (V(cc) - VR - 10.6e-3);
Jm(cc+1)  = JNa(cc+1) + JK(cc+1) + JL(cc+1);

V(cc+1) = V(cc) + (dt/Cm) * (-JK(cc+1) - JNa(cc+1) - JL(cc+1) + Jext(cc+1));

end


figure(1)
clf;
set(gcf,'units','normalized');
set(gcf,'position',[0.1 0.1 0.8 0.8]);

subplot(3,1,1)
set(gca,'fontsize',fs);
title_x = 'time  t   (ms)';   title_y = 'J   (mA.cm ^{-2})';

x = t.*sf;   y = Jext.*sf;
plot(x,y,'linewidth',2);   %  Current - ext
xlabel(title_x); ylabel(title_y);

hold on
x = t.*sf;   y = JNa.*sf;
plot(x,y,'r','linewidth',2);   %  Current - Na+

x = t.*sf;   y = JK.*sf;
plot(x,y,'m','linewidth',2);   %  Current - K+

x = t.*sf;   y = JL.*sf;
plot(x,y,'c','linewidth',2);   %  Current - leakage

x = t.*sf;   y = Jm.*sf;
plot(x,y,'k','linewidth',2);   %  Current - K+

h_L = legend('J_{ext}','J_{Na}','J_K','J_L','J_m');
ax = gca;
ax.YLim = [-1 1];

subplot(3,1,2)
set(gca,'fontsize',fs);
title_x = 'time  t   (ms)';   title_y = '  g  ( mmho.cm^{-2})';

x = t.*sf;   y = gNa.*sf;
plot(x,y,'r','linewidth',2);   % conductance  Na+
hold on
x = t.*sf;   y = gK.*sf;
plot(x,y,'m','linewidth',2);   % conductance  K+

xlabel(title_x); ylabel(title_y);
grid on
legend('g_{Na}','g_K')
ax = gca;
ax.YLim = [0 30];

subplot(3,1,3)
set(gca,'fontsize',fs);
title_x = 'time  t   (ms)';   title_y = ' V_m (mV)';

x = t.*sf;   y = V.*sf;
plot(x,y,'linewidth',2);   % membrane voltage
xlabel(title_x); ylabel(title_y);
grid on
ax = gca;
ax.YLim = [-80 40];
