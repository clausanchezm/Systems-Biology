%% Assignment Questions 4-6
clear; clf;


%% parameters of the model
dt=0.1; % integration time step [ms]
tau=10; % time constant [ms]
E_L=-65; % resting potential [mV]
theta=-55; % firing threshold [mV]
RI_ext=12; % constant external input [Ohm*mA=mV]

%% Integration with Euler method
t_step=0; v=E_L;
for t=0:dt:100;
    t_step=t_step+1;
    s=v>theta;
    v=s*E_L+(1-s)*(v-dt/tau*((v-E_L)-RI_ext));
    v_rec(t_step)=v;
    t_rec(t_step)=t;
    s_rec(t_step)=s;
end

%% Plotting results
figure,plot(t_rec,v_rec,'Linewidth',2);
hold on; plot([0 100],[-55 -55],'--','Linewidth',2);
axis([0 100 -66 -54]);
xlabel('Time [ms]'); ylabel('v [mV]')

