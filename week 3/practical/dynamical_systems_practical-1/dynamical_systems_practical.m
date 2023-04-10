%% LOTKA - VOLTERRA

% First, we look at the shark-tuna population dynamics model that we built
% intuition for during the lecture. What do we need in order to simulate
% this system?

% (1) We need the model that will tell the system where to go in
% state-space. The model in matlab is written in the file "lotka_volterra.m"
% (2) We need the timeframe of the simulation (including the starting and
% end time). The time frame of the simulation is represented by the start
% and end times in the 2-element vector "tspan"
% (3) We need an inital state (the number of tuna and sharks at t=0). The
% initial state is represented by a 2-element vector "y0" containing the
% initial numbers of tuna and shark.
% (4) We need to choose values for our parameters alpha, and beta)
% representing the "predation success rate", and the "growth rate of sharks
% due to predation", respectively.

% These 4 things (model, inital values, parameters) will determine what
% happens with this system at any point in the future. 
% (the system is deterministic!)

t0 = 0;                 	% start time of simulation [in months]
tfinal = 15;                % end time of simulation [in months]
y0 = [20 20];               % initial state i.e. number of tuna and sharks at t=0
tspan = [t0 tfinal];        % time span defined by the start and end times of the simulation [in months]
alpha = 0.03;               % predation success rate parameter ("how likely that the tuna dies if it meets a shark")
beta = 0.02;                % shark population increase rate due to predation on tuna ("how much food a shark gets out of eating a tuna")

% The system of ordindary differential equations for sharks and tuna are
% solved numerically by taking tiny steps forward in time and estimating
% the solution.

[t,y] = ode23(@(t,y) lotka_volterra(t,y,alpha,beta), tspan, y0);

% plotting the population over time
figure("Name",'Shark-Tuna')
plot(t,y)
title('Shark/Tuna Populations Over Time')
xlabel('Time [months]')
ylabel('Population')
legend('Tuna','Shark','Location','North')

% plotting the state space
figure("Name",'State space')
plot(y(:,1),y(:,2))
title('State space')
xlabel('Tuna Population')
ylabel('Shark Population')

% [TASK1: (1) Describe the population dynamics in your own words i.e. what
% happens to the number of tuna and number of sharks as we observe this
% system over time? 
% 
%
% (2) What is the maximum number of tuna and sharks that you can find in
% this patch of ocean?  
% 110 tuna and 73 sharks
% (3) Approximately what is the period of oscillations? 
% approx 6 timesteps 
% (4) Assume that you are part of a crew of fisherman hunting for sharks,
% how can such a model help you in your job? (look at the plot in time,
% when should you sail out to catch sharks?)]  
%3 and 10 month
% [ANSWER1:...]   

%% How does the system change if we change the value of alpha?

% Dynamic models often involve the use of kinetic/rate parameters that
% represent the rate at which particular processes of the system take
% place. In our simple example we have two parameters alpha and beta. Let's
% look at what happens if alpha takes on different values as we used in the
% previous simulation while we keep beta fixed.

[t,y_original] = ode23(@(t,y) lotka_volterra(t,y,alpha,beta), tspan, y0); % simulate with original parameters

figure("Name",'sensitivity_alpha')

original_alpha = plot(t,y_original, 'Color','red', 'LineWidth', 2); % plot the original simulation 

n = 201; % specify number of different parameter values we will generate
Salpha = linspace(alpha*0.5,alpha*1.5,n); % generate a vector of n parameters going from 50% to 150% of original value with equidistant spacing
colors = colormap(parula(n)); % generate n colors from a colormap to allow coloring them according to parameter value later on

% simulate and plot the model with the n paramneter values that we created
for i=1:n
    hold on
    [t,y] = ode23(@(t,y) lotka_volterra(t,y,Salpha(i),beta), tspan, y0); % simulate, note: Salpha(1) changes in each iteration!
    plot(t,y, 'Color', colors(i,:)) % plot the simulation in a color that represent the parameter value at this iteration
end

cb = colorbar('Ticks',[0,0.5,1],'TickLabels',[Salpha(1), Salpha(101), Salpha(n)]); % add colorbar to plot
ylabel(cb, 'alpha [predation success rate]');
uistack(original_alpha,'top'); % shows the original simulation on top

title('effect of changing \alpha')
xlabel('Time [months]')
ylabel('Population [# animals]')

% [TASK2: (1) Describe in your own words what happens to the shark-tuna
% dynamics if the value of alpha changes?] 
%By having a lower rate of succes in predation, tuna's population increases
% to 190 and sharks to 140, comparingto initial value, (110, 70) ; more or
% less keepig the approwx the same ratio, but when alpha increases teh
% population of sharks actually reduces to 45 and tunas are at 70; showing
% that even if sharks eat more tunas, they are less because they wont have
% enouh food
% [ANSWER2:...] 

%% How does the system change if we change the value of beta?

[t,y_original] = ode23(@(t,y) lotka_volterra(t,y,alpha,beta), tspan, y0); % simulate with original parameters

figure("Name",'sensitivity_beta')
original_beta = plot(t,y_original, 'Color','red', 'LineWidth', 2); % plot the original simulation 

n = 201; % specify number of different parameter values we will generate
Sbeta = linspace(beta*0.5,beta*1.5,n); % generate a vector of n parameters going from 50% to 150% of original value with equidistant spacing
colors = colormap(parula(n)); % generate n colors from a colormap to allow coloring them according to parameter value later on

% simulate and plot the model with the n paramneter values that we created
for i=1:n
    hold on
    [t,y] = ode23(@(t,y) lotka_volterra(t,y,alpha,Sbeta(i)), tspan, y0); % simulate, note: Sbeta(1) changes in each iteration!
    plot(t,y, 'Color', colors(i,:)) % plot the simulation in a color that represent the parameter value at this iteration
end

cb = colorbar('Ticks',[0,0.5,1],'TickLabels', [Sbeta(1), Sbeta(101) Sbeta(n)]); % add colorbar to plot
ylabel(cb, 'beta [shark pop. growth due to catching a tuna]');
uistack(original_beta,'top'); % shows the original simulation on top

title('effect of changing \beta')
xlabel('Time [months]')
ylabel('Population [# animals]')


% [TASK3: (1) Describe in your own words what happens to the shark-tuna
% dynamics if the value of beta changes?] 
%if we increase how  filling tunas are for sharks, there are going to be
%less sharks eating tunas, whereas if you decrease it, tunas will not be
%eaten that luch increasing theri pop, and sharks popo will also increase
%as tehre is a lot of available tiùme
% [ANSWER3:...]

%% How does the state-space look like as we change alpha?

[t,y_original] = ode23(@(t,y) lotka_volterra(t,y,alpha,beta), tspan, y0); % simulate with original parameters

figure("Name",'Sensitivity \alpha state-space')
original_alpha = plot(y_original(:,1),y_original(:,2), 'Color','red', 'LineWidth', 2); % plot the original simulation 

n = 201; % specify number of different parameter values we will generate
Salpha = linspace(alpha*0.5,alpha*1.5,n); % generate a vector of n parameters going from 50% to 150% of original value with equidistant spacing
colors = colormap(parula(n)); % generate n colors from a colormap to allow coloring them according to parameter value later on

% simulate and plot the model with the n paramneter values that we created
for i=1:n
    hold on
    [t,y] = ode23(@(t,y) lotka_volterra(t,y,Salpha(i),beta), tspan, y0); % simulate, note: Sbeta(1) changes in each iteration!
    plot(y(:,1),y(:,2),'Color',colors(i,:))
end

cb = colorbar('Ticks',[0,0.5,1],'TickLabels', [Salpha(1), Salpha(101) Salpha(n)]); % add colorbar to plot
ylabel(cb, '\alpha [shark pop. growth due to catching a tuna]');
uistack(original_alpha,'top'); % shows the original simulation on top

title('State space')
xlabel('Tuna Population')
ylabel('Shark Population')

%% How does the state-space look like as we change beta?

[t,y_original] = ode23(@(t,y) lotka_volterra(t,y,alpha,beta), tspan, y0); % simulate with original parameters

figure("Name",'Sensitivity \beta state-space')
original_beta = plot(y_original(:,1),y_original(:,2), 'Color','red', 'LineWidth', 2); % plot the original simulation 

n = 201; % specify number of different parameter values we will generate
Sbeta = linspace(beta*0.5,beta*1.5,n); % generate a vector of n parameters going from 50% to 150% of original value with equidistant spacing
colors = colormap(parula(n)); % generate n colors from a colormap to allow coloring them according to parameter value later on

% simulate and plot the model with the n paramneter values that we created
for i=1:n
    hold on
    [t,y] = ode23(@(t,y) lotka_volterra(t,y,alpha,Sbeta(i)), tspan, y0); % simulate, note: Sbeta(1) changes in each iteration!
    plot(y(:,1),y(:,2),'Color',colors(i,:))
end

cb = colorbar('Ticks',[0,0.5,1],'TickLabels', [Sbeta(1), Sbeta(101) Sbeta(n)]); % add colorbar to plot
ylabel(cb, '\beta [shark pop. growth due to catching a tuna]');
uistack(original_beta,'top'); % shows the original simulation on top

title('State space')
xlabel('Tuna Population')
ylabel('Shark Population')
%changing alpha , making teh predation succes increases the tuna nad shark
%pop or decreases both as sharks get goof at predation being mor scaled in
%teh vertical axis 
%whereas chaning beta, how fulling a tuna is, by increasing it we see taht
%the shark pop increase faster than the tuna but both pops decrease,
%decreasing beta as we saw amkes thetuna pop increase very quick along with
%teh sharks, an horizontal strech 

%% Shark fishing

% A fishing boat goes out to our patch of ocean to hunt sharks at t=10. We
% will assume that the fisherman will catch 40 sharks (that is the capacity
% of their fishingboat). What do you think will happen to the tuna and
% shark populations? 

alpha = 0.01;               % predation success rate parameter ("how likely that the tuna dies if it meets a shark")
beta = 0.02;                % shark population increase rate due to predation on tuna ("how much food a shark gets out of eating a tuna")

t0 = 0;                 	% start time of simulation
tboat = 10;                 % time of shark fishing boat arriving
y0 = [20 20];               % initial state i.e. number of tuna and sharks at t=0
tspan1 = [t0 tboat];        % simulation time from start to arrival of fishing boat

% solve the model from start to fishing:
[t1,y1] = ode23(@(t,y) lotka_volterra(t,y,alpha,beta), tspan1, y0);

tfinal = 15;                % end time of simulation
boat_capacity = 40;         % capacity of the fishing boat (=number of sharks that will be caught)
yboat = [y1(end,1) y1(end,2)-boat_capacity];    % initial condition at the time of fishing (where the system is at tboat minus the sharks that are caught)
tspan2 = [0 tfinal];        % simulation time from boat to end

[t2,y2] = ode23(@(t,y) lotka_volterra(t,y,alpha,beta), tspan2, yboat);

t = vertcat(t1,t2+t1(end));
y = vertcat(y1,y2);

% plotting the population over time
figure("Name",'Shark-Tuna')
plot(t,y)
title('Shark/Tuna Populations Over Time')
xlabel('Time [months]')
ylabel('Population')
legend('Tuna','Shark','Location','North')

% plotting the state space
figure("Name",'State space')
plot(y(:,1),y(:,2))
title('State space')
xlabel('Tuna Population')
ylabel('Shark Population')

% [TASK4: Describe in words, what happens to the shark-tuna dynamics after
% the fisherman catch the sharks at t=10.]
%after the fisherman caught those 40 sharks, (seen as a vertical line drop)
%teh tuna pop increasies and tehn sharks pop skyrockets as there is much
%more available tuna and less sharks so thery can reproduce mor 
% [ANSWER4:...]

%% Pen and paper practice exercise:
% [TASK 5:
% Look at the time-series and the state-space representation of the
% shark-tuna model. These representations of change in a system are
% conceptually interchangeble. Take look at a state-space and try to draw
% the time-series of the sharks and tuna. Then look at the time series and
% make a drawing of the state space! By the end of this practical you
% should be able to translate simple dynamics (like the Lotka-Volterra)
% from time domain representation to stat-space representation! If you have
% trouble, ask the instructors!]

%% BACKGROUND: The oral glucose minimal model

% The oral glucose minimal model allows us to approximate insulin
% sensitivity, your cell's ability to take glucose from circulation.
% Insulin action is a key player in the development of T2DM which decreases
% as insulin sensitivity declines over the years. Therefore, it's important
% to measure it. The oral glucose minimal model approximates insulin
% sensitivity from data from a frequeantly sampled oral glucose tolerance
% test (OGTT). We give the model the insulin measurements (or interpolated
% insulin measurements) and ask it to predict the glucose measurements. If
% it predicts the glucose measurements accurately, we can consider the
% model parameters to be a sort of "summary" of the OGTT response of the
% participant. We will look at a scenario where we cannot interpret the
% parameter values and another scenario where we can later on.

% But first, let's look at what an OGTT test looks like by plotting  
% the data over time.

ti = [0,15,30,45,60,90,120,180,240]; % Timepoints
datIns = [0.04339 0.14100 0.41646 0.51149 0.41065 0.25279 0.14389 0.07045 0.04854]; % Insulin measurements
datGlu = [98.8121, 104.9252, 117.6979, 124.3251, 116.0987, 105.7905, 100.6042, 94.3310, 87.5437]; % Glucose measurements

figure('pos',[200 200 1000 400]);
subplot(121) 
plot(ti,datGlu,'rx','markersize',12,'linewidth',2); % plot the simulated glucose concentration
xlabel('time [min]','fontsize',16)
ylabel('Glucose [mg/dL]','fontsize',16)
set(gca,'linewidth',2,'fontsize',16,'xlim',[-5 245],'xtick',0:60:240)
legend('raw data')
legend(gca,'boxoff')

subplot(122) 
plot(ti,datIns,'rx','markersize',12,'linewidth',2); % plot the insulin measurements
hold on;
t_interp = linspace(0,240,241);
i_interp = interp1(ti,datIns,t_interp,'spline'); % Interpolation
plot(t_interp,i_interp,'k--','linewidth',2); % plot the interpolated data
xlabel('time [min]','fontsize',16)
ylabel('Insulin [nmol/L]','fontsize',16)
set(gca,'linewidth',2,'fontsize',16,'xlim',[-5 245],'xtick',0:60:240)
legend('raw data','interpolated data')
legend(gca,'boxoff')

ti = [0,15,30,45,60,90,120,180,240]; % Timepoints
datGlu = [157.9039, 159.4060, 179.8218, 204.3239, 218.0918, 220.0738, 205.1847, 167.7431, 141.5470];
datIns = [0.03339 0.13100 0.31646 0.41149 0.41065 0.25279 0.14389 0.07045 0.04854]; %Insulin Datapoints

figure('pos',[200 200 1000 400]);
subplot(121) 
plot(ti,datGlu,'rx','markersize',12,'linewidth',2); % plot the simulated glucose concentration
xlabel('time [min]','fontsize',16)
ylabel('Glucose [mg/dL]','fontsize',16)
set(gca,'linewidth',2,'fontsize',16,'xlim',[-5 245],'xtick',0:60:240)
legend('raw data')
legend(gca,'boxoff')

subplot(122) 
plot(ti,datIns,'rx','markersize',12,'linewidth',2); % plot the insulin measurements
hold on;
t_interp = linspace(0,240,241);
i_interp = interp1(ti,datIns,t_interp,'spline'); % Interpolation
plot(t_interp,i_interp,'k--','linewidth',2); % plot the interpolated data
xlabel('time [min]','fontsize',16)
ylabel('Insulin [nmol/L]','fontsize',16)
set(gca,'linewidth',2,'fontsize',16,'xlim',[-5 245],'xtick',0:60:240)
legend('raw data','interpolated data')
legend(gca,'boxoff')

%% SIMULATING THE ORAL GLUCOSE MINIMAL MODEL

Gb=98.8121; % Fasting Glucose value
ti = [0,15,30,45,60,90,120,180,240]; % Timepoints
datIns = [0.04339 0.14100 0.41646 0.51149 0.41065 0.25279 0.14389 0.07045 0.04854]; % Insulin measurements
datGlu = [98.8121, 104.9252, 117.6979, 124.3251, 116.0987, 105.7905, 100.6042, 94.3310, 87.5437]; % Glucose measurements

p = [5, 1, 0.07, 0.04, 2.5, 17, 40000]; % [SG, p2, p3, k, sigma, V, D] --> for ease we created this parameter vector

options = '';
[T,Y] = ode15s(@ODEoralGlucoseMinimalModel,[0,240],[Gb*p(6), 0],options,p,ti,datIns,Gb);

figure('pos',[400 400 500 400]); % Plot the data time and measured glucose
plot(ti,datGlu,'rx','markersize',12,'linewidth',2); hold on; % Glucose measurements
plot(T,Y(:,1)/p(6),'k-','linewidth',2); % Glucose prediction/model output
xlabel('time [min]','fontsize',16)
ylabel('Glucose [mg/dL]','fontsize',16)
set(gca,'linewidth',2,'fontsize',16,'xlim',[-5 245],'xtick',0:60:240)
legend('raw data','simulated curve')
legend(gca,'boxoff')

% calculation of SI, insulin sensitivity:
SI_nonoptimized = (p(3)/p(2)) * p(6)

% [TASK6: Think back to the lecture. (1) What does SI represent? why is it
% interesting/important? How much insulin
% (2) Can we trust this particular estimate for SI? Take a look at 
% the corresponding model fit that you generated with lines 214-221. 
% (3) Does the model fit your data well? no 
% (4) What does this mean for your parameter estimates for p3 and p2 and by
% extension to the value of SI? 

% [ANSWER6:...]

%% Parameter estimation (calculate SI from data of a normoglycaemic person)

% In the previous section we did not estimate the parameters from data. We
% only chose the values and simulated the model. This is why the glucose
% measurements were so far from the model simulated glucose values. In this
% section we estimate the model parameters to minimize the squared error on
% between measured and model simulated glucose values. For this we need two
% things that we haven't looked at yet: (1) we need a cost function that
% calculates the error between model simulation and measured data. This is
% contained in the file "costFun_OralGlucoseMinimalModel.m". Take a look at
% this function!

% Running this code section will print text into your Command Window, this
% text contains information about the iterative optimization process. It
% tells you which "iteration" you are on as you are trying to find the
% optimal parameter set to describe your data. It also tells you the number
% of function calls ("Func-count" i.e. how many times the solver was
% called) the value of the cost function "f(x)" etc. The optimizer runs
% until it has found the optimal parameter set or did not converge to an
% uptimum after a certain number of tries. The optimal parameter set is
% stored in the output res.x.

Gb=98.8121; % Fasting Glucose value
ti = [0,15,30,45,60,90,120,180,240]; %Timepoints
datIns = [0.13170, 0.15222, 0.29688, 0.38983, 0.39742, 0.35831, 0.27523, 0.17372, 0.11100];
datGlu = [98.8121, 104.9252, 117.6979, 124.3251, 116.0987, 105.7905, 100.6042, 94.3310, 87.5437]; %Glucose Datapoints

p = [5, 1, 0.07, 0.03, 2.5, 17, 50000];  % [SG, p2, p3, k, sigma, V, D] --> for ease we created this parameter vector
% optimize the 5 parameters (see costFun_OralGlucoseMinimalModel.m):
lb = [0.01, 0, 0, 0.001, 1]; % lower bound of parameter search
ub = [100, 10, 10, 1, 10]; % upper bound of parameter search
options = optimset('Display', 'Iter'); % we tell the optimizer to display results from every iteration
[res.x,res.resnorm,res.residual,res.exitflag,res.output,res.lambda,res.jacobian] = lsqnonlin(@costFun_OralGlucoseMinimalModel,p(1:5),lb,ub,options,p(6:7),ti,datGlu,datIns,Gb);
popt = [res.x, p(6:7)];
[T,Y] = ode15s(@ODEoralGlucoseMinimalModel,[0,240],[Gb*p(6), 0],'',popt,ti,datIns,Gb); % solve the model using the optimal parameter values to plot accurate simulation

% plot the results below
figure('pos',[400 400 500 400]);
plot(ti,datGlu,'rx','markersize',12,'linewidth',2); hold on; %Glucose measurements
plot(T,Y(:,1)/p(6),'k-','linewidth',2); % Glucose prediction/model output
xlabel('time [min]','fontsize',16)
ylabel('Glucose [mg/dL]','fontsize',16)
set(gca,'linewidth',2,'fontsize',16,'xlim',[-5 245],'xtick',0:60:240)
legend('raw data','simulated curve')
legend(gca,'boxoff')

% calculation of SI, insulin sensitivity:
SI_1 = (res.x(3)/res.x(2)) * p(6)

%% Parameter estimation (calculate SI from data of a diabetic person)

p = [0.5, 1, 0.07, 0.04, 2.5, 17, 50000];  
ti = [0,15,30,45,60,90,120,180,240]; %Timepoints
datGlu2 = [157.9039, 159.4060, 179.8218, 204.3239, 218.0918, 220.0738, 205.1847, 167.7431, 141.5470];
datIns2 = [0.03339 0.13100 0.31646 0.41149 0.41065 0.25279 0.14389 0.07045 0.04854]; %Insulin Datapoints
Gb2= datGlu2(1);
lb = [0.01, 0, 0, 0.001, 1];
ub = [100, 10, 10, 1, 10];
options = optimset('Display', 'Iter');
[res.x,res.resnorm,res.residual,res.exitflag,res.output,res.lambda,res.jacobian] = lsqnonlin(@costFun_OralGlucoseMinimalModel,p(1:5),lb,ub,options,p(6:7),ti,datGlu2,datIns2,Gb2);
popt = [res.x, p(6:7)];
[T,Y] = ode15s(@ODEoralGlucoseMinimalModel,[0,240],[Gb2*p(6), 0],'',popt,ti,datIns2,Gb2);

% plot the results below
figure('pos',[400 400 500 400]);
plot(ti,datGlu2,'rx','markersize',12,'linewidth',2); hold on; %Glucose measurements
plot(T,Y(:,1)/p(6),'k-','linewidth',2); % Glucose prediction/model output
xlabel('time [min]','fontsize',16)
ylabel('Glucose [mg/dL]','fontsize',16)
set(gca,'linewidth',2,'fontsize',16,'xlim',[-5 245],'xtick',0:60:240)
legend('raw data','simulated curve')
legend(gca,'boxoff')

% calculation of SI, insulin sensitivity:
SI_2 = (res.x(3)/res.x(2)) * p(6)

% Now that we estimated the SI from individuals' data, we can
% interpret it.
% 
% [TASK7: How does SI_1 compare to SI_2? Which one is higher?
% What does that mean physiologically? Which of these individuals could be
% diabetic?]

% [ANSWER7:...]

