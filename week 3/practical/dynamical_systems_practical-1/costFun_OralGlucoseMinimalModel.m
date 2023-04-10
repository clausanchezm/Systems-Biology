function error = costFun_OralGlucoseMinimalModel(p,pRest,ti,datGlu,datIns,Gb)

ptemp = [p,pRest];
[~,Y] = ode15s(@ODEoralGlucoseMinimalModel,ti,[Gb*ptemp(6), 0],'',ptemp,ti,datIns,Gb);
error = (Y(:,1)./ptemp(6) - datGlu')./datGlu';
end