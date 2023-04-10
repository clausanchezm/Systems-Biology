inputs = repmat(input', 1, 100);
target_outputs= repmat(target_output', 1, 100);
mynet.trainFcn = 'trainscg'
mynet.trainParam.showWindo = 0; 
trained_mynet = train(mynet, inputs, target_outputs)
trained_network_simulated_output = trained_mynet(input')'
disp('untrained network weighs (in hidden layer)');
mynet.IW{1,1}
disp('trianed network weights (in hidden layer)');
trained_mynet.IW{1,1}