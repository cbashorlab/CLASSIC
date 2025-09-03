
%Combine the assign splits and measurements
assign_split_combined = [assign_split_train; new_assign_split; assign_split_test];
fc_combined = [fc_new, ones(size(fc_new, 1), 1); fc_new_new(:, [1:2 4]); test_set, 10*ones(size(test_set, 1), 1)];
indices = [ones(size(assign_split_train, 1), 1); 2*ones(size(new_assign_split, 1), 1); 3*ones(size(assign_split_test, 1), 1)];

ohes_combined = [];

for i = 1:size(assign_split_combined, 2)
    a = categorical(assign_split_combined(:, i));
    b = onehotencode(a, 2);
    ohes_combined(:, end+1:end+size(b, 2)) = b;
end

ohes_train_base = ohes_combined(indices == 1, :); exp_train_base = fc_combined(indices == 1, :);
ohes_new = ohes_combined(indices == 2, :); exp_new = fc_combined(indices == 2, :);
ohes_test_base = ohes_combined(indices == 3, :); exp_test_base = fc_combined(indices == 3, :);

r = randsample(size(ohes_train_base, 1), floor(0.1*size(ohes_train_base, 1)), false);
ohes_val_base = ohes_train_base(r, :); ohes_train_base(r, :) = [];
exp_val_base = exp_train_base(r, :); exp_train_base(r, :) = [];

%Split ohes_new into No IDR (number 1) & NFZ (AD number 5)

ohes_noIDR = ohes_new(ohes_new(:, 10) == 1, :); exp_noIDR = exp_new(ohes_new(:, 10) == 1, :);
%Create test set now
ohes_temp = ohes_noIDR(exp_noIDR(:, 3) > 50, :); exp_temp = exp_noIDR(exp_noIDR(:, 3) > 50, :);
r = randsample(size(ohes_temp, 1), floor(0.5*size(ohes_temp, 1)), false);
ohes_noIDR_test = ohes_temp(r, :); exp_noIDR_test = exp_temp(r, :);

r = randsample(size(exp_noIDR, 1), floor(0.8*size(exp_noIDR, 1)), false);
ohes_noIDR_train = ohes_noIDR(r, :); exp_noIDR_train = exp_noIDR(r, :);
ohes_noIDR_train(ismember(exp_noIDR_train(:, 1), exp_noIDR_test(:, 1), 'rows'), :) = [];
exp_noIDR_train(ismember(exp_noIDR_train(:, 1), exp_noIDR_test(:, 1), 'rows'), :) = [];
ohes_noIDR_val = ohes_noIDR; ohes_noIDR_val(r, :) = [];
exp_noIDR_val = exp_noIDR; exp_noIDR_val(r, :) = [];
ohes_noIDR_val(ismember(exp_noIDR_val(:, 1), exp_noIDR_test(:, 1), 'rows'), :) = [];
exp_noIDR_val(ismember(exp_noIDR_val(:, 1), exp_noIDR_test(:, 1), 'rows'), :) = [];


%% Do the same splits as above for NFZ
ohes_NFZ = ohes_new(ohes_new(:, 9) == 1, :); exp_NFZ = exp_new(ohes_new(:, 9) == 1, :);
ohes_temp = ohes_NFZ(exp_NFZ(:, 3) > 50, :); exp_temp = exp_NFZ(exp_NFZ(:, 3) > 50, :);
r = randsample(size(ohes_temp, 1), floor(0.5*size(ohes_temp, 1)), false);
ohes_NFZ_test = ohes_temp(r, :); exp_NFZ_test = exp_temp(r, :);

r = randsample(size(exp_NFZ, 1), floor(0.8*size(exp_NFZ, 1)), false);
ohes_NFZ_train = ohes_NFZ(r, :); exp_NFZ_train = exp_NFZ(r, :);
ohes_NFZ_train(ismember(exp_NFZ_train(:, 1), exp_NFZ_test(:, 1), 'rows'), :) = [];
exp_NFZ_train(ismember(exp_NFZ_train(:, 1), exp_NFZ_test(:, 1), 'rows'), :) = [];
ohes_NFZ_val = ohes_noIDR; ohes_NFZ_val(r, :) = [];
exp_NFZ_val = exp_noIDR; exp_NFZ_val(r, :) = [];
ohes_NFZ_val(ismember(exp_NFZ_val(:, 1), exp_NFZ_test(:, 1), 'rows'), :) = [];
exp_NFZ_val(ismember(exp_NFZ_val(:, 1), exp_NFZ_test(:, 1), 'rows'), :) = [];

%% Train the base model on the base data

layers = [
    featureInputLayer(41)
    fullyConnectedLayer(160)
    tanhLayer
    fullyConnectedLayer(80)
    tanhLayer
    fullyConnectedLayer(40)
    tanhLayer
    fullyConnectedLayer(20)
    tanhLayer
    fullyConnectedLayer(2)
    reluLayer];

net2 = dlnetwork(layers);

options = trainingOptions("sgdm", ...
    MaxEpochs=50, ...
    Verbose=false, ...
    LearnRateSchedule="piecewise", ...
    LearnRateDropFactor=0.5, ...
    LearnRateDropPeriod=7, ...
    Plots="training-progress", ...
    Shuffle="every-epoch", ...
    Metrics="rmse", ...
    ValidationFrequency=10, ...
    InitialLearnRate=5e-2, ...
    ValidationData=table(ohes_val_base, log10(exp_val_base(:, 1:2))), ...
    ValidationPatience=50, ...
    OutputNetwork="best-validation", ...
    MiniBatchSize=256, ...
    ExecutionEnvironment = "gpu");

trainednet2 = trainnet(table(ohes_train_base, log10(exp_train_base(:, 1:2))), net2, "mse", options);

pred_base = minibatchpredict(trainednet2, ohes_test_base);
pred_noIDR = minibatchpredict(trainednet2, ohes_noIDR_test);
pred_NFZ = minibatchpredict(trainednet2, ohes_NFZ_test);

figure
subplot(1, 3, 1)
scatter(exp_test_base(:, 1), 10.^pred_base(:, 1), 'filled')
set(gca, 'YScale', 'log'); set(gca, 'XScale', 'log')
axis square; title(num2str(corr(log10(exp_test_base(:, 1)), pred_base(:, 1)).^2))
xlim([100 10^5.5]); ylim([100 10^5.5])

subplot(1, 3, 2)
scatter(exp_noIDR_test(:, 1), 10.^pred_noIDR(:, 1), 'filled')
set(gca, 'YScale', 'log'); set(gca, 'XScale', 'log')
axis square; title(num2str(corr(log10(exp_noIDR_test(:, 1)), pred_noIDR(:, 1)).^2))
xlim([100 10^5.5]); ylim([100 10^5.5])

subplot(1, 3, 3)
scatter(exp_NFZ_test(:, 1), 10.^pred_NFZ(:, 1), 'filled')
set(gca, 'YScale', 'log'); set(gca, 'XScale', 'log')
axis square; title(num2str(corr(log10(exp_NFZ_test(:, 1)), pred_NFZ(:, 1)).^2))
xlim([100 10^5.5]); ylim([100 10^5.5])

%% Fine tune the trained network

r = randsample(size(ohes_val_base, 1), size(ohes_noIDR_val, 1), false);
options = trainingOptions("sgdm", ...
    MaxEpochs=100, ...
    Verbose=false, ...
    LearnRateSchedule="piecewise", ...
    LearnRateDropFactor=0.5, ...
    LearnRateDropPeriod=7, ...
    Plots="training-progress", ...
    Shuffle="every-epoch", ...
    Metrics="rmse", ...
    ValidationFrequency=10, ...
    InitialLearnRate=5e-3, ...
    ValidationData=table([ohes_val_base(r, :); ohes_noIDR_val], [log10(exp_val_base(r, 1:2)); log10(exp_noIDR_val(:, 1:2))]), ...
    ValidationPatience=50, ...
    OutputNetwork="best-validation", ...
    MiniBatchSize=64, ...
    ExecutionEnvironment = "gpu");

r = randsample(size(ohes_train_base, 1), 5*size(ohes_noIDR_train, 1), false);
trainednet3 = trainnet(table([ohes_train_base(r, :); ohes_noIDR_train], [log10(exp_train_base(r, 1:2)); log10(exp_noIDR_train(:, 1:2))]), trainednet2, "mse", options);

pred_base = minibatchpredict(trainednet3, ohes_test_base);
pred_noIDR = minibatchpredict(trainednet3, ohes_noIDR_test);
pred_NFZ = minibatchpredict(trainednet3, ohes_NFZ_test);

figure
subplot(1, 3, 1)
scatter(exp_test_base(:, 1), 10.^pred_base(:, 1), 'filled')
set(gca, 'YScale', 'log'); set(gca, 'XScale', 'log')
axis square; title(num2str(corr(log10(exp_test_base(:, 1)), pred_base(:, 1)).^2))

subplot(1, 3, 2)
scatter(exp_noIDR_test(:, 1), 10.^pred_noIDR(:, 1), 'filled')
set(gca, 'YScale', 'log'); set(gca, 'XScale', 'log')
axis square; title(num2str(corr(log10(exp_noIDR_test(:, 1)), pred_noIDR(:, 1)).^2))

subplot(1, 3, 3)
scatter(exp_NFZ_test(:, 1), 10.^pred_NFZ(:, 1), 'filled')
set(gca, 'YScale', 'log'); set(gca, 'XScale', 'log')
axis square; title(num2str(corr(log10(exp_NFZ_test(:, 1)), pred_NFZ(:, 1)).^2))

options = trainingOptions("sgdm", ...
    MaxEpochs=100, ...
    Verbose=false, ...
    LearnRateSchedule="piecewise", ...
    LearnRateDropFactor=0.5, ...
    LearnRateDropPeriod=7, ...
    Plots="training-progress", ...
    Shuffle="every-epoch", ...
    Metrics="rmse", ...
    ValidationFrequency=10, ...
    InitialLearnRate=1e-3, ...
    ValidationData=table(ohes_NFZ_val, log10(exp_NFZ_val(:, 1:2))), ...
    ValidationPatience=50, ...
    OutputNetwork="best-validation", ...
    MiniBatchSize=64, ...
    ExecutionEnvironment = "gpu");

trainednet4 = trainnet(table(ohes_NFZ_train, log10(exp_NFZ_train(:, 1:2))), trainednet2, "mse", options);

pred_base = minibatchpredict(trainednet4, ohes_test_base);
pred_noIDR = minibatchpredict(trainednet4, ohes_noIDR_test);
pred_NFZ = minibatchpredict(trainednet4, ohes_NFZ_test);

figure
subplot(1, 3, 1)
scatter(exp_test_base(:, 1), 10.^pred_base(:, 1), 'filled')
set(gca, 'YScale', 'log'); set(gca, 'XScale', 'log')
axis square; title(num2str(corr(log10(exp_test_base(:, 1)), pred_base(:, 1)).^2))

subplot(1, 3, 2)
scatter(exp_noIDR_test(:, 1), 10.^pred_noIDR(:, 1), 'filled')
set(gca, 'YScale', 'log'); set(gca, 'XScale', 'log')
axis square; title(num2str(corr(log10(exp_noIDR_test(:, 1)), pred_noIDR(:, 1)).^2))

subplot(1, 3, 3)
scatter(exp_NFZ_test(:, 1), 10.^pred_NFZ(:, 1), 'filled')
set(gca, 'YScale', 'log'); set(gca, 'XScale', 'log')
axis square; title(num2str(corr(log10(exp_NFZ_test(:, 1)), pred_NFZ(:, 1)).^2))

%% Do a full sweep to figure out how much data to train with (for fine tuning), 
% and how much of validation data to keep to not "forget" base examples

base_set = [0, 1, 2, 5, 10, 20, floor(size(exp_train_base, 1)./size(exp_noIDR_train, 1))];
newpart_set = 0:0.1:1;

r2_basesweep_basal = zeros(length(base_set), length(newpart_set));
r2_newpart_basal = r2_basesweep_basal;
r2_basesweep_induced = zeros(length(base_set), length(newpart_set));
r2_newpart_induced = r2_basesweep_basal;

for i = 1:length(base_set)
    x = base_set(i);
    for j = 1:length(newpart_set)
        disp(i); disp(j)
        y = newpart_set(j);
        if x == 0 && y == 0
            pred_base = minibatchpredict(trainednet2, ohes_test_base);
            pred_noIDR = minibatchpredict(trainednet2, ohes_noIDR_test);
            r2_basesweep_basal(i, j) = corr(log10(exp_test_base(:, 1)), pred_base(:, 1)).^2;
            r2_newpart_basal(i, j) = corr(log10(exp_noIDR_test(:, 1)), pred_noIDR(:, 1)).^2;
            r2_basesweep_induced(i, j) = corr(log10(exp_test_base(:, 2)), pred_base(:, 2)).^2;
            r2_newpart_induced(i, j) = corr(log10(exp_noIDR_test(:, 2)), pred_noIDR(:, 2)).^2;
        else
        r = randsample(size(ohes_val_base, 1), size(ohes_noIDR_val, 1), false);
        vals = [ohes_val_base(r, :); ohes_noIDR_val];
        expv = [exp_val_base(r, :); exp_noIDR_val];
        if x > 0
            r1 = randsample(size(ohes_train_base, 1), x*size(ohes_noIDR_train, 1), false);
        end
        r2 = randsample(size(ohes_noIDR_train, 1), floor(y*size(ohes_noIDR_train, 1)), false);
        if x > 0
            trains = [ohes_train_base(r1, :); ohes_noIDR_train(r2, :)];
            expt = [exp_train_base(r1, :); exp_noIDR_train(r2, :)];
        else
            trains = ohes_noIDR_train(r2, :);
            expt = exp_noIDR_train(r2, :);
        end
        options = trainingOptions("sgdm", ...
            MaxEpochs=100, ...
            Verbose=false, ...
            LearnRateSchedule="piecewise", ...
            LearnRateDropFactor=0.5, ...
            LearnRateDropPeriod=7, ...
            Shuffle="every-epoch", ...
            Metrics="rmse", ...
            ValidationFrequency=10, ...
            InitialLearnRate=1e-3, ...
            ValidationData=table(vals, log10(expv(:, 1:2))), ...
            ValidationPatience=25, ...
            OutputNetwork="best-validation", ...
            MiniBatchSize=64, ...
            ExecutionEnvironment = "gpu");
        trainednet3 = trainnet(table(trains, log10(expt(:, 1:2))), trainednet2, "mse", options);
        pred_base = minibatchpredict(trainednet3, ohes_test_base);
        pred_noIDR = minibatchpredict(trainednet3, ohes_noIDR_test);
        r2_basesweep_basal(i, j) = corr(log10(exp_test_base(:, 1)), pred_base(:, 1)).^2;
        r2_newpart_basal(i, j) = corr(log10(exp_noIDR_test(:, 1)), pred_noIDR(:, 1)).^2;
        r2_basesweep_induced(i, j) = corr(log10(exp_test_base(:, 2)), pred_base(:, 2)).^2;
        r2_newpart_induced(i, j) = corr(log10(exp_noIDR_test(:, 2)), pred_noIDR(:, 2)).^2;
        end
    end
end

% Without any base examples in the validation set

r2_basesweep_noval_basal = zeros(length(base_set), length(newpart_set));
r2_newpart_noval_basal = r2_basesweep_basal;
r2_basesweep_noval_induced = zeros(length(base_set), length(newpart_set));
r2_newpart_noval_induced = r2_basesweep_basal;

for i = 1:length(base_set)
    x = base_set(i);
    for j = 1:length(newpart_set)
        disp(i); disp(j)
        y = newpart_set(j);
        r = randsample(size(ohes_val_base, 1), size(ohes_noIDR_val, 1), false);
        vals = ohes_noIDR_val;
        expv = exp_noIDR_val;
        if x == 0 && y == 0
            r2_basesweep_noval_basal(i, j) = r2_basesweep_basal(i, j);
            r2_newpart_noval_basal(i, j) = r2_newpart_basal(i, j);
            r2_basesweep_noval_induced(i, j) = r2_basesweep_induced(i, j);
            r2_newpart_noval_induced(i, j) = r2_newpart_induced(i, j);
        else
        if x > 0
            r1 = randsample(size(ohes_train_base, 1), x*size(ohes_noIDR_train, 1), false);
        end
        r2 = randsample(size(ohes_noIDR_train, 1), floor(y*size(ohes_noIDR_train, 1)), false);
        if x > 0
            trains = [ohes_train_base(r1, :); ohes_noIDR_train(r2, :)];
            expt = [exp_train_base(r1, :); exp_noIDR_train(r2, :)];
        else
            trains = ohes_noIDR_train(r2, :);
            expt = exp_noIDR_train(r2, :);
        end
        options = trainingOptions("sgdm", ...
            MaxEpochs=100, ...
            Verbose=false, ...
            LearnRateSchedule="piecewise", ...
            LearnRateDropFactor=0.5, ...
            LearnRateDropPeriod=7, ...
            Shuffle="every-epoch", ...
            Metrics="rmse", ...
            ValidationFrequency=10, ...
            InitialLearnRate=1e-3, ...
            ValidationData=table(vals, log10(expv(:, 1:2))), ...
            ValidationPatience=25, ...
            OutputNetwork="best-validation", ...
            MiniBatchSize=64, ...
            ExecutionEnvironment = "gpu");
        trainednet3 = trainnet(table(trains, log10(expt(:, 1:2))), trainednet2, "mse", options);
        pred_base = minibatchpredict(trainednet3, ohes_test_base);
        pred_noIDR = minibatchpredict(trainednet3, ohes_noIDR_test);
        r2_basesweep_noval_basal(i, j) = corr(log10(exp_test_base(:, 1)), pred_base(:, 1)).^2;
        r2_newpart_noval_basal(i, j) = corr(log10(exp_noIDR_test(:, 1)), pred_noIDR(:, 1)).^2;
        r2_basesweep_noval_induced(i, j) = corr(log10(exp_test_base(:, 2)), pred_base(:, 2)).^2;
        r2_newpart_noval_induced(i, j) = corr(log10(exp_noIDR_test(:, 2)), pred_noIDR(:, 2)).^2;
        end
    end
end

% One final test, with no validation and no training examples in the base
% set

r2_basesweep_novalnotrain_basal = zeros(length(base_set), length(newpart_set));
r2_newpart_novalnotrain_basal = r2_basesweep_novalnotrain_basal;
r2_basesweep_novalnotrain_induced = zeros(length(base_set), length(newpart_set));
r2_newpart_novalnotrain_induced = r2_basesweep_novalnotrain_basal;

for i = 1:length(base_set)
    x = base_set(i);
    for j = 2:length(newpart_set)
        disp(i); disp(j)
        y = newpart_set(j);
        r = randsample(size(ohes_val_base, 1), size(ohes_noIDR_val, 1), false);
        vals = ohes_noIDR_val;
        expv = exp_noIDR_val;
        if x == 0 && y == 0
            r2_basesweep_novalnotrain_basal(i, j) = r2_basesweep_basal(i, j);
            r2_newpart_novalnotrain_basal(i, j) = r2_newpart_basal(i, j);
            r2_basesweep_novalnotrain_induced(i, j) = r2_basesweep_induced(i, j);
            r2_newpart_novalnotrain_induced(i, j) = r2_newpart_induced(i, j);
        else
        r2 = randsample(size(ohes_noIDR_train, 1), floor(y*size(ohes_noIDR_train, 1)), false);
        trains = ohes_noIDR_train(r2, :);
        expt = exp_noIDR_train(r2, :);
        options = trainingOptions("sgdm", ...
            MaxEpochs=100, ...
            Verbose=false, ...
            LearnRateSchedule="piecewise", ...
            LearnRateDropFactor=0.5, ...
            LearnRateDropPeriod=7, ...
            Shuffle="every-epoch", ...
            Metrics="rmse", ...
            ValidationFrequency=10, ...
            InitialLearnRate=1e-3, ...
            ValidationData=table(vals, log10(expv(:, 1:2))), ...
            ValidationPatience=25, ...
            OutputNetwork="best-validation", ...
            MiniBatchSize=64, ...
            ExecutionEnvironment = "gpu");
        trainednet3 = trainnet(table(trains, log10(expt(:, 1:2))), trainednet2, "mse", options);
        pred_base = minibatchpredict(trainednet3, ohes_test_base);
        pred_noIDR = minibatchpredict(trainednet3, ohes_noIDR_test);
        r2_basesweep_novalnotrain_basal(i, j) = corr(log10(exp_test_base(:, 1)), pred_base(:, 1)).^2;
        r2_newpart_novalnotrain_basal(i, j) = corr(log10(exp_noIDR_test(:, 1)), pred_noIDR(:, 1)).^2;
        r2_basesweep_novalnotrain_induced(i, j) = corr(log10(exp_test_base(:, 2)), pred_base(:, 2)).^2;
        r2_newpart_novalnotrain_induced(i, j) = corr(log10(exp_noIDR_test(:, 2)), pred_noIDR(:, 2)).^2;
        end
    end
end

%% Now do the same search for NFZ

base_set = [0, 1, 2, 5, 10, 20, floor(size(exp_train_base, 1)./size(exp_NFZ_train, 1))];
newpart_set = 0:0.1:1;

r2_basesweep_basal = zeros(length(base_set), length(newpart_set));
r2_newpart_basal = r2_basesweep_basal;
r2_basesweep_induced = zeros(length(base_set), length(newpart_set));
r2_newpart_induced = r2_basesweep_basal;

for i = 1:length(base_set)
    x = base_set(i);
    for j = 1:length(newpart_set)
        disp(i); disp(j)
        y = newpart_set(j);
        if x == 0 && y == 0
            pred_base = minibatchpredict(trainednet3, ohes_test_base);
            pred_NFZ = minibatchpredict(trainednet3, ohes_NFZ_test);
            r2_basesweep_basal(i, j) = corr(log10(exp_test_base(:, 1)), pred_base(:, 1)).^2;
            r2_newpart_basal(i, j) = corr(log10(exp_NFZ_test(:, 1)), pred_NFZ(:, 1)).^2;
            r2_basesweep_induced(i, j) = corr(log10(exp_test_base(:, 2)), pred_base(:, 2)).^2;
            r2_newpart_induced(i, j) = corr(log10(exp_NFZ_test(:, 2)), pred_NFZ(:, 2)).^2;
        else
        r = randsample(size(ohes_val_base, 1), size(ohes_NFZ_val, 1), false);
        vals = [ohes_val_base(r, :); ohes_NFZ_val];
        expv = [exp_val_base(r, :); exp_NFZ_val];
        if x > 0
            r1 = randsample(size(ohes_train_base, 1), x*size(ohes_NFZ_train, 1), false);
        end
        r2 = randsample(size(ohes_NFZ_train, 1), floor(y*size(ohes_NFZ_train, 1)), false);
        if x > 0
            trains = [ohes_train_base(r1, :); ohes_NFZ_train(r2, :)];
            expt = [exp_train_base(r1, :); exp_NFZ_train(r2, :)];
        else
            trains = ohes_NFZ_train(r2, :);
            expt = exp_NFZ_train(r2, :);
        end
        options = trainingOptions("sgdm", ...
            MaxEpochs=100, ...
            Verbose=false, ...
            LearnRateSchedule="piecewise", ...
            LearnRateDropFactor=0.5, ...
            LearnRateDropPeriod=7, ...
            Shuffle="every-epoch", ...
            Metrics="rmse", ...
            ValidationFrequency=10, ...
            InitialLearnRate=1e-3, ...
            ValidationData=table(vals, log10(expv(:, 1:2))), ...
            ValidationPatience=50, ...
            OutputNetwork="best-validation", ...
            MiniBatchSize=64, ...
            ExecutionEnvironment = "gpu");
        trainednet4 = trainnet(table(trains, log10(expt(:, 1:2))), trainednet2, "mse", options);
        pred_base = minibatchpredict(trainednet4, ohes_test_base);
        pred_NFZ = minibatchpredict(trainednet4, ohes_NFZ_test);
        r2_basesweep_basal(i, j) = corr(log10(exp_test_base(:, 1)), pred_base(:, 1)).^2;
        r2_newpart_basal(i, j) = corr(log10(exp_NFZ_test(:, 1)), pred_NFZ(:, 1)).^2;
        r2_basesweep_induced(i, j) = corr(log10(exp_test_base(:, 2)), pred_base(:, 2)).^2;
        r2_newpart_induced(i, j) = corr(log10(exp_NFZ_test(:, 2)), pred_NFZ(:, 2)).^2;
        end
    end
end

% Without any base examples in the validation set
r2_basesweep_noval_basal = zeros(length(base_set), length(newpart_set));
r2_newpart_noval_basal = r2_basesweep_basal;
r2_basesweep_noval_induced = zeros(length(base_set), length(newpart_set));
r2_newpart_noval_induced = r2_basesweep_basal;

for i = 1:length(base_set)
    x = base_set(i);
    for j = 1:length(newpart_set)
        disp(i); disp(j)
        y = newpart_set(j);
        r = randsample(size(ohes_val_base, 1), size(ohes_NFZ_val, 1), false);
        vals = ohes_NFZ_val;
        expv = exp_NFZ_val;
        if x == 0 && y == 0
            r2_basesweep_noval_basal(i, j) = r2_basesweep_basal(i, j);
            r2_newpart_noval_basal(i, j) = r2_newpart_basal(i, j);
            r2_basesweep_noval_induced(i, j) = r2_basesweep_induced(i, j);
            r2_newpart_noval_induced(i, j) = r2_newpart_induced(i, j);
        else
        if x > 0
            r1 = randsample(size(ohes_train_base, 1), x*size(ohes_NFZ_train, 1), false);
        end
        r2 = randsample(size(ohes_NFZ_train, 1), floor(y*size(ohes_NFZ_train, 1)), false);
        if x > 0
            trains = [ohes_train_base(r1, :); ohes_NFZ_train(r2, :)];
            expt = [exp_train_base(r1, :); exp_NFZ_train(r2, :)];
        else
            trains = ohes_NFZ_train(r2, :);
            expt = exp_NFZ_train(r2, :);
        end
        options = trainingOptions("sgdm", ...
            MaxEpochs=100, ...
            Verbose=false, ...
            LearnRateSchedule="piecewise", ...
            LearnRateDropFactor=0.5, ...
            LearnRateDropPeriod=7, ...
            Shuffle="every-epoch", ...
            Metrics="rmse", ...
            ValidationFrequency=10, ...
            InitialLearnRate=1e-3, ...
            ValidationData=table(vals, log10(expv(:, 1:2))), ...
            ValidationPatience=25, ...
            OutputNetwork="best-validation", ...
            MiniBatchSize=64, ...
            ExecutionEnvironment = "gpu");
        trainednet4 = trainnet(table(trains, log10(expt(:, 1:2))), trainednet3, "mse", options);
        pred_base = minibatchpredict(trainednet4, ohes_test_base);
        pred_NFZ = minibatchpredict(trainednet4, ohes_NFZ_test);
        r2_basesweep_noval_basal(i, j) = corr(log10(exp_test_base(:, 1)), pred_base(:, 1)).^2;
        r2_newpart_noval_basal(i, j) = corr(log10(exp_NFZ_test(:, 1)), pred_NFZ(:, 1)).^2;
        r2_basesweep_noval_induced(i, j) = corr(log10(exp_test_base(:, 2)), pred_base(:, 2)).^2;
        r2_newpart_noval_induced(i, j) = corr(log10(exp_NFZ_test(:, 2)), pred_NFZ(:, 2)).^2;
        end
    end
end

% One final test, with no validation and no training examples in the base
% set

r2_basesweep_novalnotrain_basal = zeros(length(base_set), length(newpart_set));
r2_newpart_novalnotrain_basal = r2_basesweep_novalnotrain_basal;
r2_basesweep_novalnotrain_induced = zeros(length(base_set), length(newpart_set));
r2_newpart_novalnotrain_induced = r2_basesweep_novalnotrain_basal;

for i = 1:length(base_set)
    x = base_set(i);
    for j = 2:length(newpart_set)
        disp(i); disp(j)
        y = newpart_set(j);
        r = randsample(size(ohes_val_base, 1), size(ohes_NFZ_val, 1), false);
        vals = ohes_NFZ_val;
        expv = exp_NFZ_val;
        if x == 0 && y == 0
            r2_basesweep_novalnotrain_basal(i, j) = r2_basesweep_basal(i, j);
            r2_newpart_novalnotrain_basal(i, j) = r2_newpart_basal(i, j);
            r2_basesweep_novalnotrain_induced(i, j) = r2_basesweep_induced(i, j);
            r2_newpart_novalnotrain_induced(i, j) = r2_newpart_induced(i, j);
        else
        r2 = randsample(size(ohes_NFZ_train, 1), floor(y*size(ohes_NFZ_train, 1)), false);
        trains = ohes_NFZ_train(r2, :);
        expt = exp_NFZ_train(r2, :);
        options = trainingOptions("sgdm", ...
            MaxEpochs=100, ...
            Verbose=false, ...
            LearnRateSchedule="piecewise", ...
            LearnRateDropFactor=0.5, ...
            LearnRateDropPeriod=7, ...
            Shuffle="every-epoch", ...
            Metrics="rmse", ...
            ValidationFrequency=10, ...
            InitialLearnRate=1e-3, ...
            ValidationData=table(vals, log10(expv(:, 1:2))), ...
            ValidationPatience=25, ...
            OutputNetwork="best-validation", ...
            MiniBatchSize=64, ...
            ExecutionEnvironment = "gpu");
        trainednet4 = trainnet(table(trains, log10(expt(:, 1:2))), trainednet3, "mse", options);
        pred_base = minibatchpredict(trainednet4, ohes_test_base);
        pred_NFZ = minibatchpredict(trainednet4, ohes_NFZ_test);
        r2_basesweep_novalnotrain_basal(i, j) = corr(log10(exp_test_base(:, 1)), pred_base(:, 1)).^2;
        r2_newpart_novalnotrain_basal(i, j) = corr(log10(exp_NFZ_test(:, 1)), pred_NFZ(:, 1)).^2;
        r2_basesweep_novalnotrain_induced(i, j) = corr(log10(exp_test_base(:, 2)), pred_base(:, 2)).^2;
        r2_newpart_novalnotrain_induced(i, j) = corr(log10(exp_NFZ_test(:, 2)), pred_NFZ(:, 2)).^2;
        end
    end
end

%% Look at weight changes

x = trainednet2.Layers;
x = x(2).Weights;

y = trainednet3.Layers;
y = y(2).Weights;

%% Validation Individuals

indi_data_master = readcell('Single_input_Sublib_Indis.xlsx'); %Variable available upon request
indi_ids = [cell2mat(indi_data_master(3:14, 19:22)), ones(12, 1), cell2mat(indi_data_master(3:14, 26:28)) cell2mat(indi_data_master(3:14, 25))];
indi_exp = cell2mat(indi_data_master(3:14, [4 7]));
indi_exp(:, 3) = indi_exp(:, 2)./indi_exp(:, 1);

assign_split_indis = [assign_split_train; new_assign_split; assign_split_test; indi_ids];
indi_ohes = [];

for i = 1:size(assign_split_combined, 2)
    a = categorical(assign_split_indis(:, i));
    b = onehotencode(a, 2);
    indi_ohes(:, end+1:end+size(b, 2)) = b;
end

indi_ohes = indi_ohes(end-11:end, :);
pred_indis = minibatchpredict(trainednet4, indi_ohes);

p = polyfit(log10(indi_exp(:, 1)), pred_indis(:, 1), 1);
pred_indis(:, 1) = (pred_indis(:, 1) - p(2))./p(1);
p = polyfit(log10(indi_exp(:, 2)), pred_indis(:, 2), 1);
pred_indis(:, 2) = (pred_indis(:, 2) - p(2))./p(1);

pred_indi_fc = 10.^pred_indis(:, 2)./10.^pred_indis(:, 1);

%% Predict the whole space

%Generate the whole space first

cats = max(assign_split_combined);
assign_split_full_space = zeros(prod(cats), 9);
counter = 0;
for a = 1:max(cats(1))
    for b = 1:max(cats(2))
        for c = 1:max(cats(3))
            for d = 1:max(cats(4))
                for e = 1:max(cats(5))
                    for f = 1:max(cats(6))
                        for g = 1:max(cats(7))
                            for h = 1:max(cats(8))
                                for i = 1:max(cats(9))
                                    counter = counter + 1;
                                    assign_split_full_space(counter, :) = [a b c d e f g h i];
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

%ohes_full_space = zeros(size(assign_split_full_space, 1), sum(cats));
ohes_full_space = [];

for i = 1:size(cats, 2)
    a = onehotencode(categorical(assign_split_full_space(:, i)), 2);
    ohes_full_space(1:size(a, 1), end+1:end+size(a, 2)) = a;
end

pred_full_space = minibatchpredict(trainednet4, ohes_full_space);
pred_NFZ = minibatchpredict(trainednet4, ohes_NFZ_test);
pred_noIDR = minibatchpredict(trainednet4, ohes_noIDR_test);
pred_base = minibatchpredict(trainednet4, ohes_test_base);

completed_space = [assign_split_full_space, 10.^pred_full_space];

