%% Base Data Processing

load('Single_input_indis_isolates.mat')
load('Single_input_measurements.mat')

%Indis
indi_ids = [individuals_96(:, 4:7) (individuals_96(:, 8) - 1)*3 + individuals_96(:, 9) individuals_96(:, 10:13)];
isolate_ids = [isolates_40(:, 4:7) (isolates_40(:, 8) - 1)*3 + isolates_40(:, 9) isolates_40(:, 10:13)];

%Combine the assign splits and measurements
assign_split_combined = eu_ordered_assign_split(ordered_eu_exp(:, 1) > 0 & eu_ordered_assign_split(:, 3) > 1, :);
fc_combined = ordered_eu_exp(ordered_eu_exp(:, 1) > 0 & eu_ordered_assign_split(:, 3) > 1, :);

ohes_combined = [];

for i = 1:size(assign_split_combined, 2)
    a = categorical(assign_split_combined(:, i));
    b = onehotencode(a, 2);
    ohes_combined(:, end+1:end+size(b, 2)) = b;
end

%Create the test set
x = ohes_combined(fc_combined(:, 4) > 12, :);
y = fc_combined(fc_combined(:, 4) > 12, :);
r = randsample(size(x, 1), 2000, false);
test_exp = y(r, :); test_ohes = x(r, :);

ohes_base = ohes_combined; exp_base = fc_combined;
ohes_base(ismember(ohes_combined, test_ohes, 'rows'), :) = [];
exp_base(ismember(ohes_combined, test_ohes, 'rows'), :) = [];

%Process the individuals (individuals_96 and isolates_40)
xcomb = [assign_split_combined; indi_ids; isolate_ids];
zeta = [];

for i = 1:size(xcomb, 2)
    a = categorical(xcomb(:, i));
    b = onehotencode(a, 2);
    zeta(:, end+1:end+size(b, 2)) = b;
end

indi_ohes = zeta(end-40-95:end-40, :);
iso_ohes = zeta(end-39:end);
clear zeta; clear xcomb

%% MLP

layers = [
    featureInputLayer(39)
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
r2s = zeros(60, 6);

for i = 1:60
    disp(i)
    x = ohes_base(exp_base(:, 4) > i-1, :);
    y = exp_base(exp_base(:, 4) > i-1, :);
    r = randsample(size(x, 1), floor(0.1.*size(x, 1)), false);
    val_exp = y(r, :); val_ohes = x(r, :);

    ohes_train = x; exp_train = y;
    ohes_train(ismember(x, val_ohes, 'rows'), :) = [];
    exp_train(ismember(x, val_ohes, 'rows'), :) = [];

    options = trainingOptions("sgdm", ...
        MaxEpochs=50, ...
        Verbose=false, ...
        LearnRateSchedule="piecewise", ...
        LearnRateDropFactor=0.5, ...
        LearnRateDropPeriod=7, ...
        Shuffle="every-epoch", ...
        Metrics="rmse", ...
        ValidationFrequency=10, ...
        Plots="training-progress", ...
        InitialLearnRate=5e-2, ...
        ValidationData=table(val_ohes, log10(val_exp(:, 1:2))), ...
        ValidationPatience=50, ...
        OutputNetwork="best-validation", ...
        MiniBatchSize=256, ...
        ExecutionEnvironment = "gpu");

    trainednet2 = trainnet(table(ohes_train, log10(exp_train(:, 1:2))), net2, "mse", options);
    pred_base = minibatchpredict(trainednet2, test_ohes);
    pred_val = minibatchpredict(trainednet2, val_ohes);
    pred_iso = minibatchpredict(trainednet2, indi_ohes);
    r2s(i, 1:2) = [corr(pred_base(:, 1), log10(test_exp(:, 1))) corr(pred_base(:, 2), log10(test_exp(:, 2)))].^2;
    r2s(i, 3:4) = [corr(pred_iso(:, 1), log10(individuals_96(:, 1))) corr(pred_iso(:, 2), log10(individuals_96(:, 2)))].^2;
    r2s(i, 5:6) = [corr(pred_val(:, 1), log10(val_exp(:, 1))) corr(pred_val(:, 2), log10(val_exp(:, 2)))].^2;
end

subplot(1, 2, 1)
scatter(test_exp(:, 1), 10.^pred_base(:, 1), 'filled')
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log')
xlim([100 200000]); ylim([100 200000]); axis square
ax = gca; ax.TickDir = 'out';
title(['r^2 = ', num2str(round(corr(log10(test_exp(:, 1)), pred_base(:, 1)).^2, 2))])

subplot(1, 2, 2)
scatter(individuals_96(:, 1), 10.^pred_iso(:, 1), 'filled')
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log')
xlim([100 200000]); ylim([100 200000]); axis square
ax = gca; ax.TickDir = 'out';
title(['r^2 = ', num2str(round(corr(log10(individuals_96(:, 1)), pred_iso(:, 1)).^2, 2))])

p = polyfit(log10(individuals_96(:, 1)), pred_iso(:, 1), 1);
pred_iso(:, 1) = (pred_iso(:, 1) - p(2))./p(1);
p = polyfit(log10(individuals_96(:, 2)), pred_iso(:, 2), 1);
pred_iso(:, 2) = (pred_iso(:, 2) - p(2))./p(1);
pred_iso(:, 3) = 10.^pred_iso(:, 2)./10.^pred_iso(:, 1);

p = polyfit(log10(test_exp(:, 1)), pred_base(:, 1), 1);
pred_base(:, 1) = (pred_base(:, 1) - p(2))./p(1);
p = polyfit(log10(test_exp(:, 2)), pred_base(:, 2), 1);
pred_base(:, 2) = (pred_base(:, 2) - p(2))./p(1);

err = abs(pred_base - log10(test_exp(:, 1:2)));

subplot(1, 2, 1)
scatter(test_exp(:, 1), 10.^pred_base(:, 1), 'filled')
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log')
xlim([100 200000]); ylim([100 200000]); axis square
ax = gca; ax.TickDir = 'out';
title(['r^2 = ', num2str(round(corr(log10(test_exp_new(:, 1)), pred_base(:, 1)).^2, 2))])

subplot(1, 2, 2)
scatter(test_exp_new(:, 2), 10.^pred_base(:, 2), 'filled')
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log')
xlim([100 200000]); ylim([100 200000]); axis square
ax = gca; ax.TickDir = 'out';
title(['r^2 = ', num2str(round(corr(log10(test_exp_new(:, 2)), pred_base(:, 2)).^2, 2))])

p = polyfit(log10(test_exp_new(:, 1)), pred_base(:, 1), 1);
p(2, :) = polyfit(log10(test_exp_new(:, 2)), pred_base(:, 2), 1);

%% CNN

%Process inputs differently
%Need to reshape ohes_base, test_ohes, indi_ohes
z = zeros(size(ohes_base, 1), 1);
[~, zeta] = find(ohes_base(:, 16:27) == 1);
zeta = floor((zeta - 1)./3) + 1;
zeta = onehotencode(categorical(zeta), 2);
ohes_base_cnn = [ohes_base(:, 1:15) z zeta ohes_base(:, 28:34) z ohes_base(:, 35:37) z ohes_base(:, 38:39) z z];
ohes_base_reshape = zeros(9, 4, 1, size(ohes_base_cnn, 1));

for i = 1:size(ohes_base_cnn, 1)
    ohes_base_reshape(:, :, :, i) = reshape(ohes_base_cnn(i, :), 4, 9)';
end

z = zeros(size(test_ohes, 1), 1);
[~, zeta] = find(test_ohes(:, 16:27) == 1);
zeta = [floor((zeta - 1)./3) + 1; 1; 2; 3; 4];
zeta = onehotencode(categorical(zeta), 2);
zeta = zeta(1:end-4, :);
test_ohes_cnn = [test_ohes(:, 1:15) z zeta test_ohes(:, 28:34) z test_ohes(:, 35:37) z test_ohes(:, 38:39) z z];
test_ohes_reshape = zeros(9, 4, 1, size(test_ohes_cnn, 1));
for i = 1:size(test_ohes_cnn, 1)
    test_ohes_reshape(:, :, :, i) = reshape(test_ohes_cnn(i, :), 4, 9)';
end

z = zeros(size(indi_ohes, 1), 1);
[~, zeta] = find(indi_ohes(:, 16:27) == 1);
zeta = [floor((zeta - 1)./3) + 1; 1; 2; 3; 4];
zeta = onehotencode(categorical(zeta), 2);
zeta = zeta(1:end-4, :);
indi_ohes_cnn = [indi_ohes(:, 1:15) z zeta indi_ohes(:, 28:34) z indi_ohes(:, 35:37) z indi_ohes(:, 38:39) z z];
indi_ohes_reshape = zeros(9, 4, 1, size(indi_ohes_cnn, 1));
for i = 1:size(indi_ohes_cnn, 1)
    indi_ohes_reshape(:, :, :, i) = reshape(indi_ohes_cnn(i, :), 4, 9)';
end

inputSize = [9, 4, 1];

layers = [
    imageInputLayer(inputSize)
    convolution2dLayer(4,100)
    batchNormalizationLayer
    reluLayer
    convolution2dLayer(1,100)
    batchNormalizationLayer
    reluLayer
    flattenLayer
    fullyConnectedLayer(20)
    tanhLayer
    fullyConnectedLayer(2)
    reluLayer];

net2 = dlnetwork(layers);
r2s_cnn = zeros(40, 4);

for i = 1:40
    disp(i)
    x = ohes_base_reshape(:, :, :, exp_base(:, 4) > i-1);
    y = exp_base(exp_base(:, 4) > i-1, :);
    r = randsample(size(x, 4), floor(0.1.*size(x, 4)), false);
    val_exp = y(r, :); val_ohes = x(:, :, :, r);

    ohes_train = x; exp_train = y;
    ohes_train(:, :, :, r) = [];
    exp_train(r, :) = [];

    options = trainingOptions("sgdm", ...
        MaxEpochs=50, ...
        Verbose=false, ...
        LearnRateSchedule="piecewise", ...
        LearnRateDropFactor=0.5, ...
        LearnRateDropPeriod=5, ...
        Shuffle="every-epoch", ...
        Metrics="rmse", ...
        ValidationFrequency=10, ...
        InitialLearnRate=5e-2, ...
        ValidationData={val_ohes, log10(val_exp(:, 1:2))}, ...
        ValidationPatience=50, ...
        OutputNetwork="best-validation", ...
        MiniBatchSize=256, ...
        ExecutionEnvironment = "gpu");

    trainednet2 = trainnet({ohes_train, log10(exp_train(:, 1:2))}, net2, "mse", options);
    pred_base = minibatchpredict(trainednet2, test_ohes_reshape);
    pred_iso = minibatchpredict(trainednet2, indi_ohes_reshape);
    pred_val = minibatchpredict(trainednet2, val_ohes);
    r2s_cnn(i, 1:2) = [corr(pred_base(:, 1), log10(test_exp(:, 1))) corr(pred_base(:, 2), log10(test_exp(:, 2)))].^2;
    r2s_cnn(i, 3:4) = [corr(pred_iso(:, 1), log10(individuals_96(:, 1))) corr(pred_iso(:, 2), log10(individuals_96(:, 2)))].^2;
    r2s_cnn(i, 5:6) = [corr(pred_val(:, 1), log10(val_exp(:, 1))) corr(pred_val(:, 2), log10(val_exp(:, 2)))].^2;
end

p = polyfit(log10(individuals_96(:, 1)), pred_iso(:, 1), 1);
pred_iso(:, 1) = (pred_iso(:, 1) - p(2))./p(1);
p = polyfit(log10(individuals_96(:, 2)), pred_iso(:, 2), 1);
pred_iso(:, 2) = (pred_iso(:, 2) - p(2))./p(1);

pred_iso(:, 3) = 10.^pred_iso(:, 2)./10.^pred_iso(:, 1);

%% RF

assign_split_rf_train = assign_split_combined(ismember(ohes_combined, ohes_base, 'rows'), :);
exp_rf_train = fc_combined(ismember(ohes_combined, ohes_base, 'rows'), :);

assign_split_rf_test = zeros(size(test_ohes, 1), 9);
exp_rf_test = zeros(size(assign_split_rf_test, 1), 4);
for i = 1:size(assign_split_rf_test, 1)
    a = test_ohes(i, :);
    assign_split_rf_test(i, :) = assign_split_combined(sum(ohes_combined == a, 2) == 39, :);
    exp_rf_test(i, :) = fc_combined(sum(ohes_combined == a, 2) == 39, :);
end
assign_split_rf_indi = indi_ids;

x = [assign_split_rf_train; assign_split_rf_test; assign_split_rf_indi];
i = [ones(size(assign_split_rf_train, 1), 1); 2*ones(size(assign_split_rf_test, 1), 1); 3*ones(size(assign_split_rf_indi, 1), 1)];

x = categorical(x);
x = table(x(:, 1), x(:, 2), x(:, 3), x(:, 4), x(:, 5), x(:, 6), x(:, 7), x(:, 8), x(:, 9));
assign_split_rf_train = x(i == 1, :);
assign_split_rf_test = x(i == 2, :);
assign_split_rf_indi = x(i == 3, :);

num_comps = zeros(60, 1);
for i = 1:60
    num_comps(i) = sum(fc_combined(:, 4) > i-1);
end

r2s_rf = zeros(40, 6);
for i = 1:40
    disp(i)
    x = assign_split_rf_train(exp_rf_train(:, 4) > i-1, :);
    y = exp_rf_train(exp_rf_train(:, 4) > i-1, :);
    r = randsample(size(x, 1), floor(0.1*size(x, 1)), false);
    x_val = x(r, :); y_val = y(r, :);
    x(r, :) = []; y(r, :) = [];
    t = templateTree('NumVariablesToSample', 5,'Surrogate','on', 'predictorSelection', 'interaction-curvature', 'MinLeafSize', 14, 'MaxNumSplits', 113);
    trained_rf = fitrensemble(x, log10(y(:, 1)), 'Learners', t, 'method', 'LSBoost', 'numLearningCycles', 495, 'LearnRate', 0.016939);
    t = templateTree('NumVariablesToSample', 5,'Surrogate','on', 'predictorSelection', 'interaction-curvature', 'MinLeafSize', 14, 'MaxNumSplits', 113);
    trained_rf2 = fitrensemble(x, log10(y(:, 2)), 'Learners', t, 'method', 'LSBoost', 'numLearningCycles', 495, 'LearnRate', 0.016939);
    pred_base = predict(trained_rf, assign_split_rf_test);
    pred_val = predict(trained_rf, x_val);
    pred_indi = predict(trained_rf, assign_split_rf_indi);
    pred_base(:, 2) = predict(trained_rf2, assign_split_rf_test);
    pred_val(:, 2) = predict(trained_rf2, x_val);
    pred_indi(:, 2) = predict(trained_rf2, assign_split_rf_indi);
    r2s_rf(i, 1) = [corr(pred_base(:, 1), log10(test_exp(:, 1)))].^2;
    r2s_rf(i, 3) = [corr(pred_indi(:, 1), log10(individuals_96(:, 1)))].^2;
    r2s_rf(i, 5) = [corr(pred_val(:, 1), log10(y_val(:, 1)))].^2;
end

p = polyfit(log10(individuals_96(:, 1)), pred_iso(:, 1), 1);
pred_iso(:, 1) = (pred_indi(:, 1) - p(2))./p(1);
p = polyfit(log10(individuals_96(:, 2)), pred_iso(:, 2), 1);
pred_iso(:, 2) = (pred_indi(:, 2) - p(2))./p(1);

pred_iso(:, 3) = 10.^pred_indi(:, 2)./10.^pred_iso(:, 1);

%% Linear Regression

assign_split_rf_train = assign_split_combined(ismember(ohes_combined, ohes_base, 'rows'), :);
exp_rf_train = fc_combined(ismember(ohes_combined, ohes_base, 'rows'), :);

assign_split_rf_test = zeros(size(test_ohes, 1), 9);
exp_rf_test = zeros(size(assign_split_rf_test, 1), 4);
for i = 1:size(assign_split_rf_test, 1)
    a = test_ohes(i, :);
    assign_split_rf_test(i, :) = assign_split_combined(sum(ohes_combined == a, 2) == 39, :);
    exp_rf_test(i, :) = fc_combined(sum(ohes_combined == a, 2) == 39, :);
end
assign_split_rf_indi = indi_ids;

x = [assign_split_rf_train; assign_split_rf_test; assign_split_rf_indi];
i = [ones(size(assign_split_rf_train, 1), 1); 2*ones(size(assign_split_rf_test, 1), 1); 3*ones(size(assign_split_rf_indi, 1), 1)];

x = categorical(x);
assign_split_rf_train = x(i == 1, :);
assign_split_rf_test = x(i == 2, :);
assign_split_rf_indi = x(i == 3, :);

r2s_lm = zeros(40, 6);
r2s_logm = zeros(40, 6);
for i = 1:40
    disp(i)
    x = assign_split_rf_train(exp_rf_train(:, 4) > i-1, :);
    y = exp_rf_train(exp_rf_train(:, 4) > i-1, :);
    r = randsample(size(x, 1), floor(0.1*size(x, 1)), false);
    x_val = x(r, :); y_val = y(r, :);
    x(r, :) = []; y(r, :) = [];
    
    yfit = table(x(:, 1), x(:, 2), x(:, 3), x(:, 4), x(:, 5), x(:, 6), x(:, 7), x(:, 8), x(:, 9), log10(y(:, 1)));
    lm = fitlm(yfit);
    lm2 = fitlm(yfit, 'interactions');
    pred_base = predict(lm, assign_split_rf_test);
    pred_val = predict(lm, x_val);
    pred_indi = predict(lm, assign_split_rf_indi);
    pred_base(:, 3) = predict(lm2, assign_split_rf_test);
    pred_val(:, 3) = predict(lm2, x_val);
    pred_indi(:, 3) = predict(lm2, assign_split_rf_indi);
    yfit = table(x(:, 1), x(:, 2), x(:, 3), x(:, 4), x(:, 5), x(:, 6), x(:, 7), x(:, 8), x(:, 9), log10(y(:, 2)));
    lm = fitlm(yfit);
    lm2 = fitlm(yfit, 'interactions');
    pred_base(:, 2) = predict(lm, assign_split_rf_test);
    pred_val(:, 2) = predict(lm, x_val);
    pred_indi(:, 2) = predict(lm, assign_split_rf_indi);
    pred_base(:, 4) = predict(lm2, assign_split_rf_test);
    pred_val(:, 4) = predict(lm2, x_val);
    pred_indi(:, 4) = predict(lm2, assign_split_rf_indi);
    r2s_lm(i, 1:2) = [corr(pred_base(:, 1), log10(test_exp(:, 1))) corr(pred_base(:, 2), log10(test_exp(:, 2)))].^2;
    r2s_lm(i, 3:4) = [corr(pred_indi(:, 1), log10(individuals_96(:, 1))) corr(pred_indi(:, 2), log10(individuals_96(:, 2)))].^2;
    r2s_lm(i, 5:6) = [corr(pred_val(:, 1), log10(y_val(:, 1))) corr(pred_val(:, 2), log10(y_val(:, 2)))].^2;
    r2s_logm(i, 1:2) = [corr(pred_base(:, 3), log10(test_exp(:, 1))) corr(pred_base(:, 4), log10(test_exp(:, 2)))].^2;
    r2s_logm(i, 3:4) = [corr(pred_indi(:, 3), log10(individuals_96(:, 1))) corr(pred_indi(:, 4), log10(individuals_96(:, 2)))].^2;
    r2s_logm(i, 5:6) = [corr(pred_val(:, 3), log10(y_val(:, 1))) corr(pred_val(:, 4), log10(y_val(:, 2)))].^2;

end

p = polyfit(log10(individuals_96(:, 1)), pred_indi(:, 1), 1);
pred_indi(:, 1) = (pred_indi(:, 1) - p(2))./p(1);
p = polyfit(log10(individuals_96(:, 2)), pred_indi(:, 2), 1);
pred_indi(:, 2) = (pred_indi(:, 2) - p(2))./p(1);
pred_indi(:, 3) = 10.^pred_indi(:, 2)./10.^pred_iso(:, 1);

p = polyfit(log10(individuals_96(:, 1)), pred_indi(:, 3), 1);
pred_indi(:, 3) = (pred_indi(:, 3) - p(2))./p(1);
p = polyfit(log10(individuals_96(:, 2)), pred_indi(:, 4), 1);
pred_indi(:, 4) = (pred_indi(:, 4) - p(2))./p(1);
pred_indi(:, 5) = 10.^pred_indi(:, 4)./10.^pred_indi(:, 3);
