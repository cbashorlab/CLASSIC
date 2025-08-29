%% Random forest regression
a = readcell('Average_384_expression.csv');
expression = cell2mat(a(2:end, 1));
prom = categorical(cell2mat(a(2:end, 2)));
koz = categorical(cell2mat(a(2:end, 3)));
term = categorical(cell2mat(a(2:end, 4)));
X = table(prom, koz, term);

r = randsample(size(a, 1) - 1, floor(0.8*length(expression)));

x = expression; y = X;
etrain = x(r, :); Xtrain = y(r, :);
x(r, :) = []; y(r, :) = [];
etest = x; Xtest = y;

t = templateTree('NumVariablesToSample', 'all', 'predictorselection', 'interaction-curvature', 'surrogate', 'on');
Mdl = fitrensemble(Xtrain, etrain, 'Method','Bag', 'NumLearningCycles', 100, 'Learners', t);

yHat = resubPredict(Mdl); R2 = corr(Mdl.Y,yHat)^2;
ypred = predict(Mdl, Xtest);
R2pred = corr(ypred, etest)^2;

% Alternatively, load our pre-trained model and train-test splits 

a = readcell('RF_trainTestSplitValues.csv');
prom = categorical(cell2mat(a(:, 1)));
koz = categorical(cell2mat(a(:, 2)));
term = categorical(cell2mat(a(:, 3)));

etrain = cell2mat(a(1:307, 4)); Xtrain = table(prom(1:307), koz(1:307), term(1:307));
etest = cell2mat(a(308:end, 4)); Xtest = table(prom(308:end), koz(308:end), term(308:end));

%% Standard Regression fit on data

a = readcell('Average_384_expression.csv');
expression = [2; 2; 2; 2; cell2mat(a(2:end, 1))];
prom = [0; 0; 8; 8; cell2mat(a(2:end, 2))];
koz = [0; 5; 0; 5; cell2mat(a(2:end, 3))];
term = [0; 3; 3; 0; cell2mat(a(2:end, 4))];

syms p1 p2 p3 p4 p5 p6 p7 p8
syms k1 k2 k3 k4 k5 k6 k7 k8
syms t1 t2 t3 t4 t5 t6 t7 t8

zeta = p1*k1/t1 == expression(1);

for i = 2:24
    zeta(end+1) = eval(strcat('p', num2str(prom(i))))*eval(strcat('k', num2str(koz(i))))/eval(strcat('t', num2str(term(i)))) == expression(i);
end

x1 = prom;
x2 = koz;
x3 = term;

X = x1.*x2./x3;
b = regress(expression,X);

X = table(prom, koz, term, expression);

mdl = fitlm(X,'linear','ResponseVar','expression',...
    'PredictorVars',{'prom','koz','term'},...
    'CategoricalVar',{'prom','koz', 'term'});
