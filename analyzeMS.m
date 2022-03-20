function analyzeMS(outBase, dataFileXLS, aggregationType, non_inf_kd_int)
if (~exist('dataFileXLS', 'var') || isempty(dataFileXLS))
    dataFileXLS = 'ms_allData3.xlsx';
end
if (~exist('aggregationType', 'var'))
    aggregationType = 'mean';
end
if (~exist('non_inf_kd_int', 'var'))
    non_inf_kd_int = true;
end

D = 35.29; % high domain concentration    
D2 = 3.529; % low domain concentration
units = 10^-6; % concentration units in M
allResult = [];
show = 0;
groups = 5;

% 4096-peptide library
if ~isempty(find(groups == 1, 1))
    allResult = fitKds(struct('xlsFile', dataFileXLS, 'sheetName', 'trial3','region','B1:L2247', 'domainConc', D, 'show', show, 'format', 1), allResult);
    allResult = fitKds(struct('xlsFile', dataFileXLS, 'sheetName', 'trial3','region','S1:AC2218', 'domainConc', D, 'show', show, 'format', 1), allResult);
    allResult = fitKds(struct('xlsFile', dataFileXLS, 'sheetName', 'trial3','region','AJ1:AT2201', 'domainConc', D, 'show', show, 'format', 1), allResult);
end

%24-pep
if ~isempty(find(groups == 2, 1))
    allResult = fitKds(struct('xlsFile', dataFileXLS, 'sheetName', 'trial2','region','B1:K127','domainConc', D, 'show', show, 'format', 4), allResult);
    allResult = fitKds(struct('xlsFile', dataFileXLS, 'sheetName', 'trial2','region','M1:V109','domainConc', D, 'show', show, 'format', 4), allResult);
end

%432-pep
if ~isempty(find(groups == 3, 1))
    allResult = fitKds(struct('xlsFile', dataFileXLS, 'sheetName', 'trial9-set1','region','B1:L318', 'domainConc', D, 'show', show, 'format', 2), allResult);
    allResult = fitKds(struct('xlsFile', dataFileXLS, 'sheetName', 'trial9-set1','region','O1:Y354', 'domainConc', D, 'show', show, 'format', 2), allResult);
    allResult = fitKds(struct('xlsFile', dataFileXLS, 'sheetName', 'trial9-set1','region','AB1:AL337', 'domainConc', D, 'show', show, 'format', 2), allResult);
end

%864-pep
if ~isempty(find(groups == 4, 1))
    allResult = fitKds(struct('xlsFile', dataFileXLS, 'sheetName', 'trial9-set2','region','B1:L777', 'domainConc', D, 'show', show, 'format', 2), allResult);
    allResult = fitKds(struct('xlsFile', dataFileXLS, 'sheetName', 'trial9-set2','region','O1:Y800', 'domainConc', D, 'show', show, 'format', 2), allResult);
    allResult = fitKds(struct('xlsFile', dataFileXLS, 'sheetName', 'trial9-set2','region','AB1:AL836', 'domainConc', D, 'show', show, 'format', 2), allResult);
end

%432-pep(low domain concentration)
if ~isempty(find(groups == 5, 1))
    allResult = fitKds(struct('xlsFile', dataFileXLS, 'sheetName', 'trial8-set2','region','B1:L217', 'domainConc', D2, 'show', show,'format', 2), allResult);
    allResult = fitKds(struct('xlsFile', dataFileXLS, 'sheetName', 'trial8-set2','region','O1:Y206', 'domainConc', D2, 'show', show,'format', 2), allResult);
    allResult = fitKds(struct('xlsFile', dataFileXLS, 'sheetName', 'trial8-set2','region','AB1:AL194', 'domainConc', D2, 'show', show,'format', 2), allResult);
end


presentResults(allResult, struct('aggregationType', aggregationType, 'units', units, 'show', 1, 'non_inf_kd_int', non_inf_kd_int), outBase);
plotFpMsKd([outBase, '.csv']);


function presentResults(allResult, opts, outBase)
fprintf('Presenting results for %s\n', outBase);
if (~isfield(opts, 'show')), opts.show = 1; end

% now learning sequence-Kd mapping
M = containers.Map; % mean Kd for each unique sequence
A = containers.Map; % all Kds for each unique sequence
E = containers.Map; % error of each Kd estimate for each unique sequence
LE = containers.Map; % error of each logKd estimate for each unique sequence
F = containers.Map; % all alpha values for each unique peptide
G = containers.Map; % mean alpha values for each unique peptide
bM = containers.Map; % best estimate of Kd for each unique sequence
bE = containers.Map; % error associated with the best estimate of Kd for each unique sequence

J = containers.Map; % all Xcorr values
K = containers.Map; % mean Xcorr values

for j = 1:length(allResult)
    for i = 1:size(allResult(j).Kd, 1) 
        idx = isfinite(allResult(j).Kd(i, :)); 
        kds = allResult(j).Kd(i, idx);
        kd_err = allResult(j).KdError(i, idx);
        logkd_err = allResult(j).logKdError(i, idx);
        kd_ints = allResult(j).KdHi(i, idx) - allResult(j).KdLo(i, idx);
        f = allResult(j).a(i,idx);
        xcorr = allResult(j).xcor(i,idx);
        pass = opts.non_inf_kd_int && all(~isempty(kds)) && all(~isinf(kd_ints)) && all(~isinf(kd_err)) && all(~isnan(kds));
        pass = pass || (~opts.non_inf_kd_int && all(~isempty(kds)) && all(~isnan(kds)));
        if pass
            if (isKey(M, allResult(j).seqs{i}))
                A(allResult(j).seqs{i}) = [A(allResult(j).seqs{i}) kds];
                M(allResult(j).seqs{i}) = mean(A(allResult(j).seqs{i}));
                E(allResult(j).seqs{i}) = [E(allResult(j).seqs{i}) kd_err];
                LE(allResult(j).seqs{i}) = [LE(allResult(j).seqs{i}) logkd_err];
                F(allResult(j).seqs{i}) = [F(allResult(j).seqs{i}) f];
                G(allResult(j).seqs{i}) = mean(F(allResult(j).seqs{i}));
             
                [~, si] = min(kd_err);
                if (kd_err(si) < bE(allResult(j).seqs{i}))
                    bM(allResult(j).seqs{i}) = kds(si);
                    bE(allResult(j).seqs{i}) = kd_err(si);
                end

                J(allResult(j).seqs{i}) = [J(allResult(j).seqs{i}) xcorr];
                K(allResult(j).seqs{i}) = mean(J(allResult(j).seqs{i}));
            else
                M(allResult(j).seqs{i}) = mean(kds);
                A(allResult(j).seqs{i}) = kds;
                E(allResult(j).seqs{i}) = kd_err;
                LE(allResult(j).seqs{i}) = logkd_err;
                F(allResult(j).seqs{i}) = f;
                G(allResult(j).seqs{i}) = mean(f);

                [~, si] = min(kd_err);
                bM(allResult(j).seqs{i}) = kds(si);
                bE(allResult(j).seqs{i}) = kd_err(si);

                J(allResult(j).seqs{i}) = xcorr;
                K(allResult(j).seqs{i}) = mean(xcorr); 
            end
        end
    end
    fprintf('%d unique peptides \n', length(keys(M)));
end


% write summary
fid = fopen(sprintf('%s.csv', outBase), 'w');
fprintf(fid, 'sequence, Kd_estimate, error(obs), error (est),number of samples, alpha, xcorr\n');
keySeqs = keys(M);
if strcmp(opts.aggregationType, 'mean')
    vals = values(M);
    vals = [vals{:}];
elseif strcmp(opts.aggregationType, 'best')
    vals = values(bM);
    vals = [vals{:}];
else
    error('unrecognized aggregation type "%s"', opts.aggregationType);
end
alphaValue = values(G);
alphaValue = [alphaValue{:}];

xVal = values(K);
xVal = [xVal{:}];
errTab = zeros(length(keySeqs), 4);
for i = 1:length(keySeqs)
    err_est = mean(E(keySeqs{i}));
    err_obs = std(A(keySeqs{i}));

    errTab(i, :) = [err_est, err_obs, length(A(keySeqs{i})), vals(i)];
    fprintf(fid, '%s,%f,%f,%f,%f,%f,%f\n', keySeqs{i}, vals(i), err_obs, bE(keySeqs{i}), length(A(keySeqs{i})), alphaValue(i), xVal(i));
end
fclose(fid);

% plot estimated versus observed error
if opts.show
    figure;
    p = errTab(errTab(:, 3) > 4, 1);
    q = errTab(errTab(:, 3) > 4, 2);
    fprintf('Observed errors computed for %d peptides\n', length(p));
    plot(log10(q/1000000), log10(p/1000000), 'o','color','black');
    ylabel('log10(Estimated error/[M])');
    xlabel('log10(Observed error/[M])');

    figure;
    plot(q,p,'o','color','black')
    ylabel('Estimated error/[uM]');
    xlabel('Observed error/[uM]');
end

r = corrcoef(p,q); fprintf('R between observed and estimated errors: %f\n', r(1,2));
r = corrcoef(log10(p),log10(q)); fprintf('R between logs of observed and estimated errors: %f\n', r(1,2));

% file a site-independent model
[mm, alpha] = modelMatrix(keySeqs);
b = fitlm(mm(:, 2:end), log10(vals') + log10(opts.units)); 
ci = coefCI(b,.01); 

m1 = ci(:,1);
m2 = ci(:,2);
m_hat = table2array(b.Coefficients(1:end, 1));
m_error = (m2-m1)/2;

if opts.show
    figure;
    plot(log10(vals) + log10(opts.units), mm*m_hat, 'o');
    set(gca, 'FontSize', 14);
    xlabel('experimental log10(Kd [M])');
    ylabel('linear-model prediction');
end

cor3 = corrcoef(log10(vals) + log10(opts.units), mm*m_hat)

k = 1;
fprintf('--- parameters ---\n');
fprintf('const = %f\t%f\n', m_hat(1), m_error(1));
for i = 1:length(alpha)
    fprintf('--> site %d\n', i);
    for j = 1:length(alpha{i})
        % NOTE: we are going to absorb the intercept term into each position's energies
        if (j == 1)
            fprintf('\t\t%s\t0.0\n', alpha{i}(j));
        else
            fprintf('\t\t%s\t%.2f\t%.2f\n', alpha{i}(j), m_hat(k+1), m_error(k+1));
            k = k + 1;
        end
    end
end

function results = fitKds(inputs, results)
% read data
[num, ~, raw] = xlsread(inputs.xlsFile, inputs.sheetName, inputs.region);
allSeqs = raw(2:end, 7);
 
% collect data from the three repeats
switch inputs.format
    case 1
        xcorr = num(:, 4);
        expOut = num(:, 9);
        ctrIn = num(:, 11);
        ctrOut = num(:, 10);
    case 2
        xcorr = num(:, 4);
        expOut = num(:, 9);
        ctrIn = num(:, 10);
        ctrOut = num(:, 11);
    case 4
        % collect data from the three repeats
        xcorr = num(:, 5);% put three repeats in three columns
        expOut = num(:, 8);
        ctrIn = num(:, 10);
        ctrOut = num(:, 9);
    otherwise
        error(sprintf('unrecognized format: %d', inputs.format));
end

data = struct('seqs', {allSeqs}, 'xcorr', xcorr, 'expOut', expOut, 'ctrIn', ctrIn, 'ctrOut', ctrOut);
ret = fitKdsFromData(inputs, data);
if isempty(results)
    results = ret;
else
    results(end+1) = ret;
end


function result = fitKdsFromData(inputs, data)

% mark which sequences are okay
seqsOK = false(1, length(data.seqs));
for i = 1:length(data.seqs)
    seqsOK(i) = ~isempty(regexp(data.seqs{i}, '-\.[ACDEFGHIKLMNPQRSTVWY]{6}\.-', 'once'));
end
fprintf('%d / %d data points have valid sequences\n', length(find(seqsOK)), length(data.seqs));

% process sequences
seqs = data.seqs(seqsOK);
procSeqs = seqs;
for i = 1:length(procSeqs)
    procSeqs{i} = procSeqs{i}(3:end-2);
end

% go over each repeat
Kd = nan(length(find(seqsOK)), 1);
KdError = 0*Kd + inf; logKdError = KdError;
a = Kd; xcor = Kd;
KdHi = Kd; KdLo = Kd;
if (~isfield(inputs, 'show')), inputs.show = 1; end
x = data.ctrIn(seqsOK);
y = data.ctrOut(seqsOK);
[g, ~] = gaussianErrorModelWithBaseline(x, y);
b = g(end-1:end);

X = data.expOut(seqsOK);
Y = data.ctrOut(seqsOK);
xc = data.xcorr(seqsOK);
valid = find(~isnan(X) & ~isnan(Y)); X = X(valid); Y = Y(valid);
xc = xc(valid);

if (inputs.show)
    figure;
    subplot(1, 3, 1);
    hold off; plot(y, x*b(1) + b(2), 'k.'); hold on;
    plot(sort(y), sort(y) + 3*sqrt(myVar(sort(y), g))); 
    plot(sort(y), sort(y) - 3*sqrt(myVar(sort(y), g)));

    subplot(1, 3, 2);
    plot(sort(y), sqrt(myVar(sort(y), g))); 
end

% compute ratio alpha = [peptide outside in experiment]/[peptide outside in control]
[alpha, alphaStd] = ratioWithErrorPropagation(X, Y, sqrt(myVar(X, g)), sqrt(myVar(Y, g)));
okPoints = find(alphaStd < 0.4);
fprintf('--> %d / %d points have tolerable error\n', length(okPoints), length(alpha));
if (inputs.show)
    subplot(1, 3, 3);
    errorbar(1:length(X), alpha, alphaStd, alphaStd, 'k.');
end

% compute Kd range
[Kdr, error, errorLog, KdRange] = estKd(inputs.domainConc, alpha(okPoints), alphaStd(okPoints));
KdError(valid(okPoints)) = error;
logKdError(valid(okPoints)) = errorLog;
Kd(valid(okPoints)) = Kdr;
a(valid(okPoints)) = alpha(okPoints);
xcor(valid(okPoints)) = xc(okPoints);
KdLo(valid(okPoints)) = KdRange(:, 1);
KdHi(valid(okPoints)) = KdRange(:, 2);
result = struct('Kd', Kd, 'KdError', KdError, 'logKdError', logKdError, 'KdLo', KdLo, 'KdHi', KdHi, 'seqs', {procSeqs}, "a", a, "xcor", xcor);


% --- Gaussian error modeling --- %
% figures out a Gaussian error model of the form
% yi = xi + err(xi) = xi + A*gauss(0, B*xi + C)
function [g, o] = gaussianErrorModel(x0, y0)

ii = find(~isnan(x0) & ~isnan(y0));
x = [x0(ii); y0(ii)];
y = [y0(ii); x0(ii)];
% err = y - x;
% p = exp(-(err.^2)./(B*x + C));

opts = optimset('Display', 'off');
g = fminunc(@(p) myObj(x, y, p), [0.2 0 100], opts);
g = fminsearch(@(p) myObj(x, y, p), g, opts);
o = myObj(x, y, g);

% figures out a Gaussian error model of the form
% yi = a*xi + b + err(xi) = a*xi + b +  + A*gauss(0, B*xi + C)
function [g, o] = gaussianErrorModelWithBaseline(x0, y0)

ii = find(~isnan(x0) & ~isnan(y0));
x = x0(ii);
y = y0(ii);

% first do a linear regression
b = regress(y, [x (0*x + 1)]);

opts = optimset('Display', 'off');
g = fminunc(@(p) myObj1(x, y, p), [0.2 0 0.1 100 b'], opts);
g = fminsearch(@(p) myObj1(x, y, p), g, opts);
o = myObj1(x, y, g);

function o = myObj1(x, y, p)
b = p(end-1:end);
yp = [x (0*x + 1)]*b';
o = -mean(gauss(yp - y, 0, myVar(y, p)));

function o = myObj(x, y, p)
o = -mean(gauss(x - y, 0, myVar(x, p)));

function v = gauss(x, mu, s)
% v = log((1./(s*sqrt(2*pi))).*exp(-0.5*((x - mu)./sqrt(s)).^2));
v = log(1./(s*sqrt(2*pi))) - 0.5*((x - mu)./sqrt(s)).^2;

function v = myVar(x, p)
v = abs(abs(p(1))*x.*x + p(2)*x + p(3)*sqrt(abs(x)) + abs(p(4)));
% v = abs(p(1)*x.*x + p(2)*x + abs(p(3)));
% v = p(1)*x + p(2);

function [r, rStd] = ratioWithErrorPropagation(X, Y, stdX, stdY)
% from D.T. Holmes, K.A. Buhr / Clinical Biochemistry 40 (2007) 728?734
cvy = stdY./Y;
cvx = stdX./X;
rStd = (X./Y).*sqrt((cvx.^2) + (cvy.^2) + 3*(cvy.^2).*(cvx.^2) + 8*(cvy.^4));
r = (X./Y).*(1 + cvy.^2);

function [r, rStd] = ratioWithErrorPropagationWrong(X, Y, stdX, stdY)
% from D.T. Holmes, K.A. Buhr / Clinical Biochemistry 40 (2007) 728?734
cvy = stdY./Y;
cvx = stdX./X;
rStd = (X./Y).*sqrt((cvx.^2) + (cvy.^2) + 3*(cvy.^2).*(cvx.^2) + 8*(cvy.^2));
r = (X./Y).*(1 + cvy.^2);

function [Kd, errorKd, errorLogKd, KdRange] = estKd(D, alpha, alphaStd)
errorKd = nan(length(alpha), 1); % lower and upper bounds of Kd
errorLogKd = errorKd;
Kd = nan(length(alpha), 1);
KdRange = -1*ones(length(alpha), 2); % lower and upper bounds of Kd
alphaLo = alpha - alphaStd;
alphaHi = alpha + alphaStd;

% estimated alpha in normal range
ii = find((alpha > 0) & (alpha < 1));
Kd(ii, 1) = D*alpha(ii)./(1 - alpha(ii));
errorKd(ii) = D*alphaStd(ii)./((1-alpha(ii)).^2);
errorLogKd(ii) = (1./alpha(ii) + 1./(1 - alpha(ii))) .* alphaStd(ii);

% estimated alpha is negative
ii = find(alpha <= 0);
Kd(ii) = 0;
errorKd(ii) = inf;
errorLogKd(ii) = inf;

% estimated alpha is over 1
ii = find(alpha >= 1);
Kd(ii) = inf;
errorKd(ii) = inf;
errorLogKd(ii) = inf;

ii = find((alphaLo > 0) & (alphaLo < 1)); KdRange(ii, 1) = D*alphaLo(ii)./(1 - alphaLo(ii));
ii = find(alphaLo <= 0); KdRange(ii, 1) = 0;
ii = find(alphaLo >= 1); KdRange(ii, 1) = inf;
ii = find((alphaHi > 0) & (alphaHi < 1)); KdRange(ii, 2) = D*alphaHi(ii)./(1 - alphaHi(ii));
ii = find(alphaHi <= 0); KdRange(ii, 2) = 0;
ii = find(alphaHi >= 1); KdRange(ii, 2) = inf;

function [mm, alpha] = modelMatrix(seqs)
mm = [];
if (isempty(seqs)), return; end
L = length(seqs{1});

% build positional alphabet
alpha = cell(1, L);
for i = 1:length(seqs)
    for j = 1:L
        alpha{j} = [alpha{j}, seqs{i}(j)];
    end
end
N = 1;
for j = 1:L
    alpha{j} = unique(alpha{j});
    N = N + length(alpha{j}) - 1;
end

% build model matrix
mm = zeros(length(seqs), N);
mm(:, 1) = 1;
for i = 1:length(seqs)
    off = 0;
    for j = 1:L
        k = find(alpha{j} == seqs{i}(j));
        if (k > 1), mm(i, off + k) = mm(i, off + k) + 1; end
        off = off + length(alpha{j}) - 1;
    end
end
