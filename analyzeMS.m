%this version if for trial6 with only one set of control
function analyzeMS(outBase, dataFileXLS, setName, doCorrection)
if (~exist('dataFileXLS', 'var') || isempty(dataFileXLS))
    dataFileXLS = 'data/86252_55_86900_01TMT_yu.xlsx';
end
if (~exist('doCorrection', 'var'))
    doCorrection = 0;
end

D = 40;        % overall domain concentration in experiment (100uM*500ul/1250ul)
units = 10^-6; % concentration units in M
% change how 0 to 1 to make two figures
res = fitKds(struct('xlsFile', dataFileXLS, 'sheetName', setName, 'domainConc', D, 'show', 1, 'doCorrection', doCorrection));

% now learning sequence-Kd mapping
M = containers.Map; % mean Kd for each unique sequence
A = containers.Map; % all Kds for each unique sequence
E = containers.Map; % error (size of CI interval) of each Kd estimate for each unique sequence
for i = 1:size(res.Kd, 1)
    idx = isfinite(res.Kd(i, :));
    kds = res.Kd(i, idx);
    kd_ints = res.KdHi(i, idx) - res.KdLo(i, idx);
    if (~isempty(kds))
        if (isKey(M, res.seqs{i}))
            A(res.seqs{i}) = [A(res.seqs{i}) kds];
            M(res.seqs{i}) = mean(A(res.seqs{i}));
            E(res.seqs{i}) = [E(res.seqs{i}) kd_ints];
        else
            M(res.seqs{i}) = mean(kds);
            A(res.seqs{i}) = kds;
            E(res.seqs{i}) = kd_ints;
        end
    end
end
fprintf('%d unique peptides\n', length(keys(M)));

% write summary
fid = fopen(sprintf('%s.csv', outBase), 'w');
fprintf(fid, 'sequence, Kd estimate, Kd standard error (obs), error (est), number of samples\n');
keySeqs = keys(M);
vals = values(M); vals = [vals{:}];
for i = 1:length(keySeqs)
    err_est = mean(E(keySeqs{i}))/sqrt(length(A(keySeqs{i})));
    if isinf(err_est), err_est = 1000; end
    fprintf(fid, '%s,%f,%f,%f,%f\n', keySeqs{i}, vals(i), std(A(keySeqs{i}))/sqrt(length(A(keySeqs{i}))), err_est, length(A(keySeqs{i})));
end
fclose(fid);

% file a site-independent model
[mm, alpha] = modelMatrix(keySeqs);
b = regress(log10(vals') + log10(units), mm);
figure;%add this line to make two figures
plot(log10(vals) + log10(units), mm*b, 'o');
set(gca, 'FontSize', 14);
xlabel('experimental log10(Kd [M])');
ylabel('linear-model prediction');

% print optimal parameters
k = 1;
fprintf('--- parameters ---\n');
fprintf('const = %f\n', b(k));
for i = 1:length(alpha)
    fprintf('--> site %d\n', i);
    for j = 1:length(alpha{i})
        if (j == 1)
            fprintf('\t\t%s\t0.0\n', alpha{i}(j));
        else
            fprintf('\t\t%s\t%.2f\n', alpha{i}(j), b(k+1));
            k = k + 1;
        end
    end
end

function ret = fitKds(inputs)
% read data
[num, ~, raw] = xlsread(inputs.xlsFile, inputs.sheetName);
allSeqs = raw(2:end, 8);

% collect data from the three repeats
dcorr = num(:, 5);
xcorr = num(:, 4);
expOut = num(:, 9);
ctrIn = num(:, 11);
ctrOut = num(:, 10);

data = struct('seqs', {allSeqs}, 'dcorr', dcorr, 'xcorr', xcorr, 'expOut', expOut, 'ctrIn', ctrIn, 'ctrOut', ctrOut);
ret = fitKdsFromData(inputs, data);


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
KdLo = nan(length(find(seqsOK)), 1);
KdHi = KdLo; Kd = KdLo;
if (~exist('show', 'var')), show = 1; end
x = data.ctrIn(seqsOK);
y = data.ctrOut(seqsOK);
[g, ~] = gaussianErrorModelWithBaseline(x, y);
b = g(end-1:end);

X = data.expOut(seqsOK);
Y = data.ctrOut(seqsOK);
valid = find(~isnan(X) & ~isnan(Y)); X = X(valid); Y = Y(valid);

if inputs.doCorrection
    [Yc, p] = correctByMedianFit(X, Y);
    plot(Y, X, '.');
    hold on;
    plot([min(Y) max(Y)], p(1) + p(2)*[min(Y) max(Y)], '-');
    xlabel('Control Out');
    ylabel('Experiment Out');
    title('Signal correction');
    legend('Data', 'Correction line');
    Y = Yc;
end

if (show)
    figure;
    subplot(1, 3, 1);% first column of figures in figure1
%     hold off; plot(x, y, 'k.'); hold on;
    hold off; plot(y, x*b(1) + b(2), 'k.'); hold on;
    plot(sort(y), sort(y) + 3*sqrt(myVar(sort(y), g))); % first column of figures in figure1
    plot(sort(y), sort(y) - 3*sqrt(myVar(sort(y), g)));
%       plot(sort(x), (sort(x) + g)); %why not these???
%       plot(sort(x), (sort(x) - g));
        %set(gca, 'YScale', 'lin', 'xscale', 'log') %comment this line out
        %to give log scale in figure1 most left panels
    subplot(1, 3, 2);% second column of figures in figure1, x-axis is x, y-axis 
    plot(sort(y), sqrt(myVar(sort(y), g))); 
end

% compute ratio alpha = [peptide outside in experiment]/[peptide outside in control]
[alpha, alphaStd] = ratioWithErrorPropagation(X, Y, sqrt(myVar(X, g)), sqrt(myVar(Y, g)));
okPoints = find(alphaStd < 0.4);
fprintf('--> %d / %d points have tolerable error\n', length(okPoints), length(alpha));
if (show)
    subplot(1, 3, 3);% third column of figure1
    errorbar(1:length(X), alpha, alphaStd, alphaStd, '.');
end

% compute Kd range
[Kdr, KdRange] = estimateKd(inputs.domainConc, alpha(okPoints), alphaStd(okPoints));
KdLo(valid(okPoints)) = KdRange(:, 1);
KdHi(valid(okPoints)) = KdRange(:, 2);
Kd(valid(okPoints)) = Kdr;

result = struct('Kd', Kd, 'KdLo', KdLo, 'KdHi', KdHi, 'seqs', {procSeqs});


function [Yc, p] = correctByMedianFit(X, Y)
opts = optimset('Display', 'off');
p = fminunc(@(p) median(abs(p(1) + Y*p(2) - X)), [0 1], opts);
p = fminsearch(@(p) median(abs(p(1) + Y*p(2) - X)), p, opts);
Yc = p(1) + Y*p(2);


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
v = abs(p(1)*x.*x + p(2)*x + p(3)*sqrt(abs(x)) + abs(p(4)));
% v = abs(p(1)*x.*x + p(2)*x + abs(p(3)));
% v = p(1)*x + p(2);

function [r, rStd] = ratioWithErrorPropagation(X, Y, stdX, stdY)
% from D.T. Holmes, K.A. Buhr / Clinical Biochemistry 40 (2007) 728?734
cvy = stdY./Y;
cvx = stdX./X;
rStd = (X./Y).*sqrt((cvx.^2) + (cvy.^2) + 3*(cvy.^2).*(cvx.^2) + 8*(cvy.^2));
r = (X./Y).*(1 + cvy.^2);

function [Kd, KdRange] = estimateKd(D, alpha, alphaStd)
KdRange = -1*ones(length(alpha), 2); % lower and upper bounds of Kd
Kd = -1*ones(length(alpha), 1);
alphaLo = alpha - alphaStd;
alphaHi = alpha + alphaStd;

ii = find((alphaLo > 0) & (alphaLo < 1)); KdRange(ii, 1) = D*alphaLo(ii)./(1 - alphaLo(ii));
ii = find(alphaLo <= 0); KdRange(ii, 1) = 0;
ii = find(alphaLo >= 1); KdRange(ii, 1) = inf;
ii = find((alphaHi > 0) & (alphaHi < 1)); KdRange(ii, 2) = D*alphaHi(ii)./(1 - alphaHi(ii));
ii = find(alphaHi <= 0); KdRange(ii, 2) = 0;
ii = find(alphaHi >= 1); KdRange(ii, 2) = inf;
ii = find((alpha > 0) & (alpha < 1)); Kd(ii, 1) = D*alpha(ii)./(1 - alpha(ii));
ii = find(alpha <= 0); Kd(ii) = 0;
ii = find(alpha >= 1); Kd(ii) = inf;

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
