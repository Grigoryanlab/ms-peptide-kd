%this version if for trial6 with only one set of control
function analyzeMS_yu4(outBase)
dataFileXLS = '86252_55_86900_01TMT_yu.xlsx';
D = 40;     % overall domain concentration in experiment (100uM*500ul/1250ul)
units = 10^-6; % concentration units in M
[Kd, KdLo, KdHi, seqs] = fitKds(dataFileXLS, D, 1); %change 0 to 1 to make two figures

% now learning sequence-Kd mapping
M = containers.Map;
A = containers.Map;
for i = 1:size(Kd, 1)
    kds = Kd(i, isfinite(Kd(i, :)));
    if (~isempty(kds))
        if (isKey(M, seqs{i}))
            A(seqs{i}) = [A(seqs{i}) kds];
            M(seqs{i}) = mean(A(seqs{i}));
        else
            M(seqs{i}) = mean(kds);
            A(seqs{i}) = kds;
        end
    end
end
fprintf('%d unique peptides\n', length(keys(M)));

% write summary
fid = fopen(sprintf('%s.csv', outBase), 'w');
fprintf(fid, 'sequence, Kd estimate, Kd standard error, number of samples\n');
keySeqs = keys(M);
vals = values(M); vals = [vals{:}];
for i = 1:length(keySeqs)
    fprintf(fid, '%s,%f,%f,%f\n', keySeqs{i}, vals(i), std(A(keySeqs{i}))/sqrt(length(A(keySeqs{i}))), length(A(keySeqs{i})));
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

function [Kd, KdLo, KdHi, procSeqs] = fitKds(dataFileXLS, D, show)
% read data
[num, ~, raw] = xlsread(dataFileXLS, 'Set2', 'B2:L319');
allSeqs = raw(1:end, 7);

% collect data from the three repeats
repeats = cell(1, 1);
repeats{1}.dcorr = num(:, 5);
repeats{1}.xcorr = num(:, 4);
repeats{1}.expOut = num(:, 9);
repeats{1}.ctrIn = num(:, 11);
repeats{1}.ctrOut = num(:, 10);


% mark which sequences are okay
seqsOK = false(1, length(allSeqs));
for i = 1:length(allSeqs)
    seqsOK(i) = ~isempty(regexp(allSeqs{i}, '-\.[ACDEFGHIKLMNPQRSTVWY]{6}\.-', 'once'));
end
fprintf('%d / %d data points have valid sequences\n', length(find(seqsOK)), length(allSeqs));

% process sequences
seqs = allSeqs(seqsOK);
procSeqs = seqs;
for i = 1:length(procSeqs)
    procSeqs{i} = procSeqs{i}(3:end-2);
end

% go over each repeat
KdLo = nan(length(find(seqsOK)), length(repeats{1}));
KdHi = KdLo; Kd = KdLo;
if (~exist('show', 'var')), show = 1; end
x = repeats{1}.ctrIn(seqsOK);
y = repeats{1}.ctrOut(seqsOK);
[g, ~] = gaussianErrorModel(x, y);

if (show)
    subplot(length(repeats{1}), 3, 1);% first column of figures in figure1
    hold off; plot(x, y, 'k.'); hold on;
    plot(sort(x), sort(x) + 3*sqrt(myVar(sort(x),g)));% first column of figures in figure1
    plot(sort(x), sort(x) - 3*sqrt(myVar(sort(x),g)));
%       plot(sort(x), (sort(x) + g)); %why not these???
%       plot(sort(x), (sort(x) - g));
        %set(gca, 'YScale', 'lin', 'xscale', 'log') %comment this line out
        %to give log scale in figure1 most left panels
    subplot(length(repeats{1}), 3, 2);% second column of figures in figure1, x-axis is x, y-axis 
    plot(sort(x), sqrt(myVar(sort(x), g))); 
end

X = repeats{1}.expOut(seqsOK);
Y = repeats{1}.ctrOut(seqsOK);
valid = find(~isnan(X) & ~isnan(Y)); X = X(valid); Y = Y(valid);

% compute ratio alpha = [peptide outside in experiment]/[peptide outside in control]
[alpha, alphaStd] = ratioWithErrorPropagation(X, Y, sqrt(myVar(Y, g)), sqrt(myVar(X, g)));
okPoints = find(alphaStd < 0.4);
fprintf('--> %d / %d points have tolerable error\n', length(okPoints), length(alpha));
if (show)
    subplot(length(repeats{1}), 3, 3);% third column of figure1
    errorbar(1:length(X), alpha, alphaStd, alphaStd, '.');
end

% compute Kd range
[Kdr, KdRange] = estimateKd(D, alpha(okPoints), alphaStd(okPoints));
KdLo(valid(okPoints)) = KdRange(:, 1);
KdHi(valid(okPoints)) = KdRange(:, 2);
Kd(valid(okPoints)) = Kdr;




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

function o = myObj(x, y, p)
o = -mean(gauss(x - y, 0, myVar(x, p)));

function v = gauss(x, mu, s)
% v = log((1./(s*sqrt(2*pi))).*exp(-0.5*((x - mu)./sqrt(s)).^2));
v = log(1./(s*sqrt(2*pi))) - 0.5*((x - mu)./sqrt(s)).^2;

function v = myVar(x, p)
v = abs(p(1)*x.*x + p(2)*x + abs(p(3)));
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
