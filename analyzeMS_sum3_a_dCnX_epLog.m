%this version is for error propagation
%this version if for trial3 with only one set of control, and summary three
%repeats together
%print alpha and dcorr values
%this version has confidence intervals listed, used fitlm instead of
%regress
function analyzeMS_sum3_a_dCnX_epLog(outBase, dataFileXLS)
if (~exist('dataFileXLS', 'var'))
    dataFileXLS = '80431_33TMTintensities.xlsx';
end

D = 35.29;        % overall domain concentration in experiment (100uM*300ul/850ul)

units = 10^-6; % concentration units in M
% change how 0 to 1 to make two figures
res1 = fitKds(struct('xlsFile', dataFileXLS, 'sheetName', 'Sheet1','region','B1:L2247', 'domainConc', D, 'show', 1));
res2 = fitKds(struct('xlsFile', dataFileXLS, 'sheetName', 'Sheet1','region','S1:AC2218', 'domainConc', D, 'show', 1));
res3 = fitKds(struct('xlsFile', dataFileXLS, 'sheetName', 'Sheet1','region','AJ1:AT2201', 'domainConc', D, 'show', 1));

% now learning sequence-Kd mapping
M = containers.Map; % mean Kd for each unique sequence
A = containers.Map; % all Kds for each unique sequence
E = containers.Map; % error (size of CI interval) of each Kd estimate for each unique sequence
F = containers.Map; % all alpha values for each unique peptide
G = containers.Map; % mean alpha values for each unique peptide
H = containers.Map; % all dcorr values
I = containers.Map; % mean dcorr values
J = containers.Map; % all Xcorr values
K = containers.Map; % mean Xcorr values
allResult = [res1,res2,res3];

for j = 1:3
    for i = 1:size(allResult(j).Kd, 1) %size of allResult(1) is 1, size of AllResult(1).Kd is three columns
        idx = isfinite(allResult(j).Kd(i, :)); % Kd has multiple columns here,so this for loop will go through every row in every column
        kds = log10(allResult(j).Kd(i, idx));
        kd_ints = allResult(j).KdError(i, idx);
        f = allResult(j).a(i,idx);
        dCn = allResult(j).dcor(i,idx);
        xcorr = allResult(j).xcor(i,idx);
        %if (~isempty(kds)) %uncomment line60 if uncomment this line
        if all(~isempty(kds))&& all(~isinf(kd_ints))
        
        %if ((~isempty(kds))&(~isinf(kd_ints)))% this also works
        % '&' for 3 repeats(every element in array), '&&' for one repeat   
        %if all(~isempty(kds)&& all(~isinf(kd_ints))&& all(kds >= 4.4 & kds <= 160))
        %kds >= 4.4 & kds <= 160 check the alpha values are between 0.1 and
        %0.8 for D=40uM
       
            if (isKey(M, allResult(j).seqs{i}))
                A(allResult(j).seqs{i}) = [A(allResult(j).seqs{i}) kds];
                M(allResult(j).seqs{i}) = mean(A(allResult(j).seqs{i}));
                E(allResult(j).seqs{i}) = [E(allResult(j).seqs{i}) kd_ints];
                F(allResult(j).seqs{i}) = [F(allResult(j).seqs{i}) f];
                G(allResult(j).seqs{i}) = mean(F(allResult(j).seqs{i}));
                H(allResult(j).seqs{i}) = [H(allResult(j).seqs{i}) dCn];
                I(allResult(j).seqs{i}) = mean(H(allResult(j).seqs{i}));
                J(allResult(j).seqs{i}) = [J(allResult(j).seqs{i}) xcorr];
                K(allResult(j).seqs{i}) = mean(J(allResult(j).seqs{i}));
            else
                M(allResult(j).seqs{i}) = mean(kds);
                A(allResult(j).seqs{i}) = kds;
                E(allResult(j).seqs{i}) = kd_ints;
                F(allResult(j).seqs{i}) = f;
                G(allResult(j).seqs{i}) = mean(f);
                H(allResult(j).seqs{i}) = dCn;
                I(allResult(j).seqs{i}) = mean(dCn);
                J(allResult(j).seqs{i}) = xcorr;
                K(allResult(j).seqs{i}) = mean(xcorr); 
            end
        end
    end
    fprintf('%d unique peptides \n', length(keys(M)));
end


% write summary
fid = fopen(sprintf('%s.csv', outBase), 'w');
fprintf(fid, 'sequence, log(Kd_estimate), standard_error(logKd), error (est), number of samples, alpha, dcorr, dcorr error,xcorr\n');
keySeqs = keys(M);
vals = values(M);
vals = [vals{:}];
alphaValue = values(G);
alphaValue = [alphaValue{:}];
dVal = values(I);
dVal = [dVal{:}];
xVal = values(K);
xVal = [xVal{:}];
for i = 1:length(keySeqs)
    err_est = mean(E(keySeqs{i}))/sqrt(length(A(keySeqs{i})));
    %err_est = mean(E(keySeqs{i}));
    %if isinf(err_est), err_est = 1000; end
    fprintf(fid, '%s,%f,%f,%f,%f,%f,%f,%f,%f\n', keySeqs{i}, vals(i), std(A(keySeqs{i}))/sqrt(length(A(keySeqs{i}))), err_est, length(A(keySeqs{i})),alphaValue(i),dVal(i),std(H(keySeqs{i}))/sqrt(length(H(keySeqs{i}))),xVal(i));
end
fclose(fid);

% file a site-independent model
[mm, alpha] = modelMatrix(keySeqs);
%b = regress(log10(vals') + log10(units), mm);
%b = lscov(mm,log10(vals') + log10(units),1./(vals.^2));
%b = lscov(mm,log10(vals') + log10(units),1./vals);
b = fitlm(mm, log10(vals') + log10(units));
ci = coefCI(b,.01); % this confidence interval only works for fitlm
%ci returns two columns, mean-error, mean+error
m=ci([1 3:end],:); % 2nd row is 0, does not exist in regress, can't match mm
%turn the two columns to mean and error
m1 = m(:,1);
m2 = m(:,2);
m_hat = (m1+m2)/2;
m_error = (m2-m1)/2;
figure;%add this line to make two figures
plot(log10(vals) + log10(units), mm*m_hat, 'o');
set(gca, 'FontSize', 14);
xlabel('experimental log10(Kd [M])');
ylabel('linear-model prediction');

% print optimal parameters
k = 1;
fprintf('--- parameters ---\n');
fprintf('const = %f\t%f\n', m_hat(k),m_error(k));% was b(k)
for i = 1:length(alpha)
    fprintf('--> site %d\n', i);
    for j = 1:length(alpha{i})
        if (j == 1)
            fprintf('\t\t%s\t0.0\n', alpha{i}(j));
        else
            fprintf('\t\t%s\t%.2f\t%.2f\n', alpha{i}(j), m_hat(k+1),m_error(k+1));
            k = k + 1;
        end
    end
end

function ret = fitKds(inputs)
% read data
[num, ~, raw] = xlsread(inputs.xlsFile, inputs.sheetName, inputs.region);
allSeqs = raw(2:end, 7);
 
% collect data from the three repeats
xcorr = num(:, 4);
dcorr = num(:, 8);% put three repeats in three columns
expOut = num(:, 9);
ctrIn = num(:, 11);
ctrOut = num(:, 10);

data = struct('seqs', {allSeqs},'xcorr',xcorr, 'dcorr', dcorr, 'expOut', expOut, 'ctrIn', ctrIn, 'ctrOut', ctrOut);
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
Kd = nan(length(find(seqsOK)), 1);
KdError=Kd;
a=Kd;dcor=Kd;xcor=Kd;
if (~exist('show', 'var')), show = 1; end
x = data.ctrIn(seqsOK);
y = data.ctrOut(seqsOK);
[g, ~] = gaussianErrorModelWithBaseline(x, y);
b = g(end-1:end);

X = data.expOut(seqsOK);
Y = data.ctrOut(seqsOK);
d = data.dcorr(seqsOK);
xc = data.xcorr(seqsOK);
valid = find(~isnan(X) & ~isnan(Y)); X = X(valid); Y = Y(valid);
d = d(valid);
xc = xc(valid);
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
[Kdr, error] = estKd(inputs.domainConc, alpha(okPoints), alphaStd(okPoints));
KdError(valid(okPoints)) = error;
Kd(valid(okPoints)) = Kdr;
a(valid(okPoints)) = alpha(okPoints);
dcor(valid(okPoints)) = d(okPoints);
xcor(valid(okPoints)) = xc(okPoints);
result = struct('Kd', Kd, 'KdError', KdError, 'seqs', {procSeqs},"a", a,"dcor",dcor,"xcor",xcor);



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

function [Kd, errorKd] = estKd(D, alpha, alphaStd)
errorKd = -1*ones(length(alpha), 1); % lower and upper bounds of Kd
Kd = -1*ones(length(alpha), 1);

ii = find((alpha > 0) & (alpha < 1)); Kd(ii, 1) = D*alpha(ii)./(1 - alpha(ii));errorKd(ii)=alphaStd(ii)./((alpha(ii)-alpha(ii).^2)*log(10));
ii = find(alpha <= 0); Kd(ii) = 0;errorKd(ii)=inf;
ii = find(alpha >= 1); Kd(ii) = inf;errorKd(ii)=inf;

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
