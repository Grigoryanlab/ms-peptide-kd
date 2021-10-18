function plotFpMsKd(pedal)
if (~exist('pedal', 'var'))
    pedal = 'trial3-ep2.csv';
end
fp_result = 'Kd_FP_sum.xlsx';
Range = 'A1:C34';
[num, ~, raw] = xlsread(fp_result,Range);
seq = raw(2:end,1);
x=num(:, 1);%fp
xE=num(:, 2);% fp errors
y=nan(length(x),1);
yE=nan(length(y),1);
T = readtable(pedal);
seq3 = T.sequence;
kd = T.Kd_estimate; % pedal kd
est = T.error_est_; % pedal estimated error
for i = 1:length(seq)
    for j = 1: length(seq3)
        if (seq{i}==seq3{j})
            y(i)=kd(j);
            yE(i)=est(j);
        end
    end
end
x=x./1000000;
xE=xE./1000000;
y=y./1000000;
yE=yE./1000000;
figure;
errorbar(log10(x), log10(y), yE./(y*log(10)), yE./(y*log(10)), xE./(x*log(10)),xE./(x*log(10)), 'o','MarkerSize',5,'MarkerEdgeColor','black','MarkerFaceColor','black','LineStyle','none', 'Color', [0.5,0.5,0.5],'linewidth', 0.5)

xlabel('FP log10(Ki or Kd/[M])')
ylabel('PEDAL log10(Kd/[M])')
% xlim([-7.5,-3])
% ylim([-7,-3])
%title('trial3(4096_pep_lib)')



