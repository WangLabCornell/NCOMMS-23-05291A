%% Aligns all traces from repeated winding experiment for pre-buckled activity
% Run alignOneTrace.m first for all traces to generate .mat files for each trace

files = dir("*.mat");

figure()
hold on
allBinnedExts = [];
bucklings = [];
centers = [];
bucklingZs = [];
rfits = [];
binningTimes = [];
binningExts = [];
allhatturns = [];
allhatZs = [];
tbins = [-4.95:.1:14.95];
for i = 1:numel(files)
    load(files(i).name)
    centers = [centers rfit1(2)];
    rfit1(2) = 0;
    for j = 1:length(tbins)
        binnedExts(j) = mean(allExts(allTimes >= tbins(j)-.05 & allTimes < tbins(j)+.05), 'omitnan') + f_3piece(rfit1, rfit1(4));
    end
    binningTimes = [binningTimes allTimes];
    binningExts = [binningExts; allExts];
    allBinnedExts = [allBinnedExts; binnedExts];
    bucklings = [bucklings rfit1(4)];
    bucklingZs = [bucklingZs f_3piece(rfit1, rfit1(4))];
    allhatturns = [allhatturns; lowTurns-centers(end)];
    allhatZs = [allhatZs; lowhat];
    rfits = [rfits; rfit1];
    endExts(i) = mean(allExts(allTimes >= 10 & allTimes < 15), 'omitnan') + f_3piece(rfit1, rfit1(4));
end
for i = 1:numel(files)
    allBinnedExts(i, :) = allBinnedExts(i, :) - (endExts(i) - mean(endExts));
    plot(tbins, allBinnedExts(i, :))
end
ylim([2200 4000])
plot(tbins, mean(allBinnedExts, 1, 'omitnan'), 'k', 'LineWidth', 2)
plot([-5 15], mean(bucklingZs)*[1 1], 'k--')
plot([0 0], [2200 4000], 'k--')
meanBinnedExts = mean(allBinnedExts, 'omitnan');
xlabel('Time since buckling (s)', 'FontSize', 14)
ylabel('Extension (nm)', 'FontSize', 14)
rfittotal = lsqcurvefit(@f_3piece, [-.5 0 3200 20 50], allhatturns(:), allhatZs(:), [-5 0], [0 0]);
figure()
for i = 1:1
    meanBinnedTurns = interp1(f_3piece(rfittotal, 0:.001:60), 0:.001:60, smooth(meanBinnedExts, 10*i));
    plot(tbins, meanBinnedTurns)
    hold on
end
leg = legend('1 s smoothing', '2 s smoothing', '3 s smoothing', '4 s smoothing', 'Buckling');
leg.Location = 'southwest';
leg.FontSize = 12;
xlabel('Turn state', 'FontSize', 14)
ylabel('Relaxation rate (turns/s)', 'FontSize', 14)