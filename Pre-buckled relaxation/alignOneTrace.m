%% Aligns relaxation events from individual traces
% Needs the following files: 
% 1. CircleFitByPratt (https://www.mathworks.com/matlabcentral/fileexchange/22643-circle-fit-pratt-method)
% 2. f_3piece: 3-piece function to fit to naked DNA hat curve

tetherIndex = 1+1;
stuckBeadIndex = 0+1;
FPS = 40;
 
%% Read data

files = dir('*.txt');
for i = 1:numel(files)
    tempData = importdata(files(i).name);
    importedData = tempData.data(:,8:end);
    numTethers(i) = size(importedData, 2) / 3;
    for j = 1:numTethers(i)
        Xs{i, j} = importedData(:, (j-1)*3 + 1)*146.5*2;
        Ys{i, j} = importedData(:, (j-1)*3 + 2)*146.5*2;
        Zs{i, j} = importedData(:, (j-1)*3 + 3)*1000*1.33*.5;
    end
    turns{i} = tempData.data(:, 3);
    taskIndex{i} = tempData.data(:, 2);
end

%% Geometry correction

geomIndex = 2;

for j = 1:numTethers(1)
    for k = 1:8
        tracegX = Xs{geomIndex, j};
        tracegY = Ys{geomIndex, j};
        Rxs(k) = mean(tracegX(turns{geomIndex} > (k-1) * 0.125-.02 & turns{geomIndex} < (k-1) * 0.125+.02), "omitnan");
        Rys(k) = mean(tracegY(turns{geomIndex} > (k-1) * 0.125-.02 & turns{geomIndex} < (k-1) * 0.125+.02), "omitnan");
    end
    circleFit = CircleFitByPratt([(Rxs-mean(Rxs))' (Rys-mean(Rys))']);
    Rpara = max([min([circleFit(3) 500]) 0]);
    Rperp = sqrt(500^2 - Rpara^2);
    geomCorr(j) = 500 - Rperp;
    magnetAngle(j, 1) = Rxs(1)-mean(tracegX);
    magnetAngle(j, 2) = Rys(1)-mean(tracegY);
end

figure()
plot(magnetAngle(:, 1), magnetAngle(:, 2), '.')
angleFit = polyfit(magnetAngle(:, 2), magnetAngle(:, 1), 1);
hold on
plot(angleFit(1) * [-500:500] + angleFit(2), [-500:500], 'k--')
theta = pi/2 - atan(angleFit(1));
title(['\theta = ' num2str(theta * 180/pi) ' degrees'])

%% Wind to surface correction

surfIndex = 3;

for j = tetherIndex
    windCorr(j) = mean(Zs{surfIndex, j}(end-39:end)) - mean(Zs{surfIndex, stuckBeadIndex}(end-39:end));
    extension(j) = mean(Zs{surfIndex, j}(turns{surfIndex} == 0)) - mean(Zs{surfIndex, j}(turns{surfIndex} > max(turns{surfIndex}) - 0.02));
end

%% Process hats

lowhatIndex = 4;

for j = tetherIndex
    hatTurns = round(min(turns{lowhatIndex})):round(max(turns{lowhatIndex}));
    for k = 1:numel(hatTurns)
        lowhat(k) = mean(Zs{lowhatIndex, j}(abs(turns{lowhatIndex} - hatTurns(k)) < 0.02)) - mean(Zs{lowhatIndex, stuckBeadIndex}(abs(turns{lowhatIndex} - hatTurns(k)) < 0.02)) - windCorr(j) + geomCorr(j);
    end
    lowTurns = hatTurns;
end

%% Fit hats

[maxZ, maxTurnIndex] = max(lowhat);
r0 = [-.5 hatTurns(maxTurnIndex) maxZ 20 50];
try
    rfit1 = lsqcurvefit(@f_3piece, r0, lowTurns(lowhat > 1000), lowhat(lowhat > 1000));
catch
    rfit1 = [nan nan nan nan nan];
end


%% Force calibrations

forceIndex = 5;

for j = tetherIndex
    extension = mean(Zs{forceIndex, j}) - mean(Zs{forceIndex, stuckBeadIndex}) + geomCorr(j) - windCorr(j);
    Xpara = Xs{forceIndex, j} * cos(theta) + Ys{forceIndex, j} * sin(theta);
    Xperp = Xs{forceIndex, j} * -sin(theta) + Ys{forceIndex, j} * cos(theta);
    FparaLow = 4.114 / var(Xpara) * extension;
    FperpLow = 4.114 / var(Xperp) * (extension + 500);
end

%% Plot trace

clampedExt = Zs{1, tetherIndex}-Zs{1, stuckBeadIndex} - windCorr(tetherIndex) + geomCorr(tetherIndex);
clampedExt = clampedExt - mean(clampedExt(20:60), "omitnan") + rfit1(3);
instantaneousRates = diff(turns{1})*FPS;
time = [1:length(turns{1})]/FPS;

smoothedExt = smooth(clampedExt, 60);
bucklingZ = f_3piece(rfit1, rfit1(2)+rfit1(4));
FPS = 40;

figure()
subplot(4, 4, [1 5])
plot(lowTurns, lowhat)
xlabel('Magnet Turns', 'FontSize', 14)
ylabel('Extension (nm)', 'FontSize', 14)

subplot(4, 4, [2 3 4 6 7 8])
plot(time, clampedExt, 'Color', [200 200 200]/255)
hold on
plot(time, smooth(clampedExt, 20), 'Color', [0, 0.4470, 0.7410])
ylabel('Extension (nm)', 'FontSize', 14)
title(['Tether ' num2str(tetherIndex-1) ' clamped on stuck bead ' num2str(stuckBeadIndex-1) ', R-R_|_| = ' num2str(geomCorr(tetherIndex)) ' nm, F_\perp = ' num2str(FperpLow, 2) ' pN'], 'FontSize', 18)
subplot(4, 4, [10 11 12])
plot(time, turns{1})
ylabel('Magnet turns', 'FontSize', 14)

subplot(4, 4, [14 15 16])
plot(time, instantaneousRates, '.', 'Color', [200 200 200]/255)
hold on
plot(time, smooth(instantaneousRates, 200), 'Color', [0, 0.4470, 0.7410])
ylabel('Rate (magnet turns/s)', 'FontSize', 14)
xlabel('Time (s)', 'FontSize', 14)

%% Align trace

figure()
hold on
allTimes = [];
allExts = [];
magTurns = 100:35:floor(max(turns{1}));

for t = magTurns
    try
        tempExt = clampedExt(abs(turns{1} - t) < 0.05)-bucklingZ;
        tempSmoothedExt = smoothedExt(abs(turns{1} - t) < 0.05)-bucklingZ;
        crossTime = find(tempSmoothedExt > 0, 1)/FPS;
        time = time/FPS - crossTime;
        plot(time, tempExt, '.')
        allTimes = [allTimes time];
        allExts = [allExts; tempExt];
    catch
    end
end

tbins = -4.95:.1:14.95;
for i = 1:length(tbins)
    meanExt(i) = mean(allExts(allTimes >= tbins(i)-.05 & allTimes < tbins(i)+.05));
end
plot(tbins, meanExt, 'k', 'LineWidth', 2)
plot([-5 15], [0 0], 'k--')
xlim([-5 15])
xlabel('\Deltat from reaching buckling (s)', 'FontSize', 14)
ylabel('\DeltaExtension from buckling (nm)', 'FontSize', 14)

minZ = min(meanExt);
maxZ = max(meanExt);
for i = 1:length(tbins)
    meanTurns(i) = interp1(f_3piece(rfit1, [rfit1(2):.001:(rfit1(2)+50)]),  [rfit1(2):.001:(rfit1(2)+50)], meanExt(i) + bucklingZ)-(rfit1(2)+rfit1(4));
end
x = 1:length(meanTurns);
meanTurns(isnan(meanTurns)) = -(rfit1(4));
%meanTurns(isnan(meanTurns)) = interp1(x(~isnan(meanTurns)),meanTurns(~isnan(meanTurns)),x(isnan(meanTurns)));
figure()
hold on
plot(tbins, meanTurns, 'k', 'LineWidth', 2)
plot([-5 15], [0 0], 'k--')
plot([0 0], [min(ylim) max(ylim)], 'k--')
xlim([-5 15])
xlabel('\Deltat from reaching buckling (s)', 'FontSize', 14)
ylabel('\Deltaturns from buckling (turns)', 'FontSize', 14)

figure()
plot(smooth(meanTurns(2:end), 40), -diff(smooth(meanTurns, 40))./diff(tbins), 'k', 'LineWidth', 2)
hold on
%xlim([-10 10])
ylim([-.5 5])
plot([0 0], [-.5 5], 'k--')
xlabel('\Deltaturns from buckling (turns)', 'FontSize', 14)
ylabel('Relaxation speed (turns/s)', 'FontSize', 14)

save([num2str(tetherIndex-1) ' turns'], 'allTimes', 'allExts', 'meanTurns', 'rfit1', 'lowTurns', 'lowhat')