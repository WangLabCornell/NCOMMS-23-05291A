%% Convert raw data from pre-buckled relaxation experiment and plot
% Needs the following files: 
% 1. CircleFitByPratt (https://www.mathworks.com/matlabcentral/fileexchange/22643-circle-fit-pratt-method)
% 2. f_3piece: 3-piece function to fit to naked DNA hat curve

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
        Rxs(k) = nanmean(tracegX(turns{geomIndex} > (k-1) * 0.125-.02 & turns{geomIndex} < (k-1) * 0.125+.02));
        Rys(k) = nanmean(tracegY(turns{geomIndex} > (k-1) * 0.125-.02 & turns{geomIndex} < (k-1) * 0.125+.02));
    end
    circleFit = CircleFitByPratt([(Rxs-mean(Rxs))' (Rys-mean(Rys))']);
    Rpara = max([min([circleFit(3) 500]) 0]);
    Rperp = sqrt(500^2 - Rpara^2);
    geomCorr(j) = 500 - Rperp;
    magnetAngle(j, 1) = Rxs(1)-mean(tracegX);
    magnetAngle(j, 2) = Rys(1)-mean(tracegY);
end

tempAngles = [];
for j = 1:numTethers(1)
    if(sqrt(magnetAngle(j, 1)^2+magnetAngle(j, 2)^2) < 500 && sqrt(magnetAngle(j, 1)^2+magnetAngle(j, 2)^2) > 25) 
        tempAngles = [tempAngles; magnetAngle(j, :)];
        display(j)
    end
end
magnetAngle = tempAngles;

figure()
plot(magnetAngle(:, 1), magnetAngle(:, 2), '.')
angleFit = polyfit(magnetAngle(:, 1), magnetAngle(:, 2), 1);
hold on
plot([-500:500], angleFit(1) * [-500:500] + angleFit(2), 'k--')
theta = atan(angleFit(1));
title(['\theta = ' num2str(theta * 180/pi) ' degrees'])

%% Wind to surface correction

surfIndex = 3;

for j = tetherIndex
    windCorr(j) = mean(Zs{surfIndex, j}(end-39:end)) - mean(Zs{surfIndex, stuckBeadIndex}(end-39:end));
    extension(j) = mean(Zs{surfIndex, j}(turns{surfIndex} == 0)) - mean(Zs{surfIndex, j}(turns{surfIndex} > max(turns{surfIndex}) - 0.02));
end

%% Process hats

hatIndex = 4;

for j = tetherIndex
    hatTurnsLow = round(min(turns{hatIndex})):round(max(turns{hatIndex}));
    for k = 1:numel(hatTurnsLow)
        prehatLow(k) = mean(Zs{hatIndex, j}(abs(turns{hatIndex} - hatTurnsLow(k)) < 0.02)) - mean(Zs{hatIndex, stuckBeadIndex}(abs(turns{hatIndex} - hatTurnsLow(k)) < 0.02)) - windCorr(j) + geomCorr(j);
    end
end

hatIndex = 5;

for j = tetherIndex
    hatTurns = round(min(turns{hatIndex})):round(max(turns{hatIndex}));
    for k = 1:numel(hatTurns)
        prehat(k) = mean(Zs{hatIndex, j}(abs(turns{hatIndex} - hatTurns(k)) < 0.02)) - mean(Zs{hatIndex, stuckBeadIndex}(abs(turns{hatIndex} - hatTurns(k)) < 0.02)) - windCorr(j) + geomCorr(j);
    end
end

%% Fit hats

[maxZ, maxTurnIndex] = max(prehatLow);
r0 = [-.5 hatTurnsLow(maxTurnIndex) maxZ 20 50];
try
    rfitLow = lsqcurvefit(@f_3piece, r0, hatTurnsLow(prehatLow > 1000), prehatLow(prehatLow > 1000));
    rfitLow_c = rfitLow;
    rfitLow_c(2) = 0;
catch
    rfitLow = [nan nan nan nan nan];
    rfitLow_c = rfitLow;
end

[maxZ, maxTurnIndex] = max(prehat);
r0 = [-.5 hatTurns(maxTurnIndex) maxZ 20 50];
try
    rfit = lsqcurvefit(@f_3piece, r0, hatTurns(prehat > 1000), prehat(prehat > 1000));
    rfit_c = rfit;
    rfit_c(2) = rfit_c(2)-rfitLow(2);
catch
    rfit = [nan nan nan nan nan];
    rfit_c = rfit;
end

%% Force calibrations

forceIndex = 6;

for j = tetherIndex
    extension = mean(Zs{forceIndex, j}) - mean(Zs{forceIndex, stuckBeadIndex}) + geomCorr(j) - windCorr(j);
    Xpara = Xs{forceIndex, j} * cos(theta) + Ys{forceIndex, j} * sin(theta);
    Xperp = Xs{forceIndex, j} * -sin(theta) + Ys{forceIndex, j} * cos(theta);
    FparaLow = 4.114 / var(Xpara) * extension;
    FperpLow = 4.114 / var(Xperp) * (extension + 500);
end

forceIndex = 7;

for j = tetherIndex
    extension = mean(Zs{forceIndex, j}) - mean(Zs{forceIndex, stuckBeadIndex}) + geomCorr(j) - windCorr(j);
    Xpara = Xs{forceIndex, j} * cos(theta) + Ys{forceIndex, j} * sin(theta);
    Xperp = Xs{forceIndex, j} * -sin(theta) + Ys{forceIndex, j} * cos(theta);
    Fpara = 4.114 / var(Xpara) * extension;
    Fperp = 4.114 / var(Xperp) * (extension + 500);
end

topoExt = Zs{1, tetherIndex}-Zs{1, stuckBeadIndex} - windCorr(tetherIndex) + geomCorr(tetherIndex);
instantaneousRates = diff(turns{1})*FPS;
time = [1:length(turns{1})]/FPS;

%% Plot trace

figure()
subplot(3, 4, [1 5])
plot(hatTurns, prehat)
hold on
plot([-60:.1:60], f_3piece(rfit, [-60:.1:60]), 'k--')
ylim([500 4000])
xlim([-15 100])
plot([1 1]*rfitLow(2), [min(ylim) max(ylim)], 'k--')
xlabel('Magnet Turns', 'FontSize', 14)
ylabel('Extension (nm)', 'FontSize', 14)

subplot(3, 4, [2 3 4 6 7 8])
plot(time, topoExt, 'Color', [200 200 200]/255)
hold on
plot(time, smooth(topoExt, 5), 'Color', [0, 0.4470, 0.7410])
plot([min(xlim) max(xlim)], [1 1]*f_3piece(rfit_c, rfit_c(2)+rfit_c(4)), 'k--')
plot([min(xlim) max(xlim)], [1 1]*f_3piece(rfit_c, 45), 'k--')
plot([min(xlim) max(xlim)], [1 1]*f_3piece(rfit_c, rfit_c(2)+rfit_c(4)+45), 'k--')
ylim([500 4000])
ylabel('Extension (nm)', 'FontSize', 14)
title(['Tether ' num2str(tetherIndex-1) ' clamped on stuck bead ' num2str(stuckBeadIndex-1) ', R-R_|_| = ' num2str(geomCorr(tetherIndex)) ' nm, F_\perp = ' num2str(Fperp, 3) ' pN'], 'FontSize', 18)

subplot(3, 4, [10 11 12])
plot(time, turns{1})
ylabel('Magnet turns', 'FontSize', 14)

%% Calculate pre-buckled turns

delTurns = 35;
relTurns = [100:delTurns:(round(max(turns{1})))];

waitTimes = [0 2 4 6 8 10 15 2 4 6 8 10 15 2 4 6 8 10 15 2 4 6 8 10 15 2 4 6 8 10 15];
if(numel(waitTimes) ~= numel(relTurns))
    relTurns(end) = [];
    waitTimes = waitTimes(1:numel(relTurns));
end

for n = 1:numel(relTurns)
    relZs{n} = smooth(topoExt(abs(turns{1}-relTurns(n)) < 0.1), 20);
    relTs{n} = (1:numel(relZs{n}))/FPS;
    try
        buckledRate(n) = 40*(10)/(find(relZs{n} < f_3piece(rfit_c, rfit_c(2)+rfit_c(4) + 2), 1, 'last')-find(relZs{n} > f_3piece(rfit_c, rfit_c(2)+rfit_c(4)+12), 1, 'first'));
    catch
        buckledRate(n) = nan;
    end
    try
        prebuckledTurns(n) = (interp1(f_3piece(rfit_c, [0:.01:100]), 0:.01:100, relZs{n}(19))-delTurns+buckledRate(n-1)*(0.5+delTurns/40));
    catch
        prebuckledTurns(n) = nan;
    end
end

settings = [2 4 6 8 10 15];
waitTimes = waitTimes(2:(numel(prebuckledTurns(~isnan(prebuckledTurns)))+1));
prebuckledTurns = prebuckledTurns(~isnan(prebuckledTurns));

for i = 1:numel(settings)
    try
        meanTurnStateReached(i) = nanmean(prebuckledTurns(settings(i) == waitTimes));
        stdTurnStateReached(i) = nanstd(prebuckledTurns(settings(i) == waitTimes));
    catch
        meanTurnStateReached(i) = nan;
        stdTurnStateReached(i) = nan;
    end
end
meanRate = 10*numel(buckledRate(2:end-1)) / sum(10 ./ buckledRate(2:end-1));

figure()
buckling = rfit_c(2)+rfit_c(4);
plot(waitTimes, prebuckledTurns, 'x')
hold on
errorbar(settings, meanTurnStateReached, stdTurnStateReached)
plot([0 18], [1 1]*(buckling), 'k--');
plot([0 buckling/nanmean(buckledRate(2:end-1))], [buckling 0], 'g--')
ylim([0 40])
leg = legend('Data', '\mu \pm \sigma', 'Buckling', 'Buckled rate');
leg.FontSize = 14;
leg.Location = 'northeast';
xlabel('Time waited after reaching buckling (s)', 'FontSize', 14)
ylabel('Turn state reached', 'FontSize', 14)
title(['Tether ' num2str(tetherIndex-1) ', F_\perp = ' num2str(Fperp, 3) ' pN, buckled topo speed = ' num2str(nanmean(buckledRate(2:end-1)), 3) ' turns/s'], 'FontSize', 14)
