%% Convert raw data from constant extension clamp and plot, for chromatin
% Needs the following files: 
% 1. CircleFitByPratt (https://www.mathworks.com/matlabcentral/fileexchange/22643-circle-fit-pratt-method)
% 2. fit5piece & f_5piece: functions for fitting chromatin hat curve to 5-piece function
% 3. HCpara: calculates hat curve parameters from fit5piece results
% 4. getNucQuality: determines chromatin quality from HCpara results

tetherIndex = 1+1; %index of clamped tether
stuckBeadIndex = 0+1; %index of stuck used for drift correction (fiducial marker of surface)
FPS = 40; %data frame rate
 
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
        Rxs(k) = mean(tracegX(turns{geomIndex} > (k-1) * 0.125-.02 & turns{geomIndex} < (k-1) * 0.125+.02), 'omitnan');
        Rys(k) = mean(tracegY(turns{geomIndex} > (k-1) * 0.125-.02 & turns{geomIndex} < (k-1) * 0.125+.02), 'omitnan');
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

prehatIndex = 4;
posthatIndex = 5;

for j = tetherIndex
    hatTurns = round(min(turns{prehatIndex})):round(max(turns{prehatIndex}));
    postTurns = turns{posthatIndex} - round(turns{posthatIndex}(1));
    for k = 1:numel(hatTurns)
        prehat(k) = mean(Zs{prehatIndex, j}(abs(turns{prehatIndex} - hatTurns(k)) < 0.02)) - mean(Zs{prehatIndex, stuckBeadIndex}(abs(turns{prehatIndex} - hatTurns(k)) < 0.02)) - windCorr(j) + geomCorr(j);
        posthat(k) = mean(Zs{posthatIndex, j}(abs(postTurns - hatTurns(k)) < 0.02)) - mean(Zs{posthatIndex, stuckBeadIndex}(abs(postTurns - hatTurns(k)) < 0.02)) - windCorr(j) + geomCorr(j);
    end
    rfitpre = fit5piece(hatTurns, prehat/1000);
    [height, hatCenterLow, bucklingm, bucklingp, slopem, slopep] = HCpara(rfitpre);
    [nucQualityPre,n_nucPre] = getNucQuality(rfitpre, 0.5);
    rfitpost = fit5piece(hatTurns, posthat/1000);
    [heightPost, hatCenterPost, bucklingm2, bucklingp2, slopem2, slopep2] = HCpara(rfitpost);
    [nucQualityPost, n_nucPost] = getNucQuality(rfitpost, 0.5);
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

%% Calculate rates

clampedExt = Zs{1, tetherIndex}-Zs{1, stuckBeadIndex} - windCorr(tetherIndex) + geomCorr(tetherIndex);
instantaneousRates = diff(turns{1})*FPS;
time = [1:length(turns{1})]/FPS;


%% Plot trace

figure()
subplot(4, 4, [1 5])
plot(hatTurns, prehat)
xlim([-40 70])
ylim([0 1500])
xlabel('Magnet Turns', 'FontSize', 18)
ylabel('Extension (nm)', 'FontSize', 18)
hold on
plot([-40:.1:70], f_5piece(rfitpre, [-40:.1:70]), 'k')
plot([1 1]* (bucklingm+hatCenterLow), [1 1] * f_5piece(rfitpre, bucklingm+hatCenterLow)*1000 + [-100 100], 'k--')
plot([1 1]* (bucklingp+hatCenterLow), [1 1] * f_5piece(rfitpre, bucklingp+hatCenterLow)*1000 + [-100 100], 'r--')

subplot(4, 4, [2 3 4 6 7 8])
plot(time, clampedExt, 'Color', [200 200 200]/255)
hold on
plot(time, smooth(clampedExt, 20), 'Color', [0, 0.4470, 0.7410])
ylabel('Extension (nm)', 'FontSize', 18)
title(['Tether ' num2str(tetherIndex-1) ' clamped on stuck bead ' num2str(stuckBeadIndex-1) ', R-R_|_| = ' num2str(geomCorr(tetherIndex)) ' nm, F_\perp = 0.522 pN'], 'FontSize', 18)
plot([min(xlim) max(xlim)], [1 1]*f_5piece(rfitpre, hatCenterLow+bucklingp)*1000, 'r--')
plot([min(xlim) max(xlim)], [1 1]*f_5piece(rfitpre, hatCenterLow+bucklingm)*1000, 'k--')
ylim([0 1500])

subplot(4, 4, [10 11 12])
plot(time, turns{1})
ylabel('Magnet turns', 'FontSize', 18)

subplot(4, 4, [14 15 16])
plot(time, instantaneousRates, '.', 'Color', [200 200 200]/255)
hold on
plot(time, smooth(instantaneousRates, 400), 'Color', [0, 0.4470, 0.7410])
ylim([-1 6])
ylabel('Rate (turns/s)', 'FontSize', 18)
xlabel('Time (s)', 'FontSize', 18)

%% Plot hats before and after clamp

figure()
plot(hatTurns, prehat)
hold on
plot(hatTurns, posthat)
ylim([0 1500])
plot([1 1] * (bucklingm+hatCenterLow), [min(ylim) max(ylim)], 'k--')
plot([1 1] * (bucklingp+hatCenterLow), [min(ylim) max(ylim)], 'k--')
xlabel('Magnet turns', 'FontSize', 18)
ylabel('Extension (nm)', 'FontSize', 18)
leg = legend(['Before clamp @ F_\perp = ' num2str(FperpLow) ' pN, ' nucQualityPre ', ' num2str(n_nucPre, 4) ' nucs'], ['After clamp, ' nucQualityPost ', ' num2str(n_nucPost, 4) ' nucs']);
leg.Location = 'south';