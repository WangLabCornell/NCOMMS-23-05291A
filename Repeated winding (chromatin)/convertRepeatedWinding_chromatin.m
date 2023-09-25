%% Convert raw data from repeated winding experiment and plot, for chromatin
% Needs the following files: 
% 1. CircleFitByPratt (https://www.mathworks.com/matlabcentral/fileexchange/22643-circle-fit-pratt-method)
% 2. fit5piece & f_5piece: functions for fitting chromatin hat curve to 5-piece function
% 3. HCpara: calculates hat curve parameters from fit5piece results
% 4. getNucQuality: determines chromatin quality from HCpara results

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

tetherNumbers = 1:numTethers(1);

%% Geometry correction

geomIndex = 2;

for j = 1:numTethers(1)
    for k = 1:8
        tracegX = Xs{geomIndex, j}-Xs{geomIndex, stuckBeadIndex};
        tracegY = Ys{geomIndex, j}-Ys{geomIndex, stuckBeadIndex};
        Rxs(j, k) = nanmean(tracegX(turns{geomIndex} > (k-1) * 0.125-.02 & turns{geomIndex} < (k-1) * 0.125+.02));
        Rys(j, k) = nanmean(tracegY(turns{geomIndex} > (k-1) * 0.125-.02 & turns{geomIndex} < (k-1) * 0.125+.02));
    end
    circleFit = CircleFitByPratt([(Rxs(j, :)-mean(Rxs(j, :)))' (Rys(j, :)-mean(Rys(j, :)))']);
    Rpara(j) = max([min([circleFit(3) 500]) 0]);
    Rperp = sqrt(500^2 - Rpara(j)^2);
    geomCorr(j) = 500 - Rperp;
    magnetAngle(j, 1) = Rxs(j, 1)-mean(tracegX);
    magnetAngle(j, 2) = Rys(j, 1)-mean(tracegY);
end

figure()
plot(magnetAngle(tetherNumbers, 1), magnetAngle(tetherNumbers, 2), '.')
fitXs = magnetAngle(tetherNumbers, 2);
fitYs = magnetAngle(tetherNumbers, 1);
fitXs(sqrt(magnetAngle(tetherNumbers, 2).^2+magnetAngle(tetherNumbers, 1).^2) > 550) = [];
fitYs(sqrt(magnetAngle(tetherNumbers, 2).^2+magnetAngle(tetherNumbers, 1).^2) > 550) = [];
angleFit = polyfit(fitXs, fitYs, 1);
hold on
plot(angleFit(1)*[-500:500] + angleFit(2), [-500:500], 'k--')
theta = pi/2-atan(angleFit(1));
title(['\theta = ' num2str(theta * 180/pi) ' degrees'])

%% Wind to surface correction

surfIndex = 3;

figure()
for j = tetherNumbers
    windCorr(j) = mean(Zs{surfIndex, j}(abs(turns{surfIndex}) > max(abs(turns{surfIndex})) - 0.02)) - mean(Zs{surfIndex, stuckBeadIndex}(abs(turns{surfIndex}) > max(abs(turns{surfIndex})) - 0.02));
end

%% Corrections for hat curves
prehatIndex = 4;
posthatIndex = 5;

for j = tetherNumbers
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
    prehats{j} = prehat;
    posthats{j} = posthat;
    prenucQs{j} = nucQualityPre;
    postnucQs{j} = nucQualityPost;
    prenucNs(j) = n_nucPre;
    postnucNs(j) = n_nucPost;
    rfitpres{j} = rfitpre;
    rfitposts{j} = rfitpost;
end

%% Force calibrations

forceIndex = 6;
for j = tetherNumbers
    extension(j) = mean(Zs{forceIndex, j}) - mean(Zs{forceIndex, stuckBeadIndex}) + geomCorr(j) - windCorr(j);
    Xpara = Xs{forceIndex, j} * cos(theta) + Ys{forceIndex, j} * sin(theta);
    Xperp = Xs{forceIndex, j} * -sin(theta) + Ys{forceIndex, j} * cos(theta);
    Fpara(j) = 4.114 / var(Xpara) * extension(j);
    Fperp(j) = 4.114 / var(Xperp) * (extension(j) + 500);
end

%% Corrections for topoZs

for j = tetherNumbers
    topoZs{j} = Zs{1, j} - Zs{1, stuckBeadIndex} + geomCorr(j) - windCorr(j);
end
time = [1:numel(topoZs{tetherNumbers(1)})]/FPS;

for j = tetherNumbers
    if(strcmp(prenucQs{j}, 'good') == 1 && prenucNs(j) <= 55 && prenucNs(j) >= 45)
        figure('units', 'normalized', 'outerposition', [0 0 1 1]);
        subplot(1, 3, 1)
        hold on
        plot(hatTurns, prehats{j}, 'LineWidth', 2)
        plot(hatTurns, posthats{j}, 'LineWidth', 2)
        leg = legend(['Pre, ' prenucQs{j} ', N_n_u_c = ' num2str(prenucNs(j), 3)], ['Post, ' postnucQs{j} ', N_n_u_c = ' num2str(postnucNs(j), 3)]);
        ylim([0 1500])
        ylabel('Extension (nm)', 'FontSize', 16)
        xlabel('Magnet Turns', 'FontSize', 16)
        subplot(1, 3, [2:3])
        plot(time, topoZs{j})
        hold on
        ylim([0 1500])
        xlabel('Time (s)', 'FontSize', 12)
        ylabel('Extension (nm)', 'FontSize', 14)
        suptitle(['Tether ' num2str(j-1) ', F = ' num2str(Fpara(j), 3) ' pN'])
        print(['tether ' num2str(j-1) '.png'],'-dpng','-r0');
        close;
    end
end