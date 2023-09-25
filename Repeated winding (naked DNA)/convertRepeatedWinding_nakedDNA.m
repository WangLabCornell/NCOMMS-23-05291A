%% Convert raw data from repeated winding experiment and plot, for naked DNA
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

tetherNumbers = 1:numTethers(1);

%% Geometry correction

geomIndex = 2;

for j = 1:numTethers(1)
    for k = 1:8
        tracegX = Xs{geomIndex, j}-Xs{geomIndex, stuckBeadIndex};
        tracegY = Ys{geomIndex, j}-Ys{geomIndex, stuckBeadIndex};
        Rxs(j, k) = mean(tracegX(turns{geomIndex} > (k-1) * 0.125-.02 & turns{geomIndex} < (k-1) * 0.125+.02), 'omitnan');
        Rys(j, k) = mean(tracegY(turns{geomIndex} > (k-1) * 0.125-.02 & turns{geomIndex} < (k-1) * 0.125+.02), 'omitnan');
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
angleFit = polyfit(magnetAngle(tetherNumbers, 2), magnetAngle(tetherNumbers, 1), 1);
hold on
plot(angleFit(1)*[-500:500] + angleFit(2), [-500:500], 'k--')
theta = pi/2-atan(angleFit(1));
title(['\theta = ' num2str(theta * 180/pi) ' degrees'])

%% Wind to surface correction

surfIndex = 3;

for j = tetherNumbers
    windCorr(j) = mean(Zs{surfIndex, j}(abs(turns{surfIndex}) > max(abs(turns{surfIndex})) - 0.02)) - mean(Zs{surfIndex, stuckBeadIndex}(abs(turns{surfIndex}) > max(abs(turns{surfIndex})) - 0.02));
end

%% Corrections for hat curves

hatIndex = 4;
i = hatIndex;
hatTurns = round(min(turns{hatIndex})):1:round(max(turns{hatIndex}));

for j = tetherNumbers
    for k = 1:numel(hatTurns)
        hats{j}(k) = mean(Zs{i, j}(abs(hatTurns(k) - turns{i}) <= 0.01) - Zs{i, stuckBeadIndex}(abs(hatTurns(k) - turns{i}) <= 0.01), 'omitnan') + geomCorr(j) - windCorr(j);
    end
    [maxZ, maxTurnIndex] = max(hats{j});
    r0 = [-.1 hatTurns(maxTurnIndex) maxZ 18 50];
    try
        rfit{j} = lsqcurvefit(@f_3piece, r0, hatTurns(hats{j} > 500), hats{j}(hats{j} > 500));
        buckling(j) = rfit{j}(4)+rfit{j}(2);
    catch
        rfit{j} = [nan nan nan nan nan];
    end
end

%% Force calibrations

forceIndex = 5;
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

%% Plot
time = [1:numel(topoZs{tetherNumbers(1)})]/FPS;

for j = tetherNumbers
    figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    subplot(1, 3, 1)
    hold on
    plot(hatTurns, hats{j}, 'LineWidth', 2)
    plot([-60:.1:60], f_3piece(rfit{j}, [-60:.1:60]))
    ylim([0 4000])
    ylabel('Extension (nm)', 'FontSize', 16)
    xlabel('Magnet Turns', 'FontSize', 16)
    subplot(1, 3, [2:3])
    plot(time, topoZs{j})
    hold on
    plot([min(xlim) max(xlim)], [1 1]*f_3piece(rfit{j}, 30), 'k--')
    ylim([0 4000])
    xlabel('Time (s)', 'FontSize', 12)
    ylabel('Extension (nm)', 'FontSize', 14)

    sgtitle(['Tether ' num2str(j-1) ', F = ' num2str(Fpara(j), 3) ' pN, ' num2str(meanRates(j), 3) ' turns/s'])
    print(['X:\jl3452\MATLAB\MT3 repeated winding for prebuckled processivity\plots\tether ' num2str(j-1) '.png'],'-dpng','-r0');
    close;
end