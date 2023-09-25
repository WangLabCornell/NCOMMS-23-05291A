function r = fit5piece(turn, zBead)
    
    height = max(zBead);
    indexmax = find(zBead == max(zBead)); indexmax = indexmax(1);
    peakTurn = mean(turn(indexmax));
    F_positive = @(para, x) (height + para(1) * (x-peakTurn)) .* (x <= para(3)) + ((height + para(1) * (para(3) - peakTurn)) + para(2) * (x - para(3))) .* (x > para(3));
    para0 = [-0.004, -0.020 , 40+peakTurn];
    [para_positive,~,~,~,~, ~, ~] = lsqcurvefit(F_positive,para0,turn(turn >= peakTurn), zBead(turn >= peakTurn));
    turn_buckling_positive = para_positive(3) ;
    
    F_negative = @(para, x) (height + para(1) * (x - peakTurn)) .* (x >= para(3)) + ((height + para(1) * (para(3)-peakTurn)) + para(2) * (x - para(3))) .* (x < para(3));
    para0 = [0.001, 0.020 , -10+peakTurn];
    [para_negative,~,~,~,~, ~, ~] = lsqcurvefit(F_negative,para0, turn(turn < peakTurn), zBead(turn < peakTurn));
    turn_buckling_negative = para_negative(3);
    
     r0 = [ turn_buckling_negative F_negative(para_negative,turn_buckling_negative);
                peakTurn F_negative(para_negative,peakTurn);
                turn_buckling_positive-1 F_positive(para_positive,turn_buckling_positive-1);
                turn_buckling_positive+1 F_positive(para_positive,turn_buckling_positive+1)];
    
%     r0 = [-10 1200; 
%            0 1280;
%            45 1180;
%            55 1050];
    r = lsqcurvefit(@f_5piece,r0,turn,zBead);
    
end