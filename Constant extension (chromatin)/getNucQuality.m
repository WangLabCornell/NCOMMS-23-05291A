function [nucQuality,n_nuc] = getNucQuality(r0,F)
%obtains goodness of a nucleosome array based on Height and widht of hat curve.
% N Vs Height and Width Relationship is obtained from MT data and is summarized in 210108 summary;
% Width goodness range is calculated and summarized on 210220. 

r = reshape(r0,4,2);
slope = -0.0414;
intercept = 3.2306;
[Height, ~, ~, buckling_positive, ~,~] = HCpara(r);

n_nuc =  (Height - intercept)/slope;
% n_nucRange = [n_nuc-5 n_nuc+5];
%wFit = [0.5973 22.7936];
pLow = [0.5212 19.26];
pHigh = [0.6734 26.33];
widthRange = [polyval(pLow,n_nuc) polyval(pHigh,n_nuc)];

if buckling_positive > widthRange(1) && buckling_positive < widthRange(2) %abs(n_nuc  - n_nuc_width) < dn_nuc_width
    nucQuality = 'good';
else
    nucQuality = 'bad';
end

end