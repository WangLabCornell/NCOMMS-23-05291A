function [Height,hatcenter, buckling_negative, buckling_positive, slope_negative, slope_positive] = HCpara(r0)
%obtains hat curve parameters from the result of a 5 piece fit
% buckling positive is buckling point relative to the hat center

r = reshape(r0,4,2);

M_tl = [2*r(2,1) 1 0;  r(1,1)^2  r(1,1) 1; r(2,1)^2  r(2,1) 1];
C_tl = [(r(3,2)-r(2,2))/(r(3,1)-r(2,1)) ;  r(1,2);  r(2,2)];
para_tl= M_tl\C_tl;

M_tr = [2*r(3,1) 1 0;  r(4,1)^2  r(4,1) 1;  r(3,1)^2  r(3,1) 1];
C_tr = [(r(3,2)-r(2,2))/(r(3,1)-r(2,1)) ;  r(4,2);  r(3,2)];
para_tr = M_tr\C_tr;
alpha_n = 2 * para_tl(1) * r(1,1) + para_tl(2);
alpha_p = 2 * para_tr(1) * r(4,1) + para_tr(2);
slope_negative = alpha_n;
slope_positive = alpha_p;
alpha_c = (r(3,2)-r(2,2))/(r(3,1)-r(2,1));

buckling_negative = (r(2,2) - r(1,2) + alpha_n * r(1,1) - alpha_c * r(2,1))/(alpha_n-alpha_c);
buckling_positive = (r(4,2) - r(3,2) - alpha_p * r(4,1) + alpha_c * r(3,1))/(alpha_c-alpha_p);

%% top left curve
M_tl = [2*r(2,1) 1 0;
        r(1,1)^2  r(1,1) 1;
        r(2,1)^2  r(2,1) 1];
C_tl = [(r(3,2)-r(2,2))/(r(3,1)-r(2,1)) ; 
        r(1,2);
        r(2,2)];
para_tl= M_tl\C_tl;
Height = para_tl(3) - para_tl(2)^2/(4*para_tl(1));
hatcenter =  -para_tl(2)/2/para_tl(1);
buckling_positive = buckling_positive - hatcenter;
buckling_negative = buckling_negative - hatcenter;
end