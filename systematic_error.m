function [sys_err_density,sys_err_AMA]=systematic_error(diam_km,Dmin_fit,b_at_Dmin_fit,isCalibrated,chronology)
D=diam_km(diam_km>Dmin_fit);
N=length(D);
if isCalibrated
    sys_err_diameter=0.02;
else
    sys_err_diameter=0.08;
end
sys_err_density=(1+sys_err_diameter)^abs(b_at_Dmin_fit);