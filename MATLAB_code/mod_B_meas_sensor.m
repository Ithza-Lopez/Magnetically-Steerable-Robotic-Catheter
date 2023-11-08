function [B_meas_mod] = mod_B_meas_sensor(O,B_meas)
% input:
%        O: 3x3 matrix
%        B_meas: measurement vector (3xN)x1

% output:
%        B_meas_mod: modified   (3xN)x1
N = length(B_meas)/3;
B_meas_mod = O*reshape(B_meas,[3,N]);
B_meas_mod = reshape(B_meas_mod,[3*N,1]);




end