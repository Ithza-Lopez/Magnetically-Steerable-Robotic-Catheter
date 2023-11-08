function [B_meas_mod] = mode_B_meas_loc(O_all,B_meas)
% input:
%        O_all: 3 x 3 x N_sensor
%        B_meas: 3*N_sensor x 1
N = length(B_meas)/3;
B_meas_mod = zeros(length(B_meas),1);
for i = 1:N
    
    B_meas_mod((i*3-2):(i*3)) = O_all(:,:,i)*B_meas((i*3-2):(i*3));
    
end

end 