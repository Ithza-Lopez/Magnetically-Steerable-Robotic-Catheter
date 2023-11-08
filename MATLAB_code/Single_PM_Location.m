function [pm_location] = Single_PM_Location(magnet_type, col_letter, row_num, height)
% input
% magnet_type: 'A', 'B', or 'C'
% col_letter: 'A', 'B', or 'C', A is left, B is center, C is right
% row_num: row of the magnet(-8 to +8)
% location of slots flip flop from pos to neg

%output
%pm_location [x;y;z] with respect to origin at the center
if strcmp(magnet_type, 'A')
%     disp('reading magnet A')
    hor_inter = 10; %mm
    ver_inter = 10;
    num_rows = 17;
elseif strcmp(magnet_type, 'B')
%     disp('reading magnet B')
    hor_inter = 10; %mm
    ver_inter = 10;
    num_rows = 17;
elseif strcmp(magnet_type, 'C')
%     disp('reading magnet C')
    hor_inter = 10; %mm
    ver_inter = 12;
    num_rows = 13;
end

if strcmp(col_letter, 'B')
%     disp('reading column A')
    x_loc = -hor_inter;
elseif strcmp(col_letter, 'C')
%     disp('reading column A')
    x_loc = 0;
elseif strcmp(col_letter, 'D')
%     disp('reading column A')
    x_loc = hor_inter;
end

y_loc = -row_num*ver_inter;
z_loc = height;

pm_location = [x_loc; y_loc; z_loc]*1e-3;

end