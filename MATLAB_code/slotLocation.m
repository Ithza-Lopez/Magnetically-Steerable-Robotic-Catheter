function [pm_all] = slotLocation(magnet_type, which_rows, M)
%location of slots flip flop from pos to neg

switch magnet_type
    case 'A'
        %     disp('reading A')
        hor_inter = 10; %mm
        ver_inter = 10;
        num_rows = 17;
    case 'B'
        %     disp('reading B')
        hor_inter = 10; %mm
        ver_inter = 10;
        num_rows = 17;
    case 'C'
        %     disp('reading C')
        hor_inter = 10; %mm
        ver_inter = 12;
        num_rows = 13;
    case 'D'
        hor_inter = 10; %mm
        ver_inter = 12;
        num_rows = 13; 
        phi = pi/4;
end


%print out from A1 to C#
%assuming odd number columns
if strcmp(which_rows, 'all')
%     disp('reading all')
    num_slots = num_rows * 3;

elseif strcmp(which_rows, 'center')
%     disp('reading center')
    num_slots = num_rows * 1;
    
end

pm_all = zeros(3,num_slots);
x_loc = 0;
z_loc = 0;

for row = 1:num_slots
    y_loc = (((num_slots-1)/2)-row+1)*ver_inter;
    pm_loc = [x_loc; y_loc; z_loc];
    pm_all(:, row) = pm_loc;
end
