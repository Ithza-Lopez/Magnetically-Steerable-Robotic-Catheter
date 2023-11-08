%% Initialization
clear
clc
%taking average B_meas values from one sensor with multiple magnet
%positions

% Write names of files below
% sensor_of_interest = #; %from sensor layout ppt
sensor_of_interest = 1;

no_magnet_filename = "2-9-23-0.xlsx";

%Important to order filenames from top to bottom of the bridge (+8 to -8)
magnet_filenames = ["2-9-23-A-B-8.xlsx", "2-9-23-A-B-7.xlsx", "2-9-23-A-B-6.xlsx",...
    "2-9-23-A-B-5.xlsx","2-9-23-A-B-4.xlsx", "2-9-23-A-B-3.xlsx", "2-9-23-A-B-2.xlsx",...
    "2-9-23-A-B-1.xlsx","2-9-23-A-B-0.xlsx", "2-9-23-A-B--1.xlsx","2-9-23-A-B--2.xlsx",...
    "2-9-23-A-B--3.xlsx","2-9-23-A-B--4.xlsx","2-9-23-A-B--5.xlsx","2-9-23-A-B--6.xlsx",...
    "2-9-23-A-B--7.xlsx","2-9-23-A-B--8.xlsx"];
earth_field_filename = "2-9-23-A-EarthField.xlsx";

% no_magnet_filename = "2-9-23-0.xlsx";
% magnet_filenames = ["2-11-23-C-B-6.xlsx",...
%     "2-11-23-C-B-5.xlsx","2-11-23-C-B-4.xlsx", "2-11-23-C-B-3.xlsx", "2-11-23-C-B-2.xlsx",...
%     "2-11-23-C-B-1.xlsx","2-11-23-C-B-0.xlsx", "2-11-23-C-B--1.xlsx","2-11-23-C-B--2.xlsx",...
%     "2-11-23-C-B--3.xlsx","2-11-23-C-B--4.xlsx","2-11-23-C-B--5.xlsx","2-11-23-C-B--6.xlsx"];
% earth_field_filename = "2-11-23-C-B-EarthField.xlsx";

%% Getting unbias filenames
num_magnet_tests = length(magnet_filenames);
unbias_filenames = strings(1, num_magnet_tests);

% Getting unbias filenames
for magnet_position = 1:num_magnet_tests
    magnet_setting = magnet_filenames(magnet_position);
    magnet_setting = extractBefore(magnet_setting, ".");
    unbias_name = append(magnet_setting, "-Unbias.xlsx");
    unbias_filenames(magnet_position) = unbias_name;
end

%% Populating files with processed data (magnet - earth field)
% this section takes a minute to process all files,
for magnet_position = 1:num_magnet_tests
    disp(['Creating Unbias Sensing for Magnet Position #', num2str(magnet_position)])
    magnet_reading = magnet_filenames(magnet_position);
    unbias_reading = unbias_filenames(magnet_position);
    UnbiasSensing(no_magnet_filename,magnet_reading,earth_field_filename, unbias_reading )
end

disp("Files updated!")

%% Extracting B_meas from same sensor
B_measure = zeros(1,num_magnet_tests*3);
for magnet_position = 1:num_magnet_tests
    unbias_slot_name = unbias_filenames(magnet_position);
    Unbias_avg = readtable(unbias_slot_name, 'Sheet', 'Avg_vector');
    B_ave = table2array(Unbias_avg(:,(sensor_of_interest*3-2):1:(sensor_of_interest*3)))*1e-4;
    if magnet_position == 1
        B_measure(:, magnet_position) = B_ave(1);
        B_measure(:, magnet_position+1) = B_ave(2);
        B_measure(:, magnet_position+2) = B_ave(3);
    elseif magnet_position >1
        B_measure(:, magnet_position*3-2) = B_ave(1);
        B_measure(:, magnet_position*3-1) = B_ave(2);
        B_measure(:, magnet_position*3) = B_ave(3);
    end
end
save('B_measure.mat','B_measure');
save('Sensor_of_interest.mat', 'sensor_of_interest');



