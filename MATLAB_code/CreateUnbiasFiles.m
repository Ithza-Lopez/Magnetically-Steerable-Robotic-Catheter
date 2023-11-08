%% create xlsx files from csv
clear
clc

file = dir('*.csv'); 
s= size(file,1);
for i= 1:s
    Data = readtable(file(i).name);
    filename=file(i).name;
    filename= filename(1:end-4); 
    writetable(Data,[filename '.xlsx']);
end

disp('XLSX files are ready')
%% Intializing magnet file names
% Write names of files below
clear
clc

exp_date = '5-16-23';
magnet_type = 'A';
tilt_angle = 30;
% columns_labels = ["C"];
columns_labels = [ "B", "C", "D"];
elevation = [0];

switch magnet_type
    case 'A'
        rows_labels = 8:-1:-8;
    case 'B'
        rows_labels = 8:-1:-8;
    case 'C'
        rows_labels = 6:-1:-6;
end

earth_field_filename = append(exp_date, "-EarthField.xlsx");
num_slots = length(columns_labels) * length(rows_labels);
magnet_filenames = [];
no_magnet_filename = append(exp_date, "-0.xlsx");
for height = 1:length(elevation)
    height_num = elevation(height);
for col = 1:length(columns_labels)
    col_letter = columns_labels(col);
    for row = rows_labels
%         if tilt_angle >0
%             filename = append(exp_date,'-',magnet_type, '-',col_letter,'-',num2str(row),'-', num2str(tilt_angle),'-',num2str(height_num), '.xlsx');
%         else
            filename = append(exp_date,'-',magnet_type, '-',col_letter,'-',num2str(row),'-', num2str(tilt_angle),'-',num2str(height_num), '.xlsx');
        magnet_filenames = [magnet_filenames filename];
    end 
end 
end
% %Important to order filenames from top to bottom of the bridge (+8 to -8)
% magnet_filenames = {dir('*.xlsx').name};
% earth_field_filename = append(exp_date, "-EarthField.xlsx");
% magnet_filenames(strcmp(magnet_filenames, earth_field_filename)) = [];

% no_magnet_filename = "2-22-23-0.xlsx";
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
