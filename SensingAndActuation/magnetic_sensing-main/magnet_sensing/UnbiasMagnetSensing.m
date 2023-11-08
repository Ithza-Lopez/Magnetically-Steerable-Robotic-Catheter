%% Initializing magnet and no magnet data
clear all
clc
% Write names of files below
no_magnet_filename = '1-28-23 No Magnet 2G 2';
center_magnet_filename = '1-28-23 Small Magnet Centered 2G 2';
earth_field_filename = '1-28-23 EarthMagneticField Small 2G .xlsx';
unbias_filename = '1-28-23 UnbiasMagnetSensing Small 2G 2.xls';

%%initializing data from files
no_magnet_dict = sensorDataRead(no_magnet_filename);
no_magent_average = readtable(append(no_magnet_filename, ' Processed.xls'), 'Sheet', 'Averages');

writetable(no_magent_average, earth_field_filename);
earth_field = readtable(earth_field_filename);

center_magnet_dict = sensorDataRead(center_magnet_filename);


% sensor_tables = cell(1, length(centered_magnet_sensor_list)); %pre-setting tables for each sensor

sensor_names = keys(center_magnet_dict);
sensor_data = values(center_magnet_dict);

%% Writing unbias sensing
for sensor_step = 1:length(keys(center_magnet_dict))
    
    current_sensor = sensor_data{sensor_step};
    current_name = sensor_names{sensor_step};
    num_entries = length(sensor_data);
    avg_XYZ = table2array(earth_field(sensor_step, 2:end));
        

    UnbiasX = current_sensor.X - avg_XYZ(1);
    UnbiasY = current_sensor.Y - avg_XYZ(2);
    UnbiasZ = current_sensor.Z - avg_XYZ(3);

    UnbiasTable = table(current_sensor.Time, UnbiasX, UnbiasY, UnbiasZ, 'VariableNames', {'Time', 'X','Y','Z'});
    unbias_sensing_list{sensor_step} = UnbiasTable;
    
end

for j=1:length(keys(center_magnet_dict))
  writetable(unbias_sensing_list{j},unbias_filename,'Sheet', sensor_names{j});
end
Unbias_sensor_dict = dictionary(sensor_names',unbias_sensing_list);
[avg_Unbiased,StD_Unbiased] =  average_XYZ(Unbias_sensor_dict);
writetable(avg_Unbiased,unbias_filename,'Sheet', 'Averages');
writetable(StD_Unbiased,unbias_filename,'Sheet', 'StDeviation');



%% plotting all 16 sensors without earths magnetic field
unbias_names = keys(Unbias_sensor_dict);
unbias_data = values(Unbias_sensor_dict);
for sensor_count = 1:size(keys(Unbias_sensor_dict))
    plot_dimension = sqrt(size(unbias_names));
    current_table = unbias_data{sensor_count};
    
    time_step = current_table.Time;
    X_graph = current_table.X;
    Y_graph = current_table.Y;
    Z_graph = current_table.Z;
    
%     figure(i)
    subplot(4,4,sensor_count)
    plot(time_step,X_graph); hold on
    plot(time_step,Y_graph)
    plot(time_step,Z_graph)
    title(sensor_names(sensor_count))
%     ylim([-0.1 0.1])
%     legend('X','Y','Z')
    
end 
legend('X','Y','Z')
sgtitle('UnbiasDataCollection');
saveas(gcf,'UnbiasDataCollection.png');
