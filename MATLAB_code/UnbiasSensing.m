function UnbiasSensing(no_magnet_filename,magnet_location_filename,earth_field_filename, unbias_filename)
%this filters the earth field from the raw data
%%initializing data from files
no_magnet_dict = sensorDataRead(no_magnet_filename);
no_magent_average = readtable(no_magnet_filename, 'Sheet', 'Averages');

writetable(no_magent_average, earth_field_filename);
earth_field = readtable(earth_field_filename);

center_magnet_dict = sensorDataRead(magnet_location_filename);


% sensor_tables = cell(1, length(centered_magnet_sensor_list)); %pre-setting tables for each sensor

sensor_names = keys(center_magnet_dict);
sensor_data = values(center_magnet_dict);

% Writing unbias sensing
for sensor_step = 1:length(keys(center_magnet_dict))
    
    current_sensor = sensor_data{sensor_step};
    current_name = sensor_names{sensor_step};
    num_entries = length(sensor_data);
    avg_earthfield = table2array(earth_field(sensor_step, 2:end));
        

    UnbiasX = current_sensor.X - avg_earthfield(1);
    UnbiasY = current_sensor.Y - avg_earthfield(2);
    UnbiasZ = current_sensor.Z - avg_earthfield(3);

    UnbiasTable = table( UnbiasX, UnbiasY, UnbiasZ, 'VariableNames', {'X','Y','Z'});
    unbias_sensing_list{sensor_step} = UnbiasTable;
    
end

for j=1:length(keys(center_magnet_dict))
  writetable(unbias_sensing_list{j},unbias_filename,'Sheet', sensor_names{j});
end

Unbias_sensor_dict = dictionary(sensor_names',unbias_sensing_list);
[avg_Unbiased,StD_Unbiased] =  average_XYZ(Unbias_sensor_dict);
writetable(avg_Unbiased,unbias_filename,'Sheet', 'Averages');
writetable(StD_Unbiased,unbias_filename,'Sheet', 'StDeviation');

avg_vector = reshape((table2array(avg_Unbiased(:,2:end))).',1,[]);
writematrix(avg_vector, unbias_filename,'Sheet', 'Avg_vector');

end