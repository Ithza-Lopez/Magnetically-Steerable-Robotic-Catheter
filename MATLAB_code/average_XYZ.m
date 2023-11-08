
function [avg_MagField,StD_MagField] =  average_XYZ(sensor_dict)

sensor_names = keys(sensor_dict);
sensor_data = values(sensor_dict);

for i = 1:length(sensor_names)
    
    current_table = sensor_data{i};
    XYZ_means = mean(current_table{:,["X", "Y","Z"]});
    XYZ_std = std(current_table{:,["X", "Y","Z"]});
    sensor_XYZ_avg(i,:) = XYZ_means;
    sensor_XYZ_std(i,:) = XYZ_std;
end
% disp(sensor_XYZ_avg)
% disp(sensor_names)
avg_MagField = table(sensor_names(:),sensor_XYZ_avg(:,1), sensor_XYZ_avg(:,2),sensor_XYZ_avg(:,3), 'VariableNames', {'Sensor','X_avg','Y_avg','Z_avg'});
StD_MagField = table(sensor_names(:),sensor_XYZ_std(:,1), sensor_XYZ_std(:,2),sensor_XYZ_std(:,3), 'VariableNames', {'Sensor','X_StD','Y_StD','Z_StD'});

end 