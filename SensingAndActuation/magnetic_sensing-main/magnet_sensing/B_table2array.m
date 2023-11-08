function [B_mesured_array] = B_array(filename)

file_format = append(filename, '_by_sensor.xls');
averages_table = readtable(file_format, 'Sheet','Averages');
average_array = table2array(averages_table);
B_mesured_array = reshape(average_array',[1 size(average_array,1)*size(average_array,2)]);

end