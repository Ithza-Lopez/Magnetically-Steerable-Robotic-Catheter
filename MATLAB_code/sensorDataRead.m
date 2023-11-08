function [sensor_dictionary] = sensorDataRead(filename)
data = readtable(filename);
% time_index = data.TimeIndex;
X1_raw = data.X1;
Y1_raw = data.Y1;
Z1_raw = data.Z1;
X2_raw = data.X2;
Y2_raw = data.Y2;
Z2_raw = data.Z2;
sensor_raw = data.Sensor;

num_sensors = unique(sensor_raw); %count numbe of sensors

%Loop to find data from each sensor 0-7
for sensor_step = 1:size(num_sensors)
    i=1;
    for list_item = 1:size(X1_raw)
        if sensor_step == (sensor_raw(list_item)+1) %save data per sensor 
            
            %Separates the X1,X2 into the two multiplexers
            X_m1(i,sensor_step) = X1_raw(list_item);
            Y_m1(i,sensor_step) = Y1_raw(list_item);
            Z_m1(i,sensor_step) = Z1_raw(list_item);
%             time(i,sensor_step) = time_index(list_item);

            X_m2(i,sensor_step) = X2_raw(list_item);
            Y_m2(i,sensor_step) = Y2_raw(list_item);
            Z_m2(i,sensor_step) = Z2_raw(list_item);
            i = i+1;
        end
    end
    
    Multiplexer_1 = table(X_m1(:,sensor_step),Y_m1(:,sensor_step),Z_m1(:,sensor_step),'VariableNames', {'X','Y','Z'});
    Multiplexer_1(Multiplexer_1.X == 0, :) = [];
    Multiplexer_2 = table(X_m2(:,sensor_step),Y_m2(:,sensor_step),Z_m2(:,sensor_step), 'VariableNames', {'X','Y','Z'});
    Multiplexer_2(Multiplexer_2.X == 0, :) = [];

    
    sensor_list{sensor_step} = table({Multiplexer_1}, {Multiplexer_2}, 'VariableNames', {'Multiplexer_1', 'Multiplexer_2'});
end

%Naming each sensor to derived table
sensor_5 = sensor_list{1,1}.Multiplexer_1{1,1};
sensor_6 = sensor_list{1,2}.Multiplexer_1{1,1};
sensor_7 = sensor_list{1,3}.Multiplexer_1{1,1};
sensor_8 = sensor_list{1,4}.Multiplexer_1{1,1};
sensor_9 = sensor_list{1,5}.Multiplexer_1{1,1};
sensor_10 = sensor_list{1,6}.Multiplexer_1{1,1};
sensor_11 = sensor_list{1,7}.Multiplexer_1{1,1};
sensor_12 = sensor_list{1,8}.Multiplexer_1{1,1};

sensor_1 = sensor_list{1,1}.Multiplexer_2{1,1};
sensor_2 = sensor_list{1,2}.Multiplexer_2{1,1};
sensor_3 = sensor_list{1,3}.Multiplexer_2{1,1};
sensor_4 = sensor_list{1,4}.Multiplexer_2{1,1};
sensor_13 = sensor_list{1,5}.Multiplexer_2{1,1};
sensor_14 = sensor_list{1,6}.Multiplexer_2{1,1};
sensor_15 = sensor_list{1,7}.Multiplexer_2{1,1};
sensor_16 = sensor_list{1,8}.Multiplexer_2{1,1};

%order of sensor_list and sensor_names must match
sensor_list = {sensor_1,sensor_2, sensor_3,...
    sensor_4, sensor_5, sensor_6, sensor_7,sensor_8,sensor_9,...
    sensor_10,sensor_11, sensor_12, sensor_13,sensor_14, sensor_15, sensor_16};

sensor_names = ["sensor_1","sensor_2","sensor_3", "sensor_4",...
    "sensor_5", "sensor_6", "sensor_7","sensor_8",...
    "sensor_9", "sensor_10","sensor_11", "sensor_12", "sensor_13",...
    "sensor_14", "sensor_15", "sensor_16"]; 


% table_name = append(filename,' Processed.xls');
sensor_dictionary = dictionary(sensor_names,sensor_list);

for j=1:length(keys(sensor_dictionary))
%   writetable(sensor_list{j},filename,'Sheet', sensor_names{j});
    writetable(sensor_list{j},filename,'Sheet', sensor_names{j});
end

[average,st_dev] =  average_XYZ(sensor_dictionary);
writetable(average, filename,'Sheet', 'Averages');
writetable(st_dev, filename,'Sheet', 'StDeviation');


% %plotting data from sensor
% for i = 1:length(sensor_list)
%     
%     plot_dimension = sqrt(length(sensor_list));
%     current_table = sensor_list{i};
%     
%     time_step = current_table.Time;
%     X_graph = current_table.X;
%     Y_graph = current_table.Y;
%     Z_graph = current_table.Z;
%     
% %     figure(i)
%     subplot(plot_dimension,plot_dimension,i)
%     plot(time_step,X_graph); hold on
%     plot(time_step,Y_graph)
%     plot(time_step,Z_graph)
%     title(sensor_names(i))
%     ylim([-2.5 2.5])
% %     legend('X','Y','Z')
%     
% end 
% legend('X','Y','Z')
% figure_name = append(filename,'_plot.png');
% sgtitle(filename,'Interpreter', 'none');
% saveas(gcf,figure_name);
% close all

end