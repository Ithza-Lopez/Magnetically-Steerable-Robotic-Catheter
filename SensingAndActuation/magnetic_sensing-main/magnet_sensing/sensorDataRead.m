function [sensor_dictionary] = sensorDataRead(filename)
data = readtable(filename);
time_index = data.TimeIndex;
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
    for list_item = 1:size(time_index)
        if sensor_step == (sensor_raw(list_item)+1) %save data per sensor 
            
            %Separates the X1,X2 into the two multiplexers
            X_m1(i,sensor_step) = X1_raw(list_item);
            Y_m1(i,sensor_step) = Y1_raw(list_item);
            Z_m1(i,sensor_step) = Z1_raw(list_item);
            time(i,sensor_step) = time_index(list_item);

            X_m2(i,sensor_step) = X2_raw(list_item);
            Y_m2(i,sensor_step) = Y2_raw(list_item);
            Z_m2(i,sensor_step) = Z2_raw(list_item);
            i = i+1;
        end
    end
    
    Multiplexer_1 = table(time(:,sensor_step),X_m1(:,sensor_step),Y_m1(:,sensor_step),Z_m1(:,sensor_step),'VariableNames', {'Time','X','Y','Z'});
    Multiplexer_1(Multiplexer_1.X == 0, :) = [];
    Multiplexer_2 = table(time(:,sensor_step), X_m2(:,sensor_step),Y_m2(:,sensor_step),Z_m2(:,sensor_step), 'VariableNames', {'Time','X','Y','Z'});
    Multiplexer_2(Multiplexer_2.X == 0, :) = [];

    
    sensor_list{sensor_step} = table({Multiplexer_1}, {Multiplexer_2}, 'VariableNames', {'Multiplexer_1', 'Multiplexer_2'});
end

%Naming each sensor to derived table
sensor0_top = sensor_list{1,1}.Multiplexer_1{1,1};
sensor1_top = sensor_list{1,2}.Multiplexer_1{1,1};
sensor2_top = sensor_list{1,3}.Multiplexer_1{1,1};
sensor3_top = sensor_list{1,4}.Multiplexer_1{1,1};
sensor4_left = sensor_list{1,5}.Multiplexer_1{1,1};
sensor5_left = sensor_list{1,6}.Multiplexer_1{1,1};
sensor6_left = sensor_list{1,7}.Multiplexer_1{1,1};
sensor7_left = sensor_list{1,8}.Multiplexer_1{1,1};

sensor0_right = sensor_list{1,1}.Multiplexer_2{1,1};
sensor1_right = sensor_list{1,2}.Multiplexer_2{1,1};
sensor2_right = sensor_list{1,3}.Multiplexer_2{1,1};
sensor3_right = sensor_list{1,4}.Multiplexer_2{1,1};
sensor4_bottom = sensor_list{1,5}.Multiplexer_2{1,1};
sensor5_bottom = sensor_list{1,6}.Multiplexer_2{1,1};
sensor6_bottom = sensor_list{1,7}.Multiplexer_2{1,1};
sensor7_bottom = sensor_list{1,8}.Multiplexer_2{1,1};

%order of sensor_list and sensor_names must match
sensor_list = {sensor0_top, sensor1_top, sensor2_top,sensor3_top,sensor4_left,...
    sensor5_left,sensor6_left, sensor7_left,sensor0_right,sensor1_right, sensor2_right,...
    sensor3_right, sensor4_bottom,sensor5_bottom, sensor6_bottom, sensor7_bottom};

sensor_names = ["sensor0_top", "sensor1_top", "sensor2_top","sensor3_top",...
    "sensor4_left", "sensor5_left","sensor6_left", "sensor7_left","sensor0_right",...
    "sensor1_right","sensor2_right", "sensor3_right", "sensor4_bottom",...
    "sensor5_bottom", "sensor6_bottom", "sensor7_bottom"]; 


table_name = append(filename,' Processed.xls');
sensor_dictionary = dictionary(sensor_names,sensor_list);

for j=1:length(keys(sensor_dictionary))
  writetable(sensor_list{j},table_name,'Sheet', sensor_names{j});
end

[average,st_dev] =  average_XYZ(sensor_dictionary);
writetable(average, table_name,'Sheet', 'Averages');
writetable(st_dev, table_name,'Sheet', 'StDeviation');


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