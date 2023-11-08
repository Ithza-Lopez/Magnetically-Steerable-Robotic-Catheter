function PM_Est_Plot()
close all


figure(1)
%plotting data
plot(pm_true_all(1,:)*1e3,pm_true_all(2,:)*1e3,'bo')
hold on
plot(pm_est_all(1,:)*1e3,pm_est_all(2,:)*1e3,'rx')
xlabel('x[mm]')
ylabel('y[mm]')
legend({'True Position', 'Estimated Position'}, "Location",'bestoutside')
title('Top View')
hold on

%plotting border
g_x = 31.75/2;
g_y = 178.68/2;
% grid_x = [-g_x, g_x, g_x, -g_x, -g_x];
% grid_y = [-g_y, -g_y, g_y, g_y, -g_y];
% % plot(grid_x, grid_y, 'LineWidth',2, 'Color','k')
rectangle('Position',[-g_x,-g_y,g_x*2, g_y*2] , 'LineWidth',1, 'EdgeColor','k')
hold on

%fixing aspect ratio
daspect([1,1,1])
axis([-30, 30, -100, 100])

%adding grey out boxes
box1_x = 31.75/2;
box2_y = 178.68/2;
box_y_edge = 45;
rectangle('Position',[-box1_x,-box2_y,2*box1_x, box_y_edge] , 'Facecolor', [0.7, 0.7, 0.7, 0.5])
rectangle('Position',[-box1_x,box_y_edge,2*box1_x, box2_y/2] , 'Facecolor', [0.7, 0.7, 0.7, 0.5])
hold on

saveas(figure(1) , 'Full Top View.pdf' )

%% Plot only middle section
figure(2)
slots_per_column = length(-8:1:8);
index_check = [-8:1:8,-8:1:8,-8:1:8];
row_range = -4:1:4;
num_rows = length(row_range);

index_A = 1 + max(row_range) :1: max(row_range)+num_rows;
index_B = 1 + slots_per_column + max(row_range) :1: slots_per_column + max(row_range) + num_rows;
index_C = 1 + slots_per_column*2 + max(row_range) :1: slots_per_column*2 + max(row_range) + num_rows;

%plotting data
plot(pm_true_all(1,[index_A,index_B, index_C])*1e3,pm_true_all(2,[index_A,index_B, index_C])*1e3,'bo')
hold on
plot(pm_est_all(1,[index_A,index_B, index_C])*1e3,pm_est_all(2,[index_A,index_B, index_C])*1e3,'rx')

xlabel('x[mm]')
ylabel('y[mm]')
legend({'True Position', 'Estimated Position'}, "Location",'bestoutside')
title('Top View')
hold on

%plotting grid border
g_x = 31.75/2;
g_y = (max(row_range)+1)/2 * 10 +20;
line([-g_x -g_x], [-g_y g_y], 'LineWidth',1, 'LineStyle','-', 'Color','k','HandleVisibility','off')
line([g_x g_x], [-g_y g_y], 'LineWidth',1, 'LineStyle','-', 'Color','k','HandleVisibility','off')
line([-g_x g_x], [g_y g_y], 'LineWidth',1, 'LineStyle','--', 'Color','k','HandleVisibility','off')
line([-g_x g_x], [-g_y -g_y], 'LineWidth',1, 'LineStyle','--', 'Color','k','HandleVisibility','off')
hold on

%fixing aspect ratio
margin_offset = 10;
daspect([1,1,1])
axis([-g_x - margin_offset, g_x + margin_offset, -g_y - margin_offset, g_y + margin_offset])

saveas(figure(2) , 'Row Range Top View.pdf' )