function XYZ_estimation_plots(measured, theoretical);
X_measured = measured(1:3:end);
Y_measured = measured(2:3:end);
Z_measured = measured(3:3:end);

X_theoretical = theoretical(1:3:end);
Y_theoretical = theoretical(2:3:end);
Z_theoretical = theoretical(3:3:end);

x_slope = min(X_theoretical):0.000001:max(X_theoretical);
y_slope = min(Y_theoretical):0.000001:max(Y_theoretical);
z_slope = min(Z_measured):0.000001:max(Z_measured);


figure(1)
subplot(1,3,1)
plot(X_measured, X_theoretical, 'ro')
hold on
plot(x_slope,x_slope,"r-")
xlabel('X')
ylabel('$\hat{X}$', 'Interpreter', 'latex');

subplot(1,3,2)
plot(Y_measured, Y_theoretical, 'bo')
hold on
plot(y_slope,y_slope,"b-")
xlabel('Y')
ylabel('$\hat{Y}$', 'Interpreter', 'latex');

subplot(1,3,3)
plot(Z_measured, Z_theoretical, "Color",[0.4660 0.6740 0.1880] ,'Marker','o', 'LineStyle','none')
hold on
plot(z_slope,z_slope,"Color",[0.4660 0.6740 0.1880],"LineStyle","-")
xlabel('Z')
ylabel('$\hat{Z}$', 'Interpreter', 'latex');

end