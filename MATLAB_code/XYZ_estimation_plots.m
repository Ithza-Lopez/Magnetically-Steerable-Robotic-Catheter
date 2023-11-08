function XYZ_estimation_plots(measured, theoretical)
X_measured = measured(1,:);
Y_measured = measured(2,:);
Z_measured = measured(3,:);

X_theoretical = theoretical(1,:);
Y_theoretical = theoretical(2,:);
Z_theoretical = theoretical(3,:);

step = 0.001;
x_slope_t = min(X_theoretical): step :max(X_theoretical);
y_slope_t = min(Y_theoretical): step :max(Y_theoretical);
z_slope_t = min(Z_theoretical): step :max(Z_theoretical);

x_slope_m = zeros(length(x_slope_t));
y_slope_m = min(Y_measured): step :max(Y_measured);
z_slope_m = zeros(length(z_slope_t));

figure()
subplot(1,3,1)
plot(X_measured, X_theoretical, 'ro')
hold on
plot(x_slope_m,x_slope_t,"r-")
xlabel('X')
ylabel('$\hat{X}$', 'Interpreter', 'latex');

subplot(1,3,2)
plot(Y_measured, Y_theoretical, 'bo')
hold on
plot(y_slope_m,y_slope_m,"b-")
xlabel('Y')
ylabel('$\hat{Y}$', 'Interpreter', 'latex');

subplot(1,3,3)
plot(Z_measured, Z_theoretical,'go')
hold on
plot(z_slope_m,z_slope_t,"g-")
xlabel('Z')
ylabel('$\hat{Z}$', 'Interpreter', 'latex');
hold on

saveas(figure(1), 'XYZ_estimation_plots.jpeg' )
end