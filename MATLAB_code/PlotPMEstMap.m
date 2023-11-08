function PlotPMEstMap(pm_true_all, pm_est_all )

figure(1)
plot3(pm_true_all(1,:)*1e3,pm_true_all(2,:)*1e3,pm_true_all(3,:)*1e3,'b-o')
hold on
plot3(pm_est_all(1,:)*1e3,pm_est_all(2,:)*1e3, pm_est_all(3,:)*1e3,'rx')
xlabel('x[mm]')
ylabel('y[mm]')
zlabel('z[mm]')
legend('True','Estimation')
title('XYZ View')
% xlim([-2 2])
% ylim([-90 90])
% zlim([-2 2])

figure(2)
plot(pm_true_all(1,:)*1e3,pm_true_all(2,:)*1e3,'bo')
hold on
plot(pm_est_all(1,:)*1e3,pm_est_all(2,:)*1e3,'rx')
xlabel('x[mm]')
ylabel('y[mm]')
legend('True','Estimation')
title('Top View')
% xlim([-2 2])
% ylim([-90 90])


figure(3)
plot(pm_true_all(2,:)*1e3,pm_true_all(3,:)*1e3,'bo')
hold on
plot(pm_est_all(2,:)*1e3,pm_est_all(3,:)*1e3,'rx')
xlabel('y[mm]')
ylabel('z[mm]')
legend('True','Estimation')
title('Side View')
% xlim([-90 90])
% ylim([-2.5 2.5])

saveas(figure(1), 'XYZ View.jpeg' )
saveas(figure(2), 'Top View.jpeg' )
saveas(figure(3), 'Side View.jpeg' )
end