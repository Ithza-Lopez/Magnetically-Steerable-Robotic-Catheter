function [Bx,By,Bz] = loop_wire_mag(R,pos_eva,I)
% input:
% R: radius of the current loop
% I: current value
% pos_eva: the position vector of the point to be evaluated (x,y,z)

% output:
% Bx,By,Bz: three axis magnetic field components

% configuration: 
% frame origin at center of the current loop
% z axis pointing perpendicular to the loop 

mu_0 = 4*pi*1e-7;
x_p = pos_eva(1);
y_p = pos_eva(2);
z_p = pos_eva(3);


% use intergral
funx = @(theta,x,y,z,R) z*R*cos(theta)./(x^2+y^2+z^2+R^2-2*x*R.*cos(theta)-2*x*R.*sin(theta)).^(3/2);
Bx = mu_0*I/4/pi*integral(@(x) funx(x,x_p,y_p,z_p,R),0,2*pi);

funy = @(theta,x,y,z,R) z*R*cos(theta)./(x^2+y^2+z^2+R^2-2*x*R.*cos(theta)-2*x*R.*sin(theta)).^(3/2);
By = mu_0*I/4/pi*integral(@(x) funy(x,x_p,y_p,z_p,R),0,2*pi);



end