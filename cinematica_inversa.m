clear; clc; close all;

l3= 0.21;
l4 = 0.360;
l5 = 0.03774;

L(1) = Link('theta', 0, 'a', 0, 'alpha', -pi/2, 'prismatic');
L(2) = Link('theta', 0, 'a', 0, 'alpha', pi/2, 'prismatic');
L(3) = Link('d', l3, 'a', l4, 'alpha', -pi/2);
L(4) = Link('d', l5, 'a', 0, 'alpha', 0); 

L(1).qlim = [0 0.470];
L(2).qlim = [0 0.5];

robot = SerialLink(L, 'name', 'links_manipulador');

user_input = input('Matriz de posições: '); % digite a matriz de posições
% ex: [0 0 -1 -0.1; 0 -1 0 0.86; -1 0 0 0.47; 0 0 0 1]

px = user_input(1,4);
py = user_input(2,4);
pz = user_input(3,4);
nz = user_input(3,1);

l1 = pz - l3
theta3 = acos(px/(sqrt(l4^2+l5^2))) - atan(l5/l4)
l2 = py - l4*sin(theta3) - l5*cos(theta3)
theta4 = asin(-nz)

q = [l1, l2, theta3, theta4];
robot.plot(q);



