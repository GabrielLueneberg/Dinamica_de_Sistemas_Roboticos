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

x_inicial = 0.6;
y_inicial = 1;
z_inicial = 1;
nz_inicial = -1;

x_final = 0.6;
y_final = 0.1;
z_final = 0.1;
nz_final = 0;


passos = 100;
t_final = 5;
t = linspace(0,t_final,passos);

function a = calcular_coeficientes(p_i, p_f, t)
    a0 = p_i;
    a1 = 0;
    a2 = 0;
    a3 = (10 * (p_f - p_i)) / (t^3);
    a4 = (-15 * (p_f - p_i)) / (t^4);
    a5 = (6 * (p_f - p_i)) / (t^5); 
    
    a = [a0, a1, a2, a3, a4, a5];

end

function [p, v, ac] = calcular_trajetoria(a,t)
    p = a(1) + a(2)*t + a(3)*t^2 + a(4)*t^3 + a(5)*t^4 + a(6)*t^5;
    v = a(2) + 2*a(3)*t + 3*a(4)*t^2 + 4*a(5)*t^3 + 5*a(6)*t^4;
    ac = 2*a(3) + 6*a(4)*t + 12*a(5)*t^2 + 20*a(6)*t^3;
end

ax = calcular_coeficientes(x_inicial,x_final,t_final);
ay = calcular_coeficientes(y_inicial,y_final,t_final);
az = calcular_coeficientes(z_inicial,z_final,t_final);
anz =calcular_coeficientes(nz_inicial,nz_final,t_final);

q_traj = zeros(passos, 4);

for i = 1:passos

    [px, vx, acx] = calcular_trajetoria(ax, t(i));
    [py, vy, acy] = calcular_trajetoria(ay, t(i));
    [pz, vz, acz] = calcular_trajetoria(az, t(i));
    [nz, nvz, nacz] = calcular_trajetoria(anz, t(i));

    l1 = pz - l3;
    arg_acos = px/(sqrt(l4^2+l5^2));
    arg_acos = max(-1, min(1, arg_acos));
    theta3 = acos(arg_acos) - atan(l5/l4);
    l2 = py - l4*sin(theta3) - l5*cos(theta3);
    theta4 = asin(-nz);
    q_traj(i,:) = [l1,l2,theta3,theta4];
   
end    
q_traj(:, 1:2) = max(0.01, q_traj(:, 1:2));
robot.plot(q_traj, 'fps', 24, 'trail', 'r-');
