clear; clc; close all;

n = 4;
syms q1 q2 q3 q4 dq1 dq2 dq3 dq4 ddq1 ddq2 ddq3 ddq4 real
syms m1 m2 m3 m4 real
syms Ixx1 Iyy1 Izz1 Ixy1 Ixz1 Iyz1 real
syms Ixx2 Iyy2 Izz2 Ixy2 Ixz2 Iyz2 real
syms Ixx3 Iyy3 Izz3 Ixy3 Ixz3 Iyz3 real
syms Ixx4 Iyy4 Izz4 Ixy4 Ixz4 Iyz4 real
syms l3 l4 l5 real
syms cx1 cy1 cz1 cx2 cy2 cz2 cx3 cy3 cz3 cx4 cy4 cz4 real
syms g_gravity real

q = [q1; q2; q3; q4];
dq = [dq1; dq2; dq3; dq4];

I1 = [Ixx1 -Ixy1 -Ixz1; -Ixy1 Iyy1 -Iyz1; -Ixz1 -Iyz1 Izz1];
I2 = [Ixx2 -Ixy2 -Ixz2; -Ixy2 Iyy2 -Iyz2; -Ixz2 -Iyz2 Izz2];
I3 = [Ixx3 -Ixy3 -Ixz3; -Ixy3 Iyy3 -Iyz3; -Ixz3 -Iyz3 Izz3];
I4 = [Ixx4 -Ixy4 -Ixz4; -Ixy4 Iyy4 -Iyz4; -Ixz4 -Iyz4 Izz4];

T01 = [cos(0) -sin(0)*cos(-pi/2) sin(0)*sin(-pi/2) 0*cos(0); sin(0) cos(0)*cos(-pi/2) -cos(0)*sin(-pi/2) 0*sin(0); 0 sin(-pi/2) cos(-pi/2) q1; 0 0 0 1];
T12 = [cos(0) -sin(0)*cos(pi/2) sin(0)*sin(pi/2) 0*cos(0); sin(0) cos(0)*cos(pi/2) -cos(0)*sin(pi/2) 0*sin(0); 0 sin(pi/2) cos(pi/2) q2; 0 0 0 1];
T23 = [cos(q3) -sin(q3)*cos(-pi/2) sin(q3)*sin(-pi/2) l4*cos(q3); sin(q3) cos(q3)*cos(-pi/2) -cos(q3)*sin(-pi/2) l4*sin(q3); 0 sin(-pi/2) cos(-pi/2) l3; 0 0 0 1];
T34 = [cos(q4) -sin(q4)*cos(0) sin(q4)*sin(0) 0*cos(q4); sin(q4) cos(q4)*cos(0) -cos(q4)*sin(0) 0*sin(q4); 0 sin(0) cos(0) l5; 0 0 0 1];

T02 = simplify(T01 * T12);
T03 = simplify(T02 * T23);
T04 = simplify(T03 * T34);

pcm1 = T01 * [cx1; cy1; cz1; 1];
pcm2 = T02 * [cx2; cy2; cz2; 1];
pcm3 = T03 * [cx3; cy3; cz3; 1];
pcm4 = T04 * [cx4; cy4; cz4; 1];

X = [pcm1(1:3), pcm2(1:3), pcm3(1:3), pcm4(1:3)];

Jv1 = jacobian(X(:,1), q);
Jv2 = jacobian(X(:,2), q);
Jv3 = jacobian(X(:,3), q);
Jv4 = jacobian(X(:,4), q);

z0 = [0; 0; 1];
z1 = T01(1:3, 3);
z2 = T02(1:3, 3);
z3 = T03(1:3, 3);

Jw1 = zeros(3,4);
Jw2 = zeros(3,4);
Jw3 = [zeros(3,1), zeros(3,1), z2, zeros(3,1)]; 
Jw4 = [zeros(3,1), zeros(3,1), z2, z3];

R1 = T01(1:3,1:3);
R2 = T02(1:3,1:3);
R3 = T03(1:3,1:3);
R4 = T04(1:3,1:3);

M_sym = m1*(Jv1.'*Jv1) + Jw1.'*R1*I1*R1.'*Jw1 + ...
        m2*(Jv2.'*Jv2) + Jw2.'*R2*I2*R2.'*Jw2 + ...
        m3*(Jv3.'*Jv3) + Jw3.'*R3*I3*R3.'*Jw3 + ...
        m4*(Jv4.'*Jv4) + Jw4.'*R4*I4*R4.'*Jw4;
M_sym = simplify(M_sym);

C_sym = sym(zeros(n,n));
for i=1:n
    for j=1:n
        for k=1:n
            c_ijk = 0.5 * (diff(M_sym(i,j), q(k)) + diff(M_sym(i,k), q(j)) - diff(M_sym(j,k), q(i)));
            C_sym(i,j) = C_sym(i,j) + c_ijk*dq(k);
        end
    end
end
C_sym = simplify(C_sym);

P_pot = m1*g_gravity*pcm1(3) + m2*g_gravity*pcm2(3) + m3*g_gravity*pcm3(3) + m4*g_gravity*pcm4(3);
g_sym = jacobian(P_pot, q).';
g_sym = simplify(g_sym);

calc_M = matlabFunction(M_sym, 'Vars', {q1, q2, q3, q4, m1, m2, m3, m4, l3, l4, l5, cx1, cy1, cz1, cx2, cy2, cz2, cx3, cy3, cz3, cx4, cy4, cz4, Ixx1, Iyy1, Izz1, Ixy1, Ixz1, Iyz1, Ixx2, Iyy2, Izz2, Ixy2, Ixz2, Iyz2, Ixx3, Iyy3, Izz3, Ixy3, Ixz3, Iyz3, Ixx4, Iyy4, Izz4, Ixy4, Ixz4, Iyz4});
calc_C = matlabFunction(C_sym, 'Vars', {q1, q2, q3, q4, dq1, dq2, dq3, dq4, m1, m2, m3, m4, l3, l4, l5, cx1, cy1, cz1, cx2, cy2, cz2, cx3, cy3, cz3, cx4, cy4, cz4, Ixx1, Iyy1, Izz1, Ixy1, Ixz1, Iyz1, Ixx2, Iyy2, Izz2, Ixy2, Ixz2, Iyz2, Ixx3, Iyy3, Izz3, Ixy3, Ixz3, Iyz3, Ixx4, Iyy4, Izz4, Ixy4, Ixz4, Iyz4});
calc_g = matlabFunction(g_sym, 'Vars', {q1, q2, q3, q4, m1, m2, m3, m4, l3, l4, l5, cx1, cy1, cz1, cx2, cy2, cz2, cx3, cy3, cz3, cx4, cy4, cz4, g_gravity});

val_l3 = 0.21; val_l4 = 0.360; val_l5 = 0.1;
val_m = [6.4, 10.231, 5.141, 0.914];
val_g = 9.81;

val_I1 = [0.169, 0.043, 0.176, 0, 0, 0];
val_I2 = [0.288, 0.285, 0.033, 0, -0.009, 0];
val_I3 = [0.114, 0.013, 0.116, 0.010, 0, -0.007];
val_I4 = [0.001, 0.001, 0.001, 0, 0, 0];

rcm1_local = [0; 0; 0]; 
rcm2_local = [0; 0; 0];
rcm3_local = [-0.05; 0; 0];
rcm4_local = [0; 0; -0.05]; 

l1_inicial = 0; l1_final = 0.3;
l2_inicial = 0; l2_final = 0.5;
theta3_inicial = 0; theta3_final = pi/3;
theta4_inicial = 0; theta4_final = pi/6;
passos = 100;
t_final = 3;
t = linspace(0,t_final,passos);

al1 = calcular_coeficientes(l1_inicial,l1_final,t_final);
al2 = calcular_coeficientes(l2_inicial,l2_final,t_final);
atheta3 = calcular_coeficientes(theta3_inicial,theta3_final,t_final);
atheta4 = calcular_coeficientes(theta4_inicial,theta4_final,t_final);

p_traj = zeros(passos, 4);
v_traj = zeros(passos, 4);
ac_traj = zeros(passos, 4);
tau_results = zeros(passos, 4);

for i = 1:passos
    [pl1, vl1, acl1] = calcular_trajetoria(al1, t(i));
    [pl2, vl2, acl2] = calcular_trajetoria(al2, t(i));
    [ptheta3, vtheta3, actheta3] = calcular_trajetoria(atheta3, t(i));
    [ptheta4, vtheta4, actheta4] = calcular_trajetoria(atheta4, t(i));
    
    q_now = [max(0, pl1), max(0, pl2), ptheta3, ptheta4];
    dq_now = [vl1, vl2, vtheta3, vtheta4];
    ddq_now = [acl1, acl2, actheta3, actheta4];
    
    p_traj(i,:) = q_now;
    v_traj(i,:) = dq_now;
    ac_traj(i,:) = ddq_now;
    
    M_num = calc_M(q_now(1), q_now(2), q_now(3), q_now(4), ...
                   val_m(1), val_m(2), val_m(3), val_m(4), ...
                   val_l3, val_l4, val_l5, ...
                   rcm1_local(1), rcm1_local(2), rcm1_local(3), ...
                   rcm2_local(1), rcm2_local(2), rcm2_local(3), ...
                   rcm3_local(1), rcm3_local(2), rcm3_local(3), ...
                   rcm4_local(1), rcm4_local(2), rcm4_local(3), ...
                   val_I1(1), val_I1(2), val_I1(3), val_I1(4), val_I1(6), val_I1(5), ...
                   val_I2(1), val_I2(2), val_I2(3), val_I2(4), val_I2(6), val_I2(5), ...
                   val_I3(1), val_I3(2), val_I3(3), val_I3(4), val_I3(6), val_I3(5), ...
                   val_I4(1), val_I4(2), val_I4(3), val_I4(4), val_I4(6), val_I4(5));

    C_num = calc_C(q_now(1), q_now(2), q_now(3), q_now(4), ...
                   dq_now(1), dq_now(2), dq_now(3), dq_now(4), ...
                   val_m(1), val_m(2), val_m(3), val_m(4), ...
                   val_l3, val_l4, val_l5, ...
                   rcm1_local(1), rcm1_local(2), rcm1_local(3), ...
                   rcm2_local(1), rcm2_local(2), rcm2_local(3), ...
                   rcm3_local(1), rcm3_local(2), rcm3_local(3), ...
                   rcm4_local(1), rcm4_local(2), rcm4_local(3), ...
                   val_I1(1), val_I1(2), val_I1(3), val_I1(4), val_I1(6), val_I1(5), ...
                   val_I2(1), val_I2(2), val_I2(3), val_I2(4), val_I2(6), val_I2(5), ...
                   val_I3(1), val_I3(2), val_I3(3), val_I3(4), val_I3(6), val_I3(5), ...
                   val_I4(1), val_I4(2), val_I4(3), val_I4(4), val_I4(6), val_I4(5));
                   
    g_num = calc_g(q_now(1), q_now(2), q_now(3), q_now(4), ...
                   val_m(1), val_m(2), val_m(3), val_m(4), ...
                   val_l3, val_l4, val_l5, ...
                   rcm1_local(1), rcm1_local(2), rcm1_local(3), ...
                   rcm2_local(1), rcm2_local(2), rcm2_local(3), ...
                   rcm3_local(1), rcm3_local(2), rcm3_local(3), ...
                   rcm4_local(1), rcm4_local(2), rcm4_local(3), ...
                   val_g);
                   
    tau = M_num * ddq_now.' + C_num * dq_now.' + g_num;
    tau_results(i,:) = tau.';
end

figure;
titulos = {'Força F1 (N)', 'Força F2 (N)', 'Torque T3 (Nm)', 'Torque T4 (Nm)'};
for k = 1:4
    subplot(4,1,k);
    plot(t, tau_results(:,k), 'LineWidth', 1.5);
    title(titulos{k}); grid on;
    xlabel('Tempo (s)');
end

figure;
titulos_p = {'q_1 (m)', 'q_2 (m)', 'q_3 (rad)', 'q_4 (rad)'};
for k = 1:4
    subplot(4,1,k);
    plot(t, p_traj(:,k), 'LineWidth', 1.5);
    title(titulos_p{k}); grid on;
    xlabel('Tempo (s)');
end
titulos = {'v_1 (m)', 'v_2 (m)', 'v_3 (rad)', 'v_4 (rad)'};
figure;
for k = 1:4
    subplot(4,1,k);
    plot(t, v_traj(:,k), 'LineWidth', 1.5);
    title(titulos{k}); grid on;
    xlabel('Tempo (s)');
end
titulos = {'a_1 (m)', 'a_2 (m)', 'a_3 (rad)', 'a_4 (rad)'};
figure;
for k = 1:4
    subplot(4,1,k);
    plot(t, ac_traj(:,k), 'LineWidth', 1.5);
    title(titulos{k}); grid on;
    xlabel('Tempo (s)');
end



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