clear; clc; close all;

l3 = 0.21;
l4 = 0.360;
l5 = 0.1;

L(1) = Link('theta', 0, 'a', 0, 'alpha', -pi/2, 'prismatic');
L(1).m = 6.4;
L(1).I = [0.169, 0.043, 0.176, 0, 0, 0];%kg*m^2
%vetor posição da origem até o centro de massa x=90mm,y=-10.04mm,z=74.73mm

L(2) = Link('theta', 0, 'a', 0, 'alpha', pi/2, 'prismatic');
L(2).m = 10.231;
L(2).I = [0.288, 0.285, 0.033, 0, -0.009, 0];%kg*m^2
%vetor posição da origem até o centro de massa x=95,19mm,y=-10mm,z=319.38mm

L(3) = Link('d', l3, 'a', l4, 'alpha', -pi/2);
L(3).m = 5.141;
L(3).I = [0.114, 0.013, 0.116, 0.010, 0, -0.007];
%vetor posição da origem até o centro de massa x=293.22mm,y=-21.84mm,z=290.61mm

L(4) = Link('d', l5, 'a', 0, 'alpha', 0); 
L(4).m = 0.914;
L(4).I = [0.001, 0.001, 0.001, 0, 0, 0];
%vetor posição da origem até o centro de massa x=300.02mm,y=-147.79mm,z=106.48mm

L(1).qlim = [0 0.470];
L(2).qlim = [0 0.5];

robot = SerialLink(L, 'name', 'links_manipulador');

l1_inicial = 0;
l1_final = 0.3;
l2_inicial = 0;
l2_final = 0.5;
theta3_inicial = 0;
theta3_final = pi/3;
theta4_inicial = 0;
theta4_final = pi/6;

passos = 100;
t_final = 3;
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

al1 = calcular_coeficientes(l1_inicial,l1_final,t_final);
al2 = calcular_coeficientes(l2_inicial,l2_final,t_final);
atheta3 = calcular_coeficientes(theta3_inicial,theta3_final,t_final);
atheta4 =calcular_coeficientes(theta4_inicial,theta4_final,t_final);

q_traj = zeros(passos, 4);

for i = 1:passos

    [pl1, vl1, acl1] = calcular_trajetoria(al1, t(i));
    [pl2, vl2, acl2] = calcular_trajetoria(al2, t(i));
    [ptheta3, vtheta3, actheta3] = calcular_trajetoria(atheta3, t(i));
    [ptheta4, vtheta4, actheta4] = calcular_trajetoria(atheta4, t(i));

    p_traj(i,:) = [pl1,pl2,ptheta3,ptheta4];
    v_traj(i,:) = [vl1,vl2,vtheta3,vtheta4];
    ac_traj(i,:) = [acl1,acl2,actheta3,actheta4];
   
end    
p_traj(:, 1:2) = max(0, p_traj(:, 1:2));

r_cm1 = [90, -10.04, 74.73] / 1000;
r_cm2 = [95.19, -10, 319.38] / 1000;
r_cm3 = [293.22, -21.84, 290.61] / 1000;
r_cm4 = [300.02, -147.79, 106.48] / 1000;

lcm1 = norm(r_cm1);
lcm2 = norm(r_cm2);
lcm3 = norm(r_cm3 - r_cm2);
lcm4 = norm(r_cm4 - (r_cm2+[l4,0,0]));

m_vals = [L(1).m, L(2).m, L(3).m, L(4).m];
geom_vals = [lcm1, lcm2, lcm3, lcm4, l3, l4, l5];

I_vals = [];
for k = 1:4
    I = L(k).I; 
    I_k = [I(1), I(2), I(3), I(4), I(6), I(5)]; 
    I_vals = [I_vals, I_k];
end

tau_results = zeros(passos, 4);

for i = 1:passos
    q_i   = p_traj(i,:);
    dq_i  = v_traj(i,:);
    ddq_i = ac_traj(i,:);
    
    entrada = [q_i, dq_i, ddq_i, m_vals, geom_vals, I_vals];
    
    tau_results(i,:) = TauFun(entrada)';
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
titulos = {'q_1 (m)', 'q_2 (m)', 'q_3 (rad)', 'q_4 (rad)'};
for k = 1:4
    subplot(4,1,k);
    plot(t, p_traj(:,k), 'LineWidth', 1.5);
    title(titulos{k}); grid on;
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