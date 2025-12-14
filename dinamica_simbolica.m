syms t l1 l2 theta3 theta4 dl1 dl2 dt3 dt4 ddl1 ddl2 ddt3 ddt4 real
syms m1 m2 m3 m4 real
syms Ixx1 Iyy1 Izz1 Ixy1 Ixz1 Iyz1 Ixx2 Iyy2 Izz2 Ixy2 Ixz2 Iyz2 real
syms Ixx3 Iyy3 Izz3 Ixy3 Ixz3 Iyz3 Ixx4 Iyy4 Izz4 Ixy4 Ixz4 Iyz4 real
syms lcm1 lcm2 lcm3 lcm4 l3 l4 l5 real

c3 = cos(theta3);
s3 = sin(theta3);
c4 = cos(theta4);
s4 = sin(theta4);

A1 = sym(zeros(4,4));
A2 = sym((zeros(4,4)));
A3 = sym(zeros(4,4));
A02 = sym(zeros(4,4));
A03 = sym(zeros(4,4));
A04 = sym(zeros(4,4));
A05 = sym(zeros(4,4));

A01 = [1,0,0,0;
      0,0,1,0;
      0,-1,0,l1;
      0,0,0,1];

A2 = [1,0,0,0;
      0,0,-1,0;
      0,1,0,l2;
      0,0,0,1];

A3 = [c3,0,-s3,l4*c3;
      s3,0,c3,l4*s3;
      0,-1,0,l3;
      0,0,0,1];

A4 = [c4,-s4,0,0;
      s4,c4,0,0;
      0,0,1,l5;
      0,0,0,1];

A02 = A01*A2;
A03 = A02*A3;
A04 = A03*A4;


b0 = sym(zeros(3,1));
b1 = sym(zeros(3,1));
b2 = sym(zeros(3,1));
b3 = sym(zeros(3,1));
b4 = sym(zeros(3,1));

b0 = [0;
      0;
      1];

b1 = A01(1:3,3);
b2 = A02(1:3,3);
b3 = A03(1:3,3);
b4 = A04(1:3,3);

R1 = A01(1:3,1:3);
R2 = A02(1:3,1:3);
R3 = A03(1:3,1:3);
R4 = A04(1:3,1:3);

%%Juntas prismáticas
%%Link 1
JP11 = sym(zeros(3,1));
JP12 = sym(zeros(3,1));
JP13 = sym(zeros(3,1));
JP14 = sym(zeros(3,1));
JR1 = sym(zeros(3,4)); 
JP1 = sym(zeros(3,4));


JP11 = b0;

JP1(:,1) = JP11;
JP1(:,2) = JP12;
JP1(:,3) = JP13;
JP1(:,4) = JP14;

%%Link 2
JP21 = sym(zeros(3,1));
JP22 = sym(zeros(3,1));
JP23 = sym(zeros(3,1));
JP24 = sym(zeros(3,1));
JR2 = sym(zeros(3,4)); 
JP2 = sym(zeros(3,4));

JP21 = b0;
JP22 = b1;

JP2(:,1) = JP21;
JP2(:,2) = JP22;
JP2(:,3) = JP23;
JP2(:,4) = JP24;

%%Juntas rotativas
O0 = sym(zeros(3,1));
O1 = sym(zeros(3,1));
O2 = sym(zeros(3,1));
O3 = sym(zeros(3,1));
O4 = sym(zeros(3,1));

pcm1 = sym(zeros(3,1));
pcm2 = sym(zeros(3,1));
pcm3 = sym(zeros(3,1));
pcm4 = sym(zeros(3,1));

O2 = [ 0;
       0;
       l1];
O3 = [0;
      l2;
      l1+l3];
O4 = [l4;
      l2;
      l1+l3];

pcm1 = R1*[0, 0, lcm1]';

pcm2 = O2 + R2*[0, lcm2, 0]';

pcm3 = O3 + R3*[lcm3, 0, 0]';

pcm4 = O4 + R4*[0, lcm4, 0]';


%%Link 3
JP31 = sym(zeros(3,1));
JP32 = sym(zeros(3,1));
JP33 = sym(zeros(3,1));
JP34 = sym(zeros(3,1));
JR31 = sym(zeros(3,1));
JR32 = sym(zeros(3,1));
JR33 = sym(zeros(3,1));
JR34 = sym(zeros(3,1));
JP3 = sym(zeros(3,4));
JR3 = sym(zeros(3,4));

JP31 = cross(b0,pcm3 - O0);
JP32 = cross(b1,pcm3 - O1);
JP33 = cross(b2,pcm3 - O2);

JP3(:,1) = JP31;
JP3(:,2) = JP32;
JP3(:,3) = JP33;
JP3(:,4) = JP34;

JR31 = b0;
JR32 = b1;
JR33 = b2;

JR3(:,1) = JR31;
JR3(:,2) = JR32;
JR3(:,3) = JR33;
JR3(:,4) = JR34;


%%Link 4
JP41 = sym(zeros(3,1));
JP42 = sym(zeros(3,1));
JP43 = sym(zeros(3,1));
JP44 = sym(zeros(3,1));
JR41 = sym(zeros(3,1));
JR42 = sym(zeros(3,1));
JR43 = sym(zeros(3,1));
JR44 = sym(zeros(3,1));
JP4 = sym(zeros(3,4));
JR4 = sym(zeros(3,4));

JP41 = cross(b0,pcm4 - O0);
JP42 = cross(b1,pcm4 - O1);
JP43 = cross(b2,pcm4 - O2);
JP44 = cross(b3,pcm4 - O2);

JP4(:,1) = JP41;
JP4(:,2) = JP42;
JP4(:,3) = JP43;
JP4(:,4) = JP44;

JR41 = b0;
JR42 = b1;
JR43 = b2;
JR44 = b3;

JR4(:,1) = JR41;
JR4(:,2) = JR42;
JR4(:,3) = JR43;
JR4(:,4) = JR44;

%%Matriz de Massa
R1 = A01(1:3,1:3);
R2 = A02(1:3,1:3);
R3 = A03(1:3,1:3);
R4 = A04(1:3,1:3);

I1 = sym(zeros(3,3));
I2 = sym(zeros(3,3));
I3 = sym(zeros(3,3));
M1 = sym(zeros(3,4));
M2 = sym(zeros(3,4));
M3 = sym(zeros(3,4));
M4 = sym(zeros(3,4));

I1=[Ixx1, Ixy1, Ixz1;
    Ixy1, Iyy1, Iyz1;
    Ixz1, Iyz1, Izz1 ];

I2=[Ixx2, Ixy2, Ixz2;
    Ixy2, Iyy2, Iyz2;
    Ixz2, Iyz2, Izz2 ];

I3=[Ixx3, Ixy3, Ixz3;
    Ixy3, Iyy3, Iyz3;
    Ixz3, Iyz3, Izz3 ];

I4=[Ixx4, Ixy4, Ixz4;
    Ixy4, Iyy4, Iyz4;
    Ixz4, Iyz4, Izz4 ];

M1 = m1 * JP1.' * JP1;
M2 = m2 * JP2.' * JP2;
M3 = (m3 * JP3' * JP3) + (JR3.'*R3*I3*R3.'*JR3);
M4 = (m4 * JP4.' * JP4) + (JR4.'*R4*I4*R4.'*JR4);
M = M1+M2+M3+M4;

%%Lagrange
q  = [l1; l2; theta3; theta4];
dq = [dl1; dl2; dt3; dt4];
ddq = [ddl1; ddl2; ddt3; ddt4];

g = [0; -9.81; 0];

V1 = m1 * g.' * pcm1;
V2 = m2 * g.' * pcm2;
V3 = m3 * g.' * pcm3;
V4 = m4 * g.' * pcm4;
V_total = V1 + V2 + V3 + V4;

K_total = 0.5 * dq.' * M * dq;

L = K_total - V_total;

dL_ddq = jacobian(L, dq).'; 
dL_dq  = jacobian(L, q).';  

ddt_dL_ddq = jacobian(dL_ddq, q)*dq + jacobian(dL_ddq, dq)*ddq;

Tau = ddt_dL_ddq - dL_dq;

Tau = simplify(Tau);

F1 = Tau(1);
F2 = Tau(2);
T3 = Tau(3);
T4 = Tau(4);

% Vetor de esforços
Tau = [F1; F2; T3; T4];

params_inercia = [Ixx1, Iyy1, Izz1, Ixy1, Ixz1, Iyz1, ...
                  Ixx2, Iyy2, Izz2, Ixy2, Ixz2, Iyz2, ...
                  Ixx3, Iyy3, Izz3, Ixy3, Ixz3, Iyz3, ...
                  Ixx4, Iyy4, Izz4, Ixy4, Ixz4, Iyz4];

params_massa = [m1, m2, m3, m4];
params_geom  = [lcm1, lcm2, lcm3, lcm4, l3, l4, l5];

vars_q   = [l1, l2, theta3, theta4];
vars_dq  = [dl1, dl2, dt3, dt4];
vars_ddq = [ddl1, ddl2, ddt3, ddt4];

input_vars = [vars_q, vars_dq, vars_ddq, params_massa, params_geom, params_inercia];

matlabFunction(Tau, 'File', 'TauFun', 'Vars', {input_vars});







        