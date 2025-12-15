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


user_input = input('Vetor de posições de junta: '); % digite no formato [ l1, l2, theta3, theta4]
                                                    %ex: [0.47, 0.5, -pi, pi]
l1 = user_input(1);
l2 = user_input(2);
theta3 = user_input(3);
theta4 = user_input(4);

c3 = cos(theta3);
s3 = sin(theta3);
c4 = cos(theta4);
s4 = sin(theta4);

T = [
        c3*c4,      -s4*c3,     -s3,    l4*c3 - l5*s3;
        
        s3*c4,      -s3*s4,      c3,    l2 + l4*s3 + l5*c3;
        
       -s4,         -c4,         0,     l1+l3;
        
        0,           0,          0,     1
    ];

x = T(1, 4)
y = T(2, 4)
z = T(3, 4)

%robot.fkine(user_input)
robot.plot(user_input);
%robot.teach(user_input)