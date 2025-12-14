import sympy as sp
import math

n = 4

theta = [0,0,sp.symbols(f'theta3'),sp.symbols(f'theta4')]
d = [sp.symbols(f'l1'),sp.symbols(f'l2'),sp.symbols(f'l3'),sp.symbols(f'l5')]
alpha = [-math.pi/2,math.pi/2,-math.pi/2, 0]
a = [0, 0, sp.symbols('l4'), 0]

def DH_matrix(theta, d, alpha, a):
    return sp.Matrix([
        [sp.cos(theta), -sp.sin(theta)*sp.cos(alpha),  sp.sin(theta)*sp.sin(alpha), a*sp.cos(theta)],
        [sp.sin(theta),  sp.cos(theta)*sp.cos(alpha), -sp.cos(theta)*sp.sin(alpha), a*sp.sin(theta)],
        [0,              sp.sin(alpha),               sp.cos(alpha),               d],
        [0,              0,                           0,                           1]
    ])


A01 = sp.nsimplify(DH_matrix(theta[0], d[0], alpha[0], a[0]), tolerance=1e-10, rational=False)
A12 = sp.nsimplify(DH_matrix(theta[1], d[1], alpha[1], a[1]), tolerance=1e-10, rational=False)
A23 = sp.nsimplify(DH_matrix(theta[2], d[2], alpha[2], a[2]), tolerance=1e-10, rational=False)
A34 = sp.nsimplify(DH_matrix(theta[3], d[3], alpha[3], a[3]), tolerance=1e-10, rational=False)


A02 = A01 * A12
A03 = A02 * A23
A04 = A03 * A34


sp.pprint(A01)
print("###########################")
sp.pprint(A12)
print("###########################")
sp.pprint(A23)
print("###########################")
sp.pprint(A34)
print("debaixo A02")
sp.pprint(A02)
print("debaixo A03")
sp.pprint(A03)
print("debaixo A04")
sp.pprint(A04)