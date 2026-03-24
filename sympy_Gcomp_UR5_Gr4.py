from sympy import symbols, Matrix, cos, sin, pi, simplify, diff

# Define symbolic joint variables
q1, q2, q3, q4, q5, q6 = symbols('q1 q2 q3 q4 q5 q6')
theta = [q1, q2, q3, q4, q5, q6]

# DH Parameters
alpha_DH = [pi/2, 0, 0, pi/2, -pi/2, 0]
d_DH     = [0.089159, 0, 0, 0.10915, 0.09465, 0.0823]
a_DH     = [0, -0.425, -0.39225, 0, 0, 0]

# Masses and gravity vector
m = [3.7, 8.393, 2.33, 1.1219, 1.1219, 0.1879]
g_v = Matrix([0, 0, 9.81])

# Center of Mass (CoM) vectors
center_of_mass = [
    Matrix([0.0, -0.02561, 0.00193, 1]),
    Matrix([0.2125, 0.0, 0.11336, 1]),
    Matrix([0.15, 0.0, 0.0265, 1]),
    Matrix([0.0, -0.0018, 0.01634, 1]),
    Matrix([0.0, 0.0018, 0.01634, 1]),
    Matrix([0.0, 0.0, -0.001159, 1])
]

# DH transformation function
def dh_matrix(a, alpha, d, theta):
    return Matrix([
        [cos(theta), -sin(theta)*cos(alpha), sin(theta)*sin(alpha), a*cos(theta)],
        [sin(theta), cos(theta)*cos(alpha), -cos(theta)*sin(alpha), a*sin(theta)],
        [0, sin(alpha), cos(alpha), d],
        [0, 0, 0, 1]
    ])

# Forward kinematics and potential energy calculation
H = Matrix.eye(4)
P = []

for i in range(6):
    A_i = dh_matrix(a_DH[i], alpha_DH[i], d_DH[i], theta[i])
    H = H * A_i
    rG_i = H * center_of_mass[i]
    rG_world = rG_i[:3]
    P.append(m[i] * g_v.dot(rG_world))

P_total = sum(P)

# Gravity compensation torques: ∂P_total/∂qi
G = [simplify(diff(P_total, theta[i])) for i in range(6)]

# Output torques
for i, g in enumerate(G, start=1):
    print(f"G[{i}] =", g)
