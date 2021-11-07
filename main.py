import scipy as sp
import numpy as np

a, b = 0, 1
N = 10
h = (b - a)/N
lmbd = 1


def K(x,s):
    return (1 + x + s) ** (-1)


def f(x):
    return 1 + x


# МЕТОД МЕХАНИЧЕСКИХ КВАДРАТУР ДЛЯ ИУФ-2 #

print("\nМетод механических квадратур для ИУФ-2\n")
A = []
for i in range(N + 1):
    A.append(0) if i == 0 else A.append(h)
A = np.array(A)
print("Вектор коэффициентов Ai:\n", A)

x = []
for i in range(N + 1):
    x.append(a + i * h)
x = np.array(x)
print("\nВектор узлов xi:\n", x)

B = []
for i in range(N + 1):
    Bi = []
    for j in range(N + 1):
        if i == j:
            Bi.append(1 - lmbd * A[i] * K(x[i],x[i]))
        else:
            Bi.append(- lmbd * A[i] * K(x[i], x[j]))
    B.append(Bi)
B = np.array(B)
print("\nМатрица коэффициентов для yi:\n", B)

f_i = []
for i in range(N + 1):
    f_i.append(f(x[i]))
f_i = np.array(f_i)
print("\nВектор значений функции f в узлах:\n", f_i)

y = np.linalg.inv(B).dot(f_i)
print("\nВектор yi:\n", y)

point = (a + b)/2.2
y_inpoint = 0
for k in range(N + 1):
    y_inpoint += A[k] * K(point,x[k]) * y[k]
y_inpoint += f(point)
print("\nРешение в точке x*:\n", y_inpoint)


def left_rectangles(ph, n, xk):
    Q = 0
    for i in range(n+1):
        if i>0:
            Q += ph/n * K(xk, x[i])
    return Q


# МЕТОД ПОСЛЕДОВАТЕЛЬНЫХ ПРИБЛИЖЕНИЙ ДЛЯ ИУФ-2 #

# максимум функции достигается в x=0, y=0. значит, М=1
print("\n\nМетод последовательных приближений для ИУФ-2")

M = 1
q = (b - a) * lmbd
# print("q =", q)
n = 5

Phi = []
for i in range(n+1):
    Phi_i = []
    if i == 0:
        for j in range(N + 1):
            Phi_i.append(f_i[j])
    else:
        for j in range(N + 1):
            Q = 0
            for k in range(1, N + 1):
                Q += Phi[i - 1][k] * K(x[j], x[k])
            Phi_i.append(h * Q)
    Phi_i = np.array(Phi_i)
    Phi.append(Phi_i)
Phi = np.array(Phi)
print("\nМатрица значений функций phi_i:\n", Phi)

phi = []
phi.append(f(point))
for i in range(1, n + 1):
    Q = 0
    for j in range(1, N + 1):
        Q += Phi[i - 1][j] * K(point, x[j])
    phi.append(h * Q)
phi = np.array(phi)
print("\nЗначения phi:\n", phi)

y2_inpoint = sum(phi)
print("\nРешение в точке x*:\n", y2_inpoint)

#max_K = 1 # x=0, y=0
#max_f = 2 # x=1
max_K = K(0.01, 0.01)
max_f = f(0.99)
q = lmbd * max_K * (b-a)
print("q = ", max_K, " < 1")
eps = max_f * q ** (n + 1)/(1 - q)
print("eps <= ", eps)

# МЕТОД МЕХАНИЧЕСКИХ КВАДРАТУР ДЛЯ ИУВ-2 #

print("\n\nМетод механических квадратур для ИУВ-2")
BVolt = []
for i in range(N + 1):
    Bj = []
    for j in range(i+1):
        if i == j:
            Bj.append(1 - lmbd * A.item(i) * K(x.item(i), x.item(i)))
        else:
            Bj.append(- lmbd * A.item(i) * K(x.item(i), x.item(j)))
    if len(Bj) < N + 1:
        for j in range(len(Bj), N + 1):
            Bj.append(0)
    BVolt.append(np.array(Bj))
BVolt = np.array(BVolt)
print("\nМатрица системы:\n", BVolt)

yVolt = np.linalg.inv(BVolt).dot(f_i)
print("\nВектор значений yi:\n", yVolt)

y3_inpoint = 0
for k in range(N + 1):
    y3_inpoint += A[k] * K(point, x[k]) * y[k]
y3_inpoint += f(point)
print("\nРешение в точке x*:\n", y3_inpoint)


# МЕТОД ПОСЛЕДОВАТЕЛЬНЫХ ПРИБЛИЖЕНИЙ ДЛЯ ИУВ-2 #

print("\n\nМетод последовательных приближений для ИУВ-2")
phi = []
phi.append(f(point))

Phi = []
for i in range(n + 1):
    Phi_i = []
    if i == 0:
        for j in range(N + 1):
            Phi_i.append(f(x[j]))
    else:
        for j in range(N + 1):
            Q = 0
            for k in range(1, j):
                Q += Phi[i - 1][k] * K(x[j], x[k])
            Phi_i.append(h * Q)
    Phi_i = np.array(Phi_i)
    Phi.append(Phi_i)
Phi = np.array(Phi)
print("\nМатрица значений функций phi_i:\n", Phi)

for i in range(1, n + 1):
    Q = 0
    for j in range(1, N + 1):
        Q += Phi[i - 1][j] * K(point,  x[j])
    phi.append(h * Q)
phi = np.array(phi)
print("\nЗначения phi:\n", phi)

y4_inpoint = sum(phi)
print("\nРешение в точке x*:\n", y4_inpoint)


