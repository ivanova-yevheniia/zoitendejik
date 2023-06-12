from methods import *
from scipy.optimize import linprog
import numpy as np
import matplotlib.pyplot as plt


def s(x, gg1, gg2, gf, active_g):
    '''
    :param active: активні обмеження
    :param x: точка
    :param gg: градієнт активного обмеження
    :param gf: градієнт функції
    :return: напрямок
    '''
    A = []
    b = []
    A.append([gf[0], -gf[0], gf[1], -gf[1], 1])
    b.append(0)
    if active_g[0] == True:
        A.append([-gg1[0], gg1[0], -gg1[1], gg1[1], 1])
        b.append(0)
    else: pass
    if active_g[1] == True:
        A.append([-gg2[0], gg2[0], -gg2[1], gg2[1], 1])
        b.append(0)
    else: pass
    c = [0, 0, 0, 0, -1] # які змінні ми максимізуємо
    bounds = [(0, 1), (0, 1), (0, 1), (0, 1), (0, None)] # межі змінних
    result = linprog(c, A_ub=A, b_ub=b, bounds=bounds, method='highs')
    if result.success:
        return Dim2(result.x[0]-result.x[1], result.x[2]-result.x[3])




h = [0.1, 0.01, 0.001]
alpha = [-0.1, 0, 0.1, 0.5, 1]
derivate = [forward_difference, central_difference, backward_difference]
error = [0.1, 0.01, 0.001, 0.0001, 0.00001]
lam_func = [dsk_paule, gold_sect_algo]
x_0 = Dim2(1, 1)

def zoitendejik(x, h, alpha, error):
    points = []
    points.append(x)
    function_values = []
    function_values.append(f(x.X_1, x.X_2))
    print('........ПОЧАТКОВІ УМОВИ........')
    print('x: ', x)
    print('f: ', f(x.X_1, x.X_2))
    print('     ITERATION 1     ')
    iteration = 0
    while in_(x):
        #get gradients
        gf = forward_difference(x, h, f)
        print('gf', gf)
        gg1 = forward_difference(x, h, g1)
        print('gg1', gg1)
        gg2 = forward_difference(x, h, g2)
        print('gg2', gg2)
        #get s
        active_g = restrictions(x, g) #індекси активних обмежень
        print('Діючі обмеження: ', active_g)
        new_s = s(x, gg1, gg2, gf, active_g)
        print('s'+str(iteration+1), new_s)
        del_lam = delta_lam(x, new_s, alpha)
        print('delta lambda: ', del_lam)
        undef_interval = sven_algo(x, del_lam, new_s)
        print('undefined interval: ', undef_interval)
        #interval = check_interval(x, new_s, undef_interval, 0.1)
        #print('Adjusted Interval: ', interval)
        lam = gold_sect_algo(undef_interval, x, new_s, error)
        l = (lam[0]+lam[1])/2
        print('lambda: ', l)
        x = new_point(x, new_s, l)
        points.append(x)
        print('x: ', x)
        function_values.append(f(x.X_1, x.X_2))
        print('f: ', f(x.X_1, x.X_2))
        iteration += 1
        print('-----------------------')
        print('     ITERATION '+str(iteration+1)+'     ')
    print('minimization completed!')
    print('x: ')
    for point in points:
        print(point)
    print('f: ', function_values)
    print('Загальна кількість ітерацій: ', iteration)
    return points, function_values

p, f = zoitendejik(Dim2(1, 1), 0.01, 0.001, 0.01)

x2 = np.linspace(-5, 5, 100)
x1 = (x2 ** 2 + 1)/2
x1_ = np.linspace(-5, 5, 100)
x2_ = (-0.8 * x1_ ** 2 + 9)/2
plt.plot(x1_, x2_)
plt.plot(x1, x2)
x_ = [el.X_1 for el in p]
y_ = [el.X_2 for el in p]
plt.plot(x_, y_)
plt.plot(2.5, 2, 'ro')
plt.show()



