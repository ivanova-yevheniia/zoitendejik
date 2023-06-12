from task import Dim2
from collections import OrderedDict

'''Calculate Rosenbock function'''


# def f(x1, x2):
#    x = Dim2(x1, x2)
#   return 100 * (x1 ** 2 + x2) ** 2 + (x1- 1) ** 2

def f(x1, x2):
    return (x1 - 3) ** 2 + (x2 - 3) ** 2


#def f(x_1, x_2):
 #   return round(3 * (x_1 - 8) ** 2 + x_1 * x_2 + 7 * x_2 ** 2, 6)


'''Calculate Gradinets'''


# default gradient
def gradient(x):
    df_dx1 = 400 * (x.X_1 ** 2 - x.X_2) * x.X_1 + 2 * (x.X_1 - 1)
    df_dx2 = -200 * (x.X_1 ** 2 - x.X_2)
    return [df_dx1, df_dx2]


# right schema
def forward_difference(x, h, f):
    df_dx1 = (f(x.X_1 + h, x.X_2) - f(x.X_1, x.X_2)) / h
    df_dx2 = (f(x.X_1, x.X_2 + h) - f(x.X_1, x.X_2)) / h
    return [df_dx1, df_dx2]


# left schema
def backward_difference(x, h, f):
    df_dx1 = (f(x.X_1, x.X_2) - f(x.X_1 - h, x.X_2)) / h
    df_dx2 = (f(x.X_1, x.X_2) - f(x.X_1, x.X_2 - h)) / h
    return [df_dx1, df_dx2]


# central schema
def central_difference(x, h, f):
    df_dx1 = (f(x.X_1 + h, x.X_2) - f(x.X_1 - h, x.X_2)) / 2 * h
    df_dx2 = (f(x.X_1, x.X_2 + h) - f(x.X_1, x.X_2 - h)) / 2 * h
    return [df_dx1, df_dx2]


def delta_lam(x, grad, alpha):
    x_norm = (x.X_1 ** 2 + x.X_2 ** 2) ** 0.5
    grad_norm = (grad.X_1 ** 2 + grad.X_2 ** 2) ** 0.5
    return round(alpha * x_norm / grad_norm, 5)


def new_point(x_old, grad, lam):
    return Dim2(x_old.X_1 + lam * grad.X_1, x_old.X_2 + lam * grad.X_2)

def g1(x1, x2):
    return 2 * x1 - x2 ** 2 - 1

def g2(x1, x2):
    return -0.8 * x1 ** 2 - 2 * x2 + 9

g = [g1, g2]


def in_(x):
    if g1(x.X_1, x.X_2) >= 0 and g2(x.X_1, x.X_2) >= 0:
        #print('Ми в допустимій області')
        return True
    print('Ми поза допустимої області')
    return False


def restrictions(x, g):
    aktiv_g = []
    for element in g:
        if element(x.X_1, x.X_2) == 0:
            aktiv_g.append(True)
        else: aktiv_g.append(False)
    return aktiv_g

def check_zone(x):
    if g1(x.X_1, x.X_2)>=0 and g2(x.X_1, x.X_2)>=0:
        return True
    else: return False

def check_interval(x, grad, interval, step):
    a = interval[0]
    b = interval[1]
    a_x = new_point(x, grad, a)
    b_x = new_point(x, grad, b)
    while check_zone(a_x) == False or check_zone(b_x) == False:
        if g1(a_x.X_1, a_x.X_2)<0 or g2(a_x.X_1, a_x.X_2)<0:
            a = a+abs(step)
            #print(a)
        else: a = interval[0]

        if g1(b_x.X_1, b_x.X_2)<0 or g2(b_x.X_1, b_x.X_2)<0:
            b = b-abs(step)
            #print(b)
        else: b = interval[1]

    return [a, b]

def sven_algo(x, del_lam, grad):
    f_values = []
    lams = []
    fvl = []

    f_values.append(f(x.X_1, x.X_2))
    lams.append(0)

    new_x = new_point(x, grad, del_lam)  # lambda_0 + delta_lambda
    f_values.append(f(new_x.X_1, new_x.X_2))

    if f_values[1] > f_values[0]:
        del_lam = -del_lam
        new_x = new_point(x, grad, del_lam)  # lambda_0 - delta_lambda
        f_values.pop()
        f_values.append(f(new_x.X_1, new_x.X_2))

    lams.append(del_lam)

    while f_values[-1] < f_values[-2]:
        del_lam = 2*del_lam
        new_lam = lams[-1] + del_lam
        x_new = new_point(x, grad, new_lam)
        f_values.append(f(x_new.X_1, x_new.X_2))
        if in_(x_new):
            lams.append(new_lam)
        else:
            f_values.pop()
            break

    if f_values[-1] > f_values[-2]:
        last_lam = (lams[-1] + lams[-2]) / 2
        lams.append(last_lam)

    sorted_lams = sorted(lams)
    min_index = f_values.index(min(f_values))
    min_lam = lams[min_index]

    if sorted_lams.index(min_lam) == len(sorted_lams)-1:
        return [sorted_lams[sorted_lams.index(min_lam) - 1], sorted_lams[sorted_lams.index(min_lam)]]
    elif sorted_lams.index(min_lam) == 0:
        return [sorted_lams[sorted_lams.index(min_lam)], sorted_lams[sorted_lams.index(min_lam) + 1]]
    else:
        a = sorted_lams[sorted_lams.index(min_lam) - 1]
        b = sorted_lams[sorted_lams.index(min_lam) + 1]
        return [a, b]


def gold_sect_algo(interval, x, grad, error):
    b = interval[1]
    a = interval[0]
    while abs(b - a) > error:
        lambda_1 = round(a + 0.382 * (b - a), 10)
        lambda_2 = round(a + 0.618 * (b - a), 10)
        if f(x.X_1 + grad.X_1 * lambda_1, x.X_2 + grad.X_2 * lambda_1) < f(x.X_1 + grad.X_1 * lambda_2, x.X_2 + grad.X_2 * lambda_2):
            b = lambda_2
        elif f(x.X_1 + grad.X_1 * lambda_1, x.X_2 + grad.X_2 * lambda_1) == f(x.X_1 + grad.X_1 * lambda_2, x.X_2 + grad.X_2 * lambda_2):
            a = lambda_1
            b = lambda_2
        else:
            a = lambda_1
    return [a, b]


def dsk_paule(interval, x, grad, error):
    x1 = interval[0]
    x2 = (abs(interval[1]) + abs(interval[0])) / 2
    x3 = interval[1]
    delta_x = x2-x1

    f1 = f(x.X_1+x1*grad.X_1, x.X_2+x1*grad.X_2)
    f2 = f(x.X_1 + x2 * grad.X_1, x.X_2 + x2 * grad.X_2)
    f3 = f(x.X_1 + x3 * grad.X_1, x.X_2 + x3 * grad.X_2)

    x_star = x2 + delta_x/2 *(f1-f3)/(f1-2*f2+f3)
    fx_star = f(x.X_1 + x_star * grad.X_1, x.X_2 + x_star * grad.X_2)
    if abs(f2-fx_star)>error or abs(x2-x_star)>error:
        f_val = [f1, f2, f3, fx_star]
        x_val = [x1, x2, x3, x_star]
        x2_index = f_val.index(min(f_val))
        x2 = x_val[x2_index]
        f_x2 = f(x.X_1+x2*grad.X_1, x.X_2+x2*grad.X_2)
        x_val = sorted(x_val)
        if x_val.index(x2) == 0 or x_val.index(x2) == len(x_val)-1:
            return x2
        else:
            x1 = x_val[x_val.index(x2)-1]
            f_x1 = f(x.X_1 + x1 * grad.X_1, x.X_2 + x1 * grad.X_2)
            x3 = x_val[x_val.index(x2) + 1]
            f_x3 = f(x.X_1 + x3 * grad.X_1, x.X_2 + x3 * grad.X_2)
            a1 = (f_x2-f_x1)/(x2-x1)
            a2 = 1/(x3-x2) * ((f_x3-f_x1)/(x3-x1) - (f_x2-f_x1)/(x2-x1))
            x_star = (x1+x2)/2 - a1/(2*a2)
    return x_star


derivate_schema = [gradient, forward_difference, backward_difference, central_difference]

