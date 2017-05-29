import numpy as np
import matplotlib.pyplot as plt
import func
import mag
n = 8
x_zn = np.zeros(n)
y_zn = np.zeros(n)
a = 0.1
b = 0.15
h = (b - a) / (n - 1)
t = a
for i in range(0, n, 1):
    x_zn[i] = t
    y_zn[i] = func.f(x_zn[i])
    t += h
gr_x = np.zeros((n-1)*5+1)
graph = np.zeros((4, (n-1)*5+1))


def pollag(x, xv, yv, n):
    lag_pol = 0.
    bas_pol = 1.
    for i in range(0, n, 1):
        for j in range(0, n, 1):
            if i != j:
                bas_pol *= (x - xv[j]) / (xv[i] - xv[j])
        lag_pol += bas_pol * yv[i]
        bas_pol = 1.
    return lag_pol


def lagr():
    print 'lagrang\n'
   # print x_zn
   # print y_zn
    t = a
    i = 0
    while abs(t - b) > func.e:
        l = pollag(t, x_zn, y_zn, n)
        eps = abs(func.f(t) - l)
        graph[0, i] = eps
        gr_x[i] = t
        i += 1
        print 'iteration: ' + str(i) + ' x= ' + str(t) + ' f(x)= ' + str(func.f(t)) + ' Lag= ' + str(l) + ' e= ' + str(eps) + '\n'
        t += h / 5

    gr_x[i] = t
    i += 1
    l = pollag(t, x_zn, y_zn, n)
    eps = abs(func.f(t) - l)

    print 'iteration: ' + str(i) + ' x= ' + str(t) + ' f(x)= ' + str(func.f(t)) + ' Lag= ' + str(l) + ' e= ' + str(eps) + '\n'


def spl(x, xv, yv, n):
    alpha = np.zeros(n - 1)
    beta = np.zeros(n - 1)
    c1 = np.zeros(n)
    d1 = np.zeros(n)
    b1 = np.zeros(n)
    for i in range(1, n-1, 1):
        h_i = xv[i] - xv[i - 1]
        h_i1 = xv[i + 1] - xv[i]
        A = h_i
        C = 2.0 * (h_i + h_i1)
        B = h_i1
        F = 6.0 * (((yv[i + 1] - yv[i]) / h_i1) - ((yv[i] - yv[i - 1]) / h_i))
        z = (A * alpha[i - 1] + C)
        alpha[i] = -B / z
        beta[i] = (F - A * beta[i - 1]) / z
    c1[n - 1] = (F - A * beta[n - 2]) / (C + A * alpha[n - 2])
    for i in range(n - 2, 0, -1):
        c1[i] = alpha[i] * c1[i + 1] + beta[i]
    for i in range(n - 1, 0, -1):
        h_i = xv[i] - xv[i - 1]
        d1[i] = (c1[i] - c1[i - 1]) / h_i
        b1[i] = h_i * (2.0 * c1[i] + c1[i - 1]) / 6.0 + (yv[i] - yv[i - 1]) / h_i

    for i in range(1, n, 1):
        if xv[i - 1] <= x <= xv[i] + func.e:
            return yv[i] + b1[i]*(x-xv[i]) + (c1[i]*(x-xv[i]) ** 2)/2. + (d1[i]*(x-xv[i]) ** 3)/6.


def splain():
    print'splain\n'
    t = a
    i = 0
    while abs(t - b) > func.e:
        l = spl(t, x_zn, y_zn, n)
        eps = abs(func.f(t) - l)
        graph[1, i] = eps
        i += 1
        print 'iteration: ' + str(i) + ' x= ' + str(t) + ' f(x)= ' + str(func.f(t)) + ' splain= ' + str(l) + ' e= ' + str(eps) + '\n'
        t += h / 5
    i += 1
    l = spl(t, x_zn, y_zn, n)
    eps = abs(func.f(t) - l)

    print 'iteration: ' + str(i) + ' x= ' + str(t) + ' f(x)= ' + str(func.f(t)) + ' splain= ' + str(l) + ' e= ' + str(eps) + '\n'


def newu(x, xv, yv, n):
    res = yv[0]
    for i in range(1, n, 1):
        F = 0
        for j in range(0, i+1, 1):
            den = 1
            for k in range(0, i+1, 1):
                if k != j:
                    den *= (xv[j] - xv[k])
            F += yv[j] / den
        for k in range(0, i, 1):
            F *= (x - xv[k])
        res += F
    return res


def newtonup():
    print'newton up\n'
    t = a
    i = 0
    while abs(t - b) > func.e:
        l = newu(t, x_zn, y_zn, n)
        eps = abs(func.f(t) - l)
        graph[2, i] = eps
        i += 1
        print 'iteration: ' + str(i) + ' x= ' + str(t) + ' f(x)= ' + str(func.f(t)) + ' newup= ' + str(l) + ' e= ' + str(eps) + '\n'
        t += h / 5
    i += 1
    l = newu(t, x_zn, y_zn, n)
    eps = abs(func.f(t) - l)
    print 'iteration: ' + str(i) + ' x= ' + str(t) + ' f(x)= ' + str(func.f(t)) + ' newup= ' + str(l) + ' e= ' + str(eps) + '\n'


def newtondown():
    print'newton down\n'
    t = b
    x_d = np.zeros(n)
    y_d = np.zeros(n)
    for i in range(0, n, 1):
        x_d[i] = x_zn[n - i - 1]
        y_d[i] = y_zn[n - i - 1]
    i = 0
    while abs(a - t) > func.e:
        l = newu(t, x_d, y_d, n)
        eps = abs(func.f(t) - l)
        graph[3, (n-1)*5 - i] = eps
        i += 1
        print 'iteration: ' + str(i) + ' x= ' + str(t) + ' f(x)= ' + str(func.f(t)) + ' newdown= ' + str(l) + ' e= ' + str(eps) + '\n'
        t -= h / 5
    i += 1
    l = newu(t, x_zn, y_zn, n)
    eps = abs(func.f(t) - l)

    print 'iteration: ' + str(i) + ' x= ' + str(t) + ' f(x)= ' + str(func.f(t)) + ' newdown= ' + str(l) + ' e= ' + str(eps) + '\n'


def sup():
    x = a
    maxi = mag.df(x)
    while abs(x - b) > func.e:
        if mag.df(x) > maxi:
            maxi = mag.df(x)
        x += h / 5
    return maxi


def mg(x, x_zn, n, maxi):
    w = 1
    for i in range(0, n, 1):
        w *= (x - x_zn[i])
    #print mag.df(x)
    return maxi * abs(w) / 720


def magoranta():
    print 'Magoranta\n'
    maxii = sup()
    ma = np.zeros((n-1)*5+1)
    t = a
    i = 0
    while abs(t - b) > func.e:
        i += 1
        t += h / 5
        m = mg(t, x_zn, n, maxii)
        ma[i] = m
        print 'Magoranta ' + str(m)
    plt.plot(gr_x, ma)
    plt.show()


def gr():
    print 'orange splain, red newtondown lagrang newtonup'
    for i in range(0, 4, 1):
        plt.plot(gr_x, graph[i])
    plt.show()
