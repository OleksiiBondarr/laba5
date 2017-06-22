import sympy


def df(c):
    x = sympy.Symbol('x')
    r = sympy.diff(f(x), x, 8)
    return r.evalf(subs={x: c})


def df1():
    x = sympy.Symbol('x')
    r = sympy.diff(f(x), x, 8)
    return r

#change here
def f(x):
    return (1 / sympy.tan(2))**(1/3) - ((sympy.cos(10*x))**2)/(20*sympy.sin(20*x))
