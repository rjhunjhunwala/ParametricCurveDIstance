#I'm  going to do this in python because, idk, but I'm going to avoid using scientific computing libraries, because I guess that's cheating.

import math # If only this also gave me my MS in math
import random

class Curve:
    # a curve is defined, by a parametric equations and it's first and second derivatives as functions of time
    def __init__(self, f_x = None, f_y = None, fp_x = None,fp_y = None, fpp_x = None, fpp_y = None, min_t = 0, max_t = 1):
        self.f_x = f_x
        self.f_y = f_y
        self.fp_x = fp_x
        self.fp_y = fp_y
        self.fpp_x = fpp_x
        self.fpp_y = fpp_y
        self.min_t = min_t
        self.max_t = max_t

def get_circle(center_x, center_y, radius):
    f_x = lambda t: center_x + radius * math.cos(t)
    f_y = lambda t: center_y + radius * math.sin(t)
    fp_x = lambda t: radius * -1 * math.sin(t)
    fp_y = lambda t: radius * math.cos(t)
    fpp_x = lambda t: radius * -1 * math.cos(t)
    fpp_y = lambda t: radius * -1 * math.sin(t)
    return Curve(f_x = f_x, f_y = f_y, fp_x = fp_x, fp_y = fp_y, fpp_x = fpp_x, fpp_y = fpp_y, min_t = 0, max_t = 2 * math.pi)

def scale_vec(scalar, vec):
    for i in range(len(vec)):
        vec[i] *= scalar
    return vec

def scale_mat(scalar, mat):
    for row in mat:
        scale_vec(scalar, row)
    return mat

def invert2(mat):
    a, b, c, d = mat[0][0], mat[0][1], mat[1][0], mat[1][1]
    det = a * d - b * c
    return scale_mat(1/det, [[d, -b],[-c, a]])

def dot(a, b):
   return sum(x * y for x, y in zip(a, b))


def add_vec(a, b):
    return [x + y for x, y in zip(a, b)]

def multiply_mat_vec2(mat, vec):
    return [dot(row, vec) for row in mat]

def get_distance(curve1, curve2):
    def function(vec):
        a, b = vec
        return (curve1.f_x(a) - curve2.f_x(b)) ** 2 + (curve1.f_y(a) - curve2.f_y(b)) ** 2
    def gradient(vec):
        a,b = vec
        f, g = curve1, curve2
        first = 2 * f.fp_x(a) * (f.f_x(a) - g.f_x(b)) + 2 * f.fp_y(a) * (f.f_y(a) - g.f_y(b))
        second = 2 * g.fp_x(b) * (g.f_x(b) - f.f_x(a)) + 2 * g.fp_y(b) * (g.f_y(b) - f.f_y(a))
        return [first, second]
    def hessian(vec):
        a,b = vec
        f, g = curve1, curve2
        H_a = 2 * f.fpp_x(a) * (f.f_x(a) - g.f_x(b)) + 2 * (f.fp_x(a))**2 + 2 * f.fpp_y(a) * (f.f_y(a) - g.f_y(b)) + 2 * (f.fp_y(a))**2
        H_b = - 2 * g.fp_x(b) * f.fp_x(a) - 2 * g.fp_y(b) * f.fp_y(a)
        H_c = H_b # Hessian is symmetric
        H_d = 2 * g.fpp_x(b) * (g.f_x(b) - f.f_x(a)) + 2 * (g.fp_x(b))**2 + 2 * g.fpp_y(b) * (g.f_y(b) - f.f_y(a)) + 2 * (g.fp_y(b))**2

        return [[H_a, H_b], [H_c, H_d]]

    guesses = [[random.uniform(curve1.min_t, curve1.max_t), random.uniform(curve1.min_t, curve2.max_t)] for i in range(10)]
    values = [optimize2dNewton(function, gradient, hessian, guess = guess, iters = 10, newton = True) for guess in guesses]
    parameters, distance = min(values, key = lambda tup: tup[1])
    return parameters, distance ** (.5)

def optimize2dNewton(function, gradient, hessian, guess = [0, 0], iters = 10, newton = True):
    alpha = 1 if newton else .03
    for i in range(iters):
        g = gradient(guess)

        if newton:
            H = hessian(guess)
            Hinv = invert2(H)
            step = multiply_mat_vec2(Hinv, g)
            guess = add_vec(guess, scale_vec(-alpha, step))
        else:
            guess = add_vec(guess, scale_vec(-alpha, g))
    return guess, function(guess)


origin = get_circle(5,12, 1)
other = get_circle(0,0, 1)
print(get_distance(origin, other)) # expect [a, b] 8