#!/usr/bin/env python3

from sympy import *

class Basis:
    def __init__(self, name):
        self.name   = 'd' + name
        self.sym    = symbols(name)
        self.seq    = None
    
    def __mul__(self, coeff):
        return Diff_Form({(self,): coeff}, 1, self.seq)
    def __rmul__(self, coeff):
        return self.__mul__(coeff)

    def __repr__(self):
        return self.name

    def __str__(self):
        return self.name

class Diff_Form:
    def __init__(self, value, degree, seq):
        self.value  = value
        self.degree = degree
        self.seq    = seq

    def __mul__(self, other):
        if self.degree:
            if isinstance(other, Diff_Form):
                if other.degree:
                    return wedge(self, other)
                else:
                    result = Diff_Form(self.value.copy(), self.degree, self.seq)
                    for basis in result.value:
                        result.value[basis] *= other.value
                    return result.simplify()
            else:
                result = Diff_Form(self.value.copy(), self.degree, self.seq)
                for basis in result.value:
                    result.value[basis] *= other
                return result.simplify()
        else:
            if isinstance(other, Diff_Form):
                if other.degree:
                    result = Diff_Form(other.value.copy(), other.degree, other.seq)
                    for basis in result.value:
                        result.value[basis] *= self.value
                    return result.simplify()
                else:
                    return Diff_Form(self.value * other.value, 0, self.seq)
            else:
                return Diff_Form(self.value * other, 0, self.seq)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        return self * (1 / other)

    def __rtruediv__(self, other):
        assert self.degree == 0
        return other / self.value

    def __add__(self, other):
        assert self.degree == other.degree
        if self.degree == 0:
            return Diff_Form(self.value + other.value, 0, self.seq)
        result = Diff_Form(self.value.copy(), self.degree, self.seq)
        for basis in other.value:
            if basis in result.value:
                result.value[basis] +=  other.value[basis]
            else:
                result.value[basis] =   other.value[basis]
        return result

    def __sub__(self, other):
        return self + (- other)

    def __neg__(self):
        if self.degree == 0:
            return Diff_Form(- self.value, 0, self.seq)
        result = Diff_Form(dict(), self.degree, self.seq)
        for basis in self.value:
            result.value[basis] = - self.value[basis]
        return result

    def d(self):
        result = Diff_Form(dict(), self.degree + 1, self.seq)
        if self.degree == 0:
            for basis in self.seq:
                result.value[(basis,)] = sympify(self.value).diff(basis.sym)
            return result.simplify()
        for basis in self.value:
            for nextb in self.seq:
                if nextb in basis:
                    continue
                result.value[(nextb,) + basis] = sympify(self.value[basis]).diff(nextb.sym)
        return result.simplify()

    def simplify(self):
        if self.degree == 0:
            return self
        result = dict()
        for basis in self.value:
            pre_seq = tuple(map(self.seq.index, basis))
            new_basis = tuple(self.seq[index] for index in sorted(pre_seq))
            inv_count = 0
            for i in range(len(basis)):
                inv_count += sum(1 for less in pre_seq[i + 1 :] if less < pre_seq[i])
            new_value = (-1 if inv_count % 2 else 1) * self.value[basis]
            if new_basis in result:
                result[new_basis] += new_value
            else:
                result[new_basis] =  new_value
        for basis in list(result.keys()):
            if not result[basis]:
                del result[basis]
        self.value = result
        return self

    def __repr__(self):
        if self.degree == 0:
            return str(self.value)
        if not len(self.value):
            return '0'
        return ' + '.join(f"({self.value[basis]}) {'∧'.join(map(str, basis))}" for basis in self.value)

    def __str__(self):
        return self.__repr__()

def wedge(form_1: Diff_Form, form_2: Diff_Form):
    assert form_1.seq == form_2.seq
    result = Diff_Form(dict(), form_1.degree + form_2.degree, form_1.seq)
    for basis_1 in form_1.value:
        for basis_2 in form_2.value:
            basis = basis_1 + basis_2
            if len(basis) != len(set(basis)):
                continue
            if basis in result.value:
                result.value[basis] +=  form_1.value[basis_1] * form_2.value[basis_2]
            else:
                result.value[basis] =   form_1.value[basis_1] * form_2.value[basis_2]
    return result.simplify()

def set_seq(*basis):
    for each in basis:
        each.seq = basis

"""
x, y, z = symbols('x y z')
dx = Basis('x')
dy = Basis('y')
dz = Basis('z')

set_seq(dx, dy, dz)
"""
