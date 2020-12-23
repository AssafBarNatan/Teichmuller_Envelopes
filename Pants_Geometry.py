import numpy as np
import sympy as sp
import time
import itertools
import matplotlib.pyplot as plt
import plotly.express as px
import pandas as pd
import re
t = sp.symbols('t')

class Curve:
    def __init__(self, name, length):
        self.name = name
        self.length = length

    def get_length(self):
        return(self.length)

    def get_name(self):
        return(self.name)

    def set_length(self, new_length):
        self.length = new_length

    def set_name(self, new_name):
        self.name = new_name

class Curve_with_twist(Curve):
    def __init__(self, name, length, twist_direction):
        self.name = name
        self.length = length
        self.twist_direction = twist_direction

    def set_twist(self, new_twist_direction):
        self.twist_direction = new_twist_direction


class Pants:
    def __init__(self, curve1, curve2, curve3):
        self.curves = [curve1,curve2,curve3]

    def get_curves(self):
        return(self.curves)

    def set_triangulation(self, basecurve, triangulation, twists):
        if not triangulation in ["3","2"]:
            print("please enter a valid triangulation tag")
            return()

        if not basecurve in set(self.get_curves()):
            print("the triangulation needs to be centered around a valid curve")
            return()

        if not set(twists.keys()) == set(self.get_curves()):
            print("bad twist dictionary: make sure names are consistent with pant cuff names")
            return()

        self.basecurve = basecurve
        self.triangulation = triangulation
        self.twists = twists

    def get_twists(self):
        return(self.twists)

    def get_basecurve(self):
        return(self.basecurve)

    def get_triangulation(self):
        return(self.triangulation)

    def print_pants(self):
        print("curves: " + str(self.get_curves()))
        print("twists: " + str(self.get_twists()))
        print("base: " + str(self.get_basecurve()))
        print("triangulation: " + str(self.get_triangulation()))

#forveery pants:
#forevery lamination:
#pants_with_lamination
#compute_twists() which does::
#forevery curve:
#find the pants containing the curve
#compute_twists()



# shears returns the shear computations, given the completion of the
# pair of pants to a triangulation, and the three cuff lengths. The
# completion tag can be "3symmetric", "2symmetric", or "asymmetric".
# "3symmetric" is the lamination with two infinite leaves converging
# to every cuff. "2symmetric" is the lamination with all leaves
# converging to alpha, and "2symmetric" is the lamination with just
# a single leaf converging to alpha, and three leaves convering to
# beta.


# I use this to loop over all combinations of laminations and twist
# directions
#

def shears(completion, l_alpha, l_beta, l_gamma):
    if completion == "3symmetric":
        return([(-l_alpha - l_beta + l_gamma)/2,
                (-l_alpha + l_beta - l_gamma)/2,
                (l_alpha - l_beta - l_gamma)/2])
        # s_ab, s_ac, s_bc

    if completion == "2symmetric":
        return([(-l_alpha + l_beta + l_gamma)/2,
                -l_beta,
                -l_gamma])
        # s_aa, s_ab, s_ac
    
    if completion == "asymmetric":
        return([(-l_beta + l_alpha + l_gamma)/2,
                -l_alpha,
                -l_gamma])
        # s_bb, s_ab, s_bc

# delta is defined to be the signed distance from q to p. p is 
# defined to be the point of intersection of the normal between 
# beta and alpha. q is defined to be the point of intersection of 
# the geodesic normal to alpha and to the incircle of one of the 
# two triangles forming the lamination.
#
# Note that this distance is uniquely defined up to adding 
# l_alpha/2 and perhaps some shears.

# orientations is an optional input vector that refers to the twisting
# directon around alpha, beta, gamma (respectively). 1 means left, 
# -1 means right.

# if completion is "asymmetric", we assume that all lamination leaves 
# terminate at beta

def get_lengths_from_pants(pants, join_curve):
    if not join_curve in pants.get_curves():
        print("curve not found in pants")
        return()
    l_a = join_curve.get_length()*pants.get_twists()[join_curve] 
    base = pants.get_basecurve()
    if join_curve == base:
        [l_b, l_c] = [c.get_length()*pants.get_twists()[c] for c in pants.get_curves() if c != join_curve]
    else:
        l_b = base.get_length()*pants.get_twists()[base] 
        [l_c] = [c.get_length()*pants.get_twists()[c] for c in pants.get_curves() if (c != join_curve and c != base)]

    return([l_a, l_b, l_c])

def compute_delta(pants, join_curve):
    [l_a, l_b, l_c] = get_lengths_from_pants(pants, join_curve)

    completion = "asymmetric"

    if pants.get_triangulation() == "3":
        completion = "3symmetric"

    if pants.get_triangulation() == "2" and join_curve == pants.get_basecurve():
        completion = "2symmetric"

    if completion == "3symmetric":
        [ab, ac, bc] = shears(completion, l_a, l_b, l_c)
        x = (1 + sp.exp(ab)) / (sp.exp(-l_a) - 1)
        pstar = x + (sp.exp(bc) + sp.exp(-l_b)) / (1 + sp.exp(bc))

    if completion == "2symmetric":
        [aa, ab, ac] = shears(completion, l_a, l_b, l_c)
        x = (1 + sp.exp(ab) + sp.exp(ab+aa) + \
                sp.exp(ab + aa + ac)) /  (sp.exp(-l_a) - 1)
        pstar = x + sp.exp(-l_b)

    if completion == "asymmetric":
        [bb, ab, bc] = shears(completion, l_a, l_b, l_c)
        x = 1 / (sp.exp(-l_a) - 1)
        pstar = x + (sp.exp(bb) + sp.exp(bb + bc) + sp.exp(2*bb + bc) \
                + sp.exp(-l_b)) / (sp.exp(bb) + sp.exp(bb+bc) \
                + sp.exp(2*bb+bc) + 1)

    return(-pants.get_twists()[join_curve]*0.5*sp.log((x+1)*pstar))

# Returns the twisting parameter, given full lamination data on each
# side

def germ(pants1, pants2, join_curve):
    twist0 = 0
    twist = compute_delta(pants1, join_curve)+compute_delta(pants2, join_curve)
    twist_t = twist0*sp.exp(t) - sp.exp(t)*twist.subs({t:0}) + twist
    return(float(sp.N(sp.diff(twist_t,t).subs({t:0}))))


