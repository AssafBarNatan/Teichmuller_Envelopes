import numpy as np
import sympy as sp
import time
import itertools
import matplotlib.pyplot as plt
import plotly.express as px
import pandas as pd
import json
import re
import os.path


t = sp.symbols('t')
from Pants_Geometry import *
from visualizer import *

from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
def comb(listoftuples,n):
    return(list(itertools.product(*listoftuples * n)))


Curve1 = Curve("alpha", 1.2*sp.exp(t));
Curve2 = Curve("beta", 0.8*sp.exp(t));
Curve3 = Curve("gamma", 1*sp.exp(t));
Curves = [Curve1,Curve2,Curve3]

Pants1 = Pants(Curve1, Curve2, Curve3)
Pants2 = Pants(Curve1, Curve2, Curve3)

if os.path.isfile('germs.json'):
    GERMS = pd.read_json('germs.json')
else:
    GERMS = pd.DataFrame(columns = ["lamination", "twist_a", "twist_b", "twist_c", "CR-lam"])

    # This bit of code deals with the 3-symmetric laminations:
    triang1 = "3"
    triang2 = "3"
    base1 = Curve1
    base2 = Curve2
    lam1 = triang1 + base1.get_name()
    lam2 = triang2 + base2.get_name()
    for dirs1 in comb([(-1,1)],3):
        for dirs2 in comb([(-1,1)],3):
            Pants1.set_triangulation(base1, triang1, dict(zip(Curves,dirs1)))
            Pants2.set_triangulation(base2, triang2, dict(zip(Curves,dirs2)))
            CR=0
            if Pants1.get_twists() == Pants2.get_twists():
                CR = 1
            germ1 = germ(Pants1,Pants2,Curve1)
            germ2 = germ(Pants1,Pants2,Curve2)
            germ3 = germ(Pants1,Pants2,Curve3)
            lam = lam1 + "," + lam2
            GERMS = GERMS.append({"lamination" : lam, "twist_a" : germ1, "twist_b" : germ2, "twist_c" : germ3, "CR-lam": CR}, ignore_index=True)

    # This bit of code deals with the 2-symmetric laminations
    for triang1, triang2 in [("2","2"),("2","3"),("3","2")]:
        print("working on: Pants1: " + triang1 +"symmetric triangulation"\
                +" and Pants2: " + triang2 + "symmetric triangulation")
        for base1 in Curves:
            for base2 in Curves:
                lam1 = triang1 + base1.get_name()
                lam2 = triang2 + base2.get_name()
                for dirs1 in comb([(-1,1)],3):
                    for dirs2 in comb([(-1,1)],3):
                        Pants1.set_triangulation(base1, triang1, dict(zip(Curves,dirs1)))
                        Pants2.set_triangulation(base2, triang2, dict(zip(Curves,dirs2)))
                        CR=0
                        if Pants1.get_twists() == Pants2.get_twists():
                            CR = 0.5
                            if lam1 == lam2:
                                CR = 1
                        germ1 = germ(Pants1,Pants2,Curve1)
                        germ2 = germ(Pants1,Pants2,Curve2)
                        germ3 = germ(Pants1,Pants2,Curve3)
                        lam = lam1 + "," + lam2
                        GERMS = GERMS.append({"lamination" : lam, "twist_a" : germ1, "twist_b" : germ2, "twist_c" : germ3, "CR-lam": CR}, ignore_index=True)
    GERMS.to_json('germs.json')

#GERMS.drop_duplicates(subset=['twist_a', 'twist_b', 'twist_c'])
#CH.drop_duplicates(subset=['twist_a', 'twist_b', 'twist_c'])
TWISTMATCH = GERMS[GERMS['CR-lam'] > 0]
CR = GERMS[GERMS['CR-lam'] > .5]

#print(len(my_df))
#print(my_df)
#visualize3DData(CR)
plot_convex_hull(GERMS)
