#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 11:33:44 2020

@author: emmanuel
"""


from bayes_opt import BayesianOptimization
from bayes_opt import SequentialDomainReductionTransformer
from bayes_opt import UtilityFunction
import numpy as np
import matplotlib.pyplot as plt

def cost_camille (x, y):
    return -(x+y-2)**2-x**2


def cost_camille_list(next_points):
    targets=[]
    for next_point in next_points:
        targets.append(cost_camille(**next_point))
    return targets

pbounds = {'x': (-5, 5), 'y': (-5, 5)}
bounds_transformer = SequentialDomainReductionTransformer()
## Reduces the bounds ###
optimizer = BayesianOptimization(f = cost_camille, pbounds=pbounds,
                                 bounds_transformer=bounds_transformer)




N_parallel=5


def random_kappa(x0,sigma):
    kappa=sigma*np.random.randn()+x0
    while kappa<0.:
        kappa=sigma*np.random.randn()+x0
    return kappa

points_list=[]
for _ in range(10):
    next_points=[]
    for i in range(N_parallel):
        kappa=random_kappa(3,2)
        print('random kappa =',kappa)
        utility = UtilityFunction(kind="ucb", kappa=kappa, xi=0.0)

        suggested_point=optimizer.suggest(utility)
        if suggested_point not in points_list:
            next_points.append(suggested_point)
            points_list.append(suggested_point)

    targets=cost_camille_list(next_points)


    for next_point,target in zip(next_points,targets):
        optimizer.register(params=next_point, target=target)

    plt.scatter(optimizer.space.params[-N_parallel:,0],optimizer.space.params[-N_parallel:,1])

    computed_points=next_points
    print()
    print(targets, next_points)

print('result')
print(optimizer.max)