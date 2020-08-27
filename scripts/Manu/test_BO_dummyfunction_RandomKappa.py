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
    return -(x+y-2)**2-x**2-3*np.cos(5*(x+y))

x=np.linspace(-5,5,101)
y=np.linspace(-5,5,101)
xv,yv=np.meshgrid(x,y)
plt.figure()
#plt.pcolor(xv,yv,cost_camille(xv, yv))
plt.contour(xv,yv,cost_camille(xv, yv))


def cost_camille_list(next_points):
    targets=[]
    for next_point in next_points:
        targets.append(cost_camille(**next_point))
    return targets

pbounds = {'x': (-5, 5), 'y': (-5, 5)}
bounds_transformer = SequentialDomainReductionTransformer()
## Reduces the bounds ###
optimizer = BayesianOptimization(f = cost_camille, pbounds=pbounds)#,
                              #   bounds_transformer=bounds_transformer)






def random_kappa(x0,sigma):
    kappa=sigma*np.random.randn()+x0
    while kappa<0.:
        kappa=sigma*np.random.randn()+x0
    kappa=5*np.random.rand()
    return kappa


def maximize(cost_camille_list,N_parallel=7,n_iter=15,x0=4,sigma=5):
    points_list=[]
    optimax=[]
    for _ in range(n_iter):
        next_points=[]
        for i in range(N_parallel):
            kappa=random_kappa(x0,sigma)
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

        print()
        print(targets, next_points)

        print()
        print('"current_max =', optimizer.max)
        optimax.append(optimizer.max['target'])
    plt.legend()

    print('result')
    print(optimizer.max)
    return optimax


optimax=maximize(cost_camille_list)
plt.figure()
plt.plot(optimax)
























