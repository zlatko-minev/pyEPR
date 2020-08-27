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
    return -(x+y-2)**2-x**2-0*3*np.cos(5*(x+y))

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
    #kappa=5*np.sqrt(np.random.rand())
    return kappa


def maximize(cost_camille_list,N_parallel=7,n_iter=30,x0=2.5,sigma=2):
    points_list=[]
    targets_list=[]
    optimax=[]
    for k in range(n_iter):
        next_points=[]
        print('k =',k)




        for i in range(N_parallel):

            optimizer = BayesianOptimization(f = cost_camille, pbounds=pbounds, bounds_transformer=bounds_transformer)#,
            if points_list is not []:
                rand_index=np.random.permutation(len(points_list))[:int(0.8*len(points_list))]
                print(len(rand_index))
                #print(len(targets_list))
                for r in rand_index:
                    optimizer.register(params=points_list[r], target=targets_list[r])

            #print(points_array)
            print('i =',i)


            kappa=random_kappa(x0,sigma)
            print('random kappa =',kappa)
            utility = UtilityFunction(kind="ucb", kappa=kappa, xi=0.0)
            suggested_point=optimizer.suggest(utility)
            print('random kappa =',kappa)

            if suggested_point not in points_list:
                if suggested_point not in next_points:
                    next_points.append(suggested_point)

        targets=cost_camille_list(next_points)

        for target in targets:
            targets_list.append(target)
        for next_point in next_points:
            points_list.append(next_point)
        #for next_point,target in zip(next_points,targets):
        #    optimizer.register(params=next_point, target=target)
        points_array=np.array([list(point.values()) for point in points_list])
        plt.scatter(points_array[-N_parallel:,0],points_array[-N_parallel:,1])

        print()
        print(targets, next_points)

        print()
        print('"current_max =', optimizer.max)
        try:
            optimax.append(optimizer.max['target'])
        except:
            print('lol')
    plt.legend()

    print('result')
    print(optimizer.max)
    return optimax


optimax=maximize(cost_camille_list)
plt.figure()
plt.plot(optimax)


















