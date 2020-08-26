#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 11:33:44 2020

@author: emmanuel
"""


from bayes_opt import BayesianOptimization
from bayes_opt import SequentialDomainReductionTransformer
from bayes_opt import UtilityFunction

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
optimizer = BayesianOptimization(f = cost_camille, pbounds=pbounds)




# optimizer.maximize(5, 20)

utility = UtilityFunction(kind="ucb", kappa=2.5, xi=0.0)


# next_points=[]
# targets=[]

# for _ in range(4):
#     next_point=optimizer.suggest(utility)
#     next_points.append(next_point)
#     targets.append(cost_camille(**next_point))

# for next_point, target in zip(next_points,targets):
#     optimizer.register(params=next_point, target=target)



# for i in range(20):
#     print(i)
#     next_point = optimizer.suggest(utility)
#     target = cost_camille(**next_point)
#     optimizer.register(params=next_point, target=target)
#     print(target, next_point)
# print(optimizer.max)

N_parallel=5

####inititialisation
next_points=[]
points_list=[]
for _ in range(N_parallel):
    next_point=optimizer.suggest(utility)
    next_points.append(next_point)
    points_list.append(next_point)


targets=cost_camille_list(next_points)

init_points=next_points
init_target=targets
print(targets, next_points)
print()
print()

next_points=[]
for _ in range(N_parallel):
    next_point=optimizer.suggest(utility)
    next_points.append(next_point)
    points_list.append(next_point)

targets=cost_camille_list(next_points)
computed_points=next_points
print(targets, next_points)
print()
print()


###feed the first round of target
for point,target in zip(init_points,init_target):
    optimizer.register(params=point, target=target)



for _ in range(10):

    next_points=[]
    for point,target in zip(computed_points,targets):
        optimizer.register(params=point, target=target)

        suggested_point=optimizer.suggest(utility)
        if suggested_point not in points_list:
            next_points.append(suggested_point)
            points_list.append(suggested_point)

    plt.scatter(optimizer.space.params[-5:,0],optimizer.space.params[-5:,1])

    targets=cost_camille_list(next_points)


    computed_points=next_points
    print('lol')
    print(targets, next_points)

print('result')

print(optimizer.max)