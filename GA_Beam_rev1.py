#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 23 14:35:05 2018

@author: alvaro
"""

import numpy as np
import matplotlib.pyplot as plt

np.random.seed(29)

""" Create Individual 
##################################################################"""
class Beam:
    def __init__(self, b, h):
        self.b = b      
        self.h = h  
        self.a_s = None
        self.cost = None        
        self.mu_max = None
        
        d = self.h - 5
        
        self.a_s = 0.00579 * (self.b * d)
        
        self.cost = (self.b * self.h) * 0.0332 + self.a_s * 4.1409 #custo por centimetro quadrado
        
        self.mu_max = (1/1.4) * (0.252 * (self.b * d) *((self.h/2) - 0.1036 * d) + 43.478 * self.a_s * ((self.h/2) - 5)) #kN*cm

""" Create Population 
##################################################################"""
def create_pop(n, b_min, b_max, h_min, h_max, mu_design):
    # mu_design is the minimum accepted mu. Any beam generated with a lower value will be discarded
    pop = {}
    for i in range(n):
        mu = mu_design
        while mu <= mu_design:
            b = np.random.random() * (b_max - b_min) + b_min
            h = np.random.random() * (h_max - h_min) + h_min
            key = "0.%s"%(i)
            beam = Beam(b, h)
            mu = beam.mu_max
        pop[key] = beam
    return pop

""" Mutation 
##################################################################"""
def mutation(beam, b_min, b_max, h_min, h_max, mu_design):
    # Cuidado porque Python altera o valor contido na memória, apontado pela variável.
    mu = mu_design
    while mu <= mu_design:
        if np.random.random()<0.5:
            b = np.random.random() * (b_max - b_min) + b_min
            h = beam.h
        else:
            b = beam.b
            h = np.random.random() * (h_max - h_min) + h_min
        #key = i
        mutant = Beam(b, h)
        mu = mutant.mu_max
    
    return mutant

""" Crossover 
##################################################################"""
# Crossover will generate a offspring which sits in the line between its parents
def crossover(parent_1, parent_2, mu_design):
    mu = mu_design
    while mu <= mu_design:
        b = np.random.random() * abs(parent_1.b - parent_2.b) + min(parent_1.b, parent_2.b)
        h = np.random.random() * abs(parent_1.h - parent_2.h) + min(parent_1.h, parent_2.h)
        offspring = Beam(b, h)
        mu = offspring.mu_max
    return offspring

""" Fitness 
##################################################################"""
# Vamos ponderar custo baixo e mu alto.
# Começar só com o custo

""" Selection 
##################################################################"""
def nat_select(pop, pop_size):
    score_total = 0
    probs = {}
    for key in pop:
        score_total = score_total + pop[key].cost
    for key in pop:
        probs[key] = pop[key].cost / score_total
    while len(pop) > pop_size:
        victim_key = np.random.choice(list(pop))
        roulette = np.random.random()
        if probs[victim_key] > roulette:
            del(pop[victim_key])
    return pop


""" Converggence 
##################################################################"""
def variability(pop):
    pop_costs = []
    for key in pop:
        pop_costs.append(pop[key].cost)
    pop_costs = np.array(pop_costs)
    cost_variability = np.max(pop_costs) - np.min(pop_costs)
    return cost_variability

""" Evolution 
##################################################################"""
def evolution(n, mu_design, b_min, b_max, h_min, h_max, iter_max, mutation_prob, delta_cost):
    
    pop = create_pop(n, b_min, b_max, h_min, h_max, mu_design)
    pop_init = pop.copy()
    cost_variability = variability(pop)
    i = 0
    while (i < iter_max) and (cost_variability > delta_cost):
        i += 1
        print(i)
        for pair in range(int(len(pop)/2)):
            parent_key1 = np.random.choice(list(pop))
            parent_key2 = np.random.choice(list(pop))
            
            parent_1 = pop[parent_key1]
            parent_2 = pop[parent_key2]
            
            offspring = crossover(parent_1, parent_2, mu_design)
            key = "%s.%s"%(i, pair)
            pop[key] = offspring
            
        for key in pop:
            mut_test = np.random.random()
            if mut_test < mutation_prob:
                beam = pop[key]
                mutant = mutation(beam, b_min, b_max, h_min, h_max, mu_design)
                pop[key] = mutant
                
        pop = nat_select(pop, n)
        cost_variability = variability(pop)
        
    best_key = np.random.choice(list(pop))
    best_cost = pop[best_key].cost
    for key in pop:
        if pop[key].cost < best_cost:
            best_cost = pop[key].cost
            best_key = key
            best_beam = pop[key]
    
    return (best_key, best_cost, best_beam, pop, pop_init)

""" get individuals coordinates
##################################################################"""
def coordinates(pop):
    coords = np.zeros((len(pop), 2))
    i=0
    for key in pop:
        
        coords[i][0] = pop[key].h
        coords[i][1] = pop[key].b
        
        i += 1
    return coords

""" MAIN 
##################################################################"""

n = 500
mu_design = 2000
b_min = 20
b_max = 60
h_min = 20
h_max = 60
iter_max = 20
mutation_prob = 0.05
delta_cost = 1

(key, cost, beam, pop, pop_init) = evolution(n, mu_design, b_min, b_max, h_min, h_max, iter_max, mutation_prob, delta_cost)


""" plot curves 
##################################################################"""
fig = plt.figure()
ax = fig.add_subplot(111)

u = np.linspace(20, 60, 100)
x, y = np.meshgrid(u,u)

cost = (x * y) * 0.0332 + 0.00579 * (y * (x-5)) * 4.1409 #custo por centimetro quadrad
lvl1 = np.linspace(10, 190, 20)
cs1 = plt.contourf(x, y, cost, cmap='plasma', levels=lvl1) 

mu_max = (1/1.4) * (0.252 * (y * (x-5)) *((x/2) - 0.1036 * (x-5)) + 43.478 * 0.00579 * (y * (x-5)) * ((x/2) - 5)) #kN*cm
lvl2 = np.linspace(1000, 10000, 10)
cs2 = plt.contour(x, y, mu_max, cmap='cool', levels=lvl2) 
plt.clabel(cs2)

c_initial = coordinates(pop_init)
initial = plt.scatter(c_initial[:, 0], c_initial[:, 1], c='blue')

c_final = coordinates(pop)
final = plt.scatter(c_final[:, 0], c_final[:, 1], c='red')

c_best = np.array([beam.h, beam.b])
best = plt.scatter(c_best[0], c_best[1], c='yellow')

plt.show()














