#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 23 14:35:05 2018

@author: alvaro
"""

import numpy as np

np.random.seed(29)

""" Create Individual 
##################################################################"""
class Beam:
    def __init__(self, b, h, a_s, cost, mu_max):
        self.b = b      
        self.h = h  
        self.a_s = a_s
        self.cost = cost        
        self.mu_max = mu_max      
        
def create_beam(b, h): #entradas em centímetros
    
    d = h-5
    
    a_s = 0.00579 * (b * d)
    
    cost = (b*h) * 0.0332 + a_s * 4.1409 #custo por centimetro quadrado
    
    mu_max = (1/1.4) * (0.252 * (b * d) *((h/2) - 0.1036 * d) + 43.478 * a_s * ((h/2) - 5)) #kN*cm
    
    beam = Beam(b, h, a_s, cost, mu_max)
    
    return beam

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
            key = i
            ind = create_beam(b, h)
            mu = ind.mu_max
        pop[key] = ind
    return pop

""" Mutation 
##################################################################"""
def mutation(beam, b_min, b_max, h_min, h_max, mu_design):
    # Cuidado porque Python altera o valor contido na memória, apontado pela variável.
    mu = mu_design
    while mu <= mu_design:
        if np.random.random()<0.5:
            beam.b = np.random.random() * (b_max - b_min) + b_min
        else:
            beam.h = np.random.random() * (h_max - h_min) + h_min
        #key = i
        mu = beam.mu_max
    return beam

""" Crossover 
##################################################################"""
# Crossover will generate a offspring which sits in the line between its parents
def crossover(parent_1, parent_2, mu_design):
    mu = mu_design
    while mu <= mu_design:
        b = np.random.random() * abs(parent_1.b - parent_2.b) + min(parent_1.b, parent_2.b)
        h = np.random.random() * abs(parent_1.h - parent_2.h) + min(parent_1.h, parent_2.h)
        #key = "$iter"."_"."$id_pais"
        # iter = época
        # id_pais = índice da dupla de pais, de 0 a n/2
        offspring = create_beam(b, h)
        mu = offspring.mu_max
    return offspring

""" Fitness 
##################################################################"""
# Vamos ponderar custo baixo e mu alto.
# Começar só com o custo

""" Selection 
##################################################################"""
def nat_select(pop):
    score_total = 0
    probs = {}
    pop_size = 3 # Use "n" from create_pop
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

""" Evolution 
##################################################################"""
def evolution(n, mu_design, b_min, b_max, h_min, h_max, iter_max, mutation_prob, delta_cost):
    
    pop = create_pop(n, b_min, b_max, h_min, h_max, mu_design)
    
    i = 0
    while i < iter_max:
        i += 1
        
        for pair in range(int(len(pop)/2)):
            parent_key1 = np.random.choice(list(pop))
            parent_key2 = np.random.choice(list(pop))
            
            parent_1 = pop[parent_key1]
            parent_2 = pop[parent_key2]
            
            offspring = crossover(parent_1, parent_2, mu_design)
            key = "%s%s"%(i, np.random.random())
            pop[key] = offspring
            
        for key in pop:
            mut_test = np.random.random()
            if mut_test < mutation_prob:
                beam = pop[key]
                mutant = mutation(beam, b_min, b_max, h_min, h_max, mu_design)
                pop[key] = mutant
                
        pop = nat_select(pop)
        
    best_key = np.random.choice(list(pop))
    best_cost = pop[best_key].cost
    for key in pop:
        if pop[key].cost < best_cost:
            best_cost = pop[key].cost
            best_key = key
            best_beam = pop[key]
    
    return (best_key, best_cost, best_beam)

""" MAIN 
##################################################################"""

n = 10
mu_design = 2000
b_min = 20
b_max = 60
h_min = 20
h_max = 60
iter_max = 10
mutation_prob = 0.005
delta_cost = 0.1

(key, cost, beam) = evolution(10, 2000, 20, 60, 20, 60, 10, 0.005, 0.1)

#pop = create_pop(5, 20, 60, 20, 60, 2000)
#
#for key in pop:
#    print(pop[key].cost)
#    
#for key in probs:
#    print(probs[key])
#
#pop["a"]=1
#pop["b"]=2
#key = random.choice(pop)
