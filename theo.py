"""
Created on Thu July 13 10:30 2017

@author: irvin.ding

"""

#!/usr/bin/python3
import os
import random
from data_reader import data_extract
from OpenCLGA import SimpleGene, SimpleChromosome, utils, OpenCLGA

# data extraction
data = data_extract('data_test_724_1.csv')
n = len(data)
m = len(data[0])
data_size = n*m

# Input Arguments
# put data into individual arrays
# change if have more inputs
theo = []
for i in range(0, n):
    for j in range(0, m):
        theo.append(data[i][j][0])
        
bid = []
for i in range(0, n):
    for j in range(0, m):
        bid.append(data[i][j][1])
        
bidv = []
for i in range(0, n):
    for j in range(0, m):
        bidv.append(data[i][j][2])
        
ask = []
for i in range(0, n):
    for j in range(0, m):
        ask.append(data[i][j][3])
        
askv = []
for i in range(0, n):
    for j in range(0, m):
        askv.append(data[i][j][4])
        
# info print
def show_generation_info(index, data_dict):
    print('{0}\t\t==> {1}'.format(index, data_dict['best']))

def run(num_chromosomes, generations):
    random.seed()
    
    # Operations and Terminals
    terminals = [('bid', 0), ('ask', 0), ('0.5', 0), ('bidv',0), ('askv',0)]
    functions1 = [('^2', 1), ('neg', 1)]
    functions2 = [('+',2), ('-',2), ('*',2), ('/',2)]
    
    # Depth of the tree
    min_depth = 4
    max_depth = 8
    depth = random.randint(min_depth, max_depth)
    
    # Tree Initialization (should not require modification if the data structure is not changed)
    i = 0
    j = 0
    n = 1
    formula = []
    setList = []
    temp = []
    tempSetList = []
    numNodes = n
    for i in range(depth):
        while 1:
            tmp_choice = random.randint(0,2)
            if tmp_choice == 0:
                k = random.randint(0, len(terminals)-1)
                temp.append((terminals[k], terminals, terminals[k][0]))
                tempSetList.append(tmp_choice)
                j += 0
            if tmp_choice == 1:
                k = random.randint(0, len(functions1)-1)
                temp.append((functions1[k], functions1, functions1[k][0]))
                tempSetList.append(tmp_choice)
                j += 1
            if tmp_choice == 2:
                k = random.randint(0, len(functions2)-1)
                temp.append((functions2[k], functions2, functions2[k][0]))
                tempSetList.append(tmp_choice)
                j += 2
            if len(temp) > n:
                temp.clear()
                tempSetList.clear()
                j = 0
            if len(temp) == n:
                if (j > 0 & i+1 == depth):
                    temp.clear()
                    tempSetList.clear()
                    j = 0
                if (j == 0 & i+1 < depth):
                    temp.clear()
                    tempSetList.clear()
                    j = 0
                else:
                    break
        formula.extend(temp)
        setList.extend(tempSetList)
        temp = []
        tempSetList = []
        n = j
        j = 0
        numNodes = numNodes + n
        i += 1

    for i in range(0, n):
        k = random.randint(0, len(terminals)-1)
        formula.append((terminals[k], terminals, terminals[k][0]))
        setList.append(0)
    #tree_size = numNodes;        
        
    # GP run
    # sample chromosome - which defines the format of following generations
    # first arg - gene, second arg - mutation set
    sample = SimpleChromosome([SimpleGene(node[0], node[1], node[2]) for node in formula])

    # path of the cl file
    self_path = os.path.dirname(os.path.abspath(__file__))
    f = open(os.path.join(self_path, 'theo.cl'), 'r')
    fstr = ''.join(f.readlines())
    f.close()

    import threading
    evt = threading.Event()
    evt.clear()
    def state_changed(state):
        if 'stopped' == state:
            evt.set()

    ga_cl = OpenCLGA({'sample_chromosome': sample,
                      'termination': { 'type': 'count',
                                       'count': generations },
                      'population': num_chromosomes,
                      'fitness_kernel_str': fstr,
                      'fitness_func': 'theo_fitness',   # match this with the fitness func name in cl file
                      'fitness_args': [{'t': 'float', 'v': theo, 'n': 'theo'}, #change here when have more inputs
                                       {'t': 'float', 'v': bid, 'n': 'bid'},
                                       {'t': 'float', 'v': ask, 'n': 'ask'},
                                       {'t': 'float', 'v': bidv, 'n': 'bidv'},
                                       {'t': 'float', 'v': askv, 'n': 'askv'},
                                       {'t': 'int', 'v': data_size, 'n': 'data_size'},
                                       {'t': 'int', 'v': numNodes, 'n':'tree_size'},
                                       {'t': 'int', 'v': setList, 'n':'setList'}],
                      'opt_for_max': 'min',
                      #'elitism_mode': {'top': 10,
                      #                 'every': 10},
                      'debug': False,
                      'generation_callback': show_generation_info},
                      action_callbacks = {'state' : state_changed})
    
    ga_cl.prepare()

    # set up the mutation/crossover probabilities
    prob_mutate = 0.6
    prob_cross = 0.3
    ga_cl.run(prob_mutate, prob_cross)
    evt.wait()

    print('run took', ga_cl.elapsed_time, 'seconds')
    best_chromosome, best_fitness, best_info = ga_cl.get_the_best()
    print(' '.join(str(g.dna) for g in best_info.genes))
    print('Fitness: %f'%(best_fitness))

    # test data output chunk, can be removed
    f = open('testdata.txt','a')
    value = str(best_fitness)
    time = str(ga_cl.elapsed_time)
    f.write(value + ',' + time + '\n')
    f.close()


if __name__ == '__main__':
    # main function call, iterations optional
    for i in range(5):
        run(num_chromosomes=1000, generations=1000)