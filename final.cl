#define OPTIMIZATION_FOR_MAX 0

#define POPULATION_SIZE 1000
#define CHROMOSOME_TYPE __SimpleChromosome

#define CHROMOSOME_SIZE SIMPLE_CHROMOSOME_GENE_SIZE
#define CALCULATE_FITNESS theo_fitness
#define FITNESS_ARGS , global float* _f_theo, global float* _f_bid, global float* _f_ask, global float* _f_bidv, global float* _f_askv, global int* _f_data_size, global int* _f_tree_size, global int* _f_setList
#define FITNESS_ARGV , _f_theo, _f_bid, _f_ask, _f_bidv, _f_askv, _f_data_size, _f_tree_size, _f_setList

#define SIMPLE_CHROMOSOME_GENE_ELEMENTS_SIZE {2, 4, 5, 4, 2, 4, 5, 2, 4, 2, 4, 2, 5, 5, 5, 5}
#define SIMPLE_CHROMOSOME_GENE_SIZE 16
#define SIMPLE_CHROMOSOME_GENE_MUTATE_FUNC simple_gene_mutate

#include "simple_gene.cl"
#include "simple_chromosome.cl"


/*
Created on Thu July 13 10:30 2017

@author: irvin.ding
*/

/* any computational operation should be defined here */
float add(float x, float y)
{
    return (x + y);
}

float minus(float x, float y)
{
    return (x - y);
}

float multiply(float x, float y)
{
    return (x * y);
}

float divide(float x, float y)
{
    if (x == 0 || y == 0)
        return 0;
    else
        return (x / y);
}

float square(float x)
{
    return (x * x);
}

float neg(float x)
{
    return ((-1) * x);
}
    
/* translate the genes into tuples (nodes) */
void build_root(__global int* genes, __global int* setList, int* tuples, int* parent, float* descendants, int size)
{
    /* initialize all place to be -1 */
    for (int p = 0; p < 2*size; p++)
    {
        tuples[p] = -1;
        parent[p] = -1;
        descendants[p] = -1.0;
    }
    for (int i = 0; i < size; i++)
    {
        tuples[2*i] = genes[i];
        tuples[2*i+1] = setList[i];
    }
    return;
}

/* build a tree structure from the nodes */
/* key observation: store parent-child info both ways */
void tree(__global int* genes, int* tuples, int* parent, float* descendants, int size)
{
    int current = 1;
    for (int i = 0; i < size; i++)  /* note that openCL does not support recursion so iteration is implemented */
    {
        if (tuples[2*i+1] == 1)
        {
            descendants[2*i] = 0;
            parent[current*2] = i;  /*store parent label in the genes (real index)*/
            parent[current*2+1] = 1;   /* store order as descendant */
            current += 1;
        }
        if (tuples[2*i+1] == 2)
        {
            descendants[2*i] = 0;
            descendants[2*i+1] = 0;
            parent[current*2] = i;
            parent[current*2+1] = 1;
            current += 1;
            parent[current*2] = i;
            parent[current*2+1] = 2;
            current += 1;
        }
    }
    return;
}

/* (operation, type of node); (parent, place in parent); (descendant1, descendant2). */
float eval(int* tuples, int* parent, float* descendants, float bid, float ask, float bidv, float askv, int size)
{
    for (int i = 0; i < size; i++)
    {
        int idx = size-1-i;
        int operation = tuples[idx*2+1];
        int tp = tuples[idx*2];
        int pt = parent[idx*2];
        int plc = parent[idx*2+1];
        switch(operation){
        /* change here when have more inputs */
            case 0:  /* all terminals */
                switch(tp){
                    case 0:
                        {if (idx == 0)
                            return bid;
                        else
                            descendants[pt*2+plc-1] = bid;}
                        break;
                    case 1:
                        {if (idx == 0)
                            return ask;
                        else
                            descendants[pt*2+plc-1] = ask;}
                        break;
                    case 2:
                        {if (idx == 0)
                            return 0.5;
                        else
                            descendants[pt*2+plc-1] = 0.5;}
                        break;
                    case 3:
                        {if (idx == 0)
                            return bidv;
                        else
                            descendants[pt*2+plc-1] = bidv;}
                        break;
                    case 4:
                        {if (idx == 0)
                            return askv;
                        else
                            descendants[pt*2+plc-1] = askv;}
                        break;
                }
                break;
            case 1:  /* all single argument operations */
                switch(tp){
                    case 0:
                        {if (idx == 0)
                            return (square(descendants[idx*2]));
                        else
                            descendants[pt*2+plc-1] = square(descendants[idx*2]);}
                        break;
                    case 1:
                        {if (idx == 0)
                            return (neg(descendants[idx*2]));
                        else
                            descendants[pt*2+plc-1] = neg(descendants[idx*2]);}
                        break;
                }
                break;
            case 2:  /* all double argument operations */
                switch(tp){
                    case 0:
                        {if (idx == 0)
                            return (add(descendants[idx*2], descendants[idx*2+1]));
                        else
                            descendants[pt*2+plc-1] = add(descendants[idx*2], descendants[idx*2+1]);}
                        break;
                    case 1:
                        {if (idx == 0)
                            return (minus(descendants[idx*2], descendants[idx*2+1]));
                        else
                            descendants[pt*2+plc-1] = minus(descendants[idx*2], descendants[idx*2+1]);}
                        break;
                    case 2:
                        {if (idx == 0)
                            return (multiply(descendants[idx*2], descendants[idx*2+1]));
                        else
                            descendants[pt*2+plc-1] = multiply(descendants[idx*2], descendants[idx*2+1]);}
                        break;
                    case 3:
                        {if (idx == 0)
                            return (divide(descendants[idx*2], descendants[idx*2+1]));
                        else
                            descendants[pt*2+plc-1] = divide(descendants[idx*2], descendants[idx*2+1]);}
                        break;
                }
                break;
            default:
                return -1.0;
        }
    }
}


void theo_fitness(global __SimpleChromosome* chromosome,
                  global float* fitnesses,
                  int chromosome_size,
                  int chromosome_count,
                  global float* theo,  /*change here when have more inputs*/
                  global float* bid,   /* these should correspond to the fitness_args from python interface */
                  global float* ask,
                  global float* bidv,
                  global float* askv,
                  global int* data_size,
                  global int* tree_size,
                  global int* setList)
{
    int length = data_size[0];
    int size = tree_size[0];
    float val = 0.0;
    float err = 0.0;
    int tuples[2*8191];
    int parent[2*8191];
    float descendants[2*8191];
    build_root(chromosome->genes, setList, tuples, parent, descendants, size);
    tree(chromosome->genes, tuples, parent, descendants, size);
    for (int i = 0; i < length; i++)
    {
        val = eval(tuples, parent, descendants, bid[i], ask[i], bidv[i], askv[i], size); /* change here when have more inputs */
        err += (val-theo[i])*(val-theo[i]);
    }
    *fitnesses = err;
    return;
    
}













































#include "ga_utils.cl"

/**
 * ocl_ga_calculate_fitness is a wrapper function to call specified fitness
 * function. To simplify the implementation of fitness calcuation, we handles
 * the multi-threading part and let implementor write the core calculation.
 * Note: this is a kernel function and will be called by python.
 *
 * @param *chromosomes (global) the chromosomes array
 * @param *fitness (global) the fitness array for each chromosomes.
 * @param FITNESS_ARGS the fitness arguments from ocl_ga options.
 */
__kernel void ocl_ga_calculate_fitness(global int* chromosomes,
                                       global float* fitness FITNESS_ARGS)
{
  int idx = get_global_id(0);
  // out of bound kernel task for padding
  if (idx >= POPULATION_SIZE) {
    return;
  }
  // calls the fitness function specified by user and gives the chromosome for
  // current thread.
  CALCULATE_FITNESS(((global CHROMOSOME_TYPE*) chromosomes) + idx,
                    fitness + idx,
                    CHROMOSOME_SIZE, POPULATION_SIZE FITNESS_ARGV);
}
