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













































