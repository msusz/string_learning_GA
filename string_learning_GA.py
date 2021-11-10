def encode_char(char):
    """
    encode_char - function encodes the char into binary number
    based on position of the char in the global ALPHABET variable

    :param char: char that needs to be encoded
    :return: list of values for binary representation of the char in the alphabet
    """ 
    alphabet_list = list(ALPHABET)
    gene = format(alphabet_list.index(char), "b")
    gene = (6-gene.__len__())*"0"+gene
    gene = list(gene)
    gene_o = []
    for i in gene:
        i=int(i)
        gene_o.append(i)
    return gene_o

def encode_string(string):
    """
    encode_string - function encodes the string into list of binary numbers
    based on position of the char in the global ALPHABET variable

    :param string: string that needs to be encoded
    :return: list of lists of values for binary representation of every char from the string in the alphabet
    """ 
    chromosome = []
    for ch in string:
        gene = encode_char(ch)
        chromosome.append(gene)
    return chromosome

def bits_to_int(bits):
    """
    bits_to_int - function changes list of bits into a number

    :param bits: list of binary values
    :return: integer represented by input binary values
    """ 
    string_ints = [str(bit) for bit in bits]
    b = ''.join(string_ints)
    i = int(b, 2)
    return i

def decode_gene(gene):
    """
    decode_gene - function decodes the gene into a character

    :param gene: list of binary values to be changed into a character
    :return: decoded character
    """ 
    number = bits_to_int(gene)
    return ALPHABET[number]

def decode_chromosome(chromosome):
    """
    decode_chromosome - function decodes the chromosome into a string of characters

    :param chromosome: list of lists of binary values to be changed into a string of characters
    :return: decoded string of characters
    """ 
    decoded = []
    for gene in chromosome:
        decoded.append(decode_gene(gene))
    return ''.join(decoded)

def selection(population):
    """
    selection - function selects top chromosomes out of the population
    based on selection rate (a fraction of the population that should stay)

    :param population: list of chromosomes to choose from
    :return: selected part of the population
    """ 
    global SELECTION_RATE
    global POPULATION_SIZE
    number_kept = ceil(SELECTION_RATE*POPULATION_SIZE)
    selected_population = population[:number_kept]
    return selected_population

def mate(chromosome1, chromosome2):
    """
    mate - function performs mating on the chromosomes,
    randomly choosing from which chromosome it should take every bit

    :param chromosome1: first parent for mating
    :param chromosome2: second parent for mating
    :return: child consisting bits from parents
    """ 
    child=[]
    for gp1, gp2 in zip(chromosome1, chromosome2):
        prob = random()
        if prob<0.5:
            child.append(gp1)
        else:
            child.append(gp2)
    return child

def mutate(pop):
    """
    mutate - function performs mutation on the population
    randomly choosing for every chromosome if it should be mutated
    and which bit should be mutated, then changing the chosen bit
    to the opposite one

    :param pop: population that will be mutated
    :return: mutated population 
    """ 
    mutated = []
    for ch in pop:
        c = random()
        new_entity = deepcopy(ch)
        if c<MUTATION_RATE:
            nr_of_gene = randint(0, N_BITS-1)
            new_entity[nr_of_gene//6][nr_of_gene%6] = abs(new_entity[nr_of_gene//6][nr_of_gene%6]-1)
            print("MUTATION!!!!!")
        mutated.append(new_entity)
    return mutated

def create_population():
    """
    create_population - function creates new population (list of chromosomes)
    based on globally specified population size

    :return: created population 
    """ 
    global TARGET
    global POPULATION_SIZE
    chromosome_len = len(TARGET)
    population = []
    for i in range (0, POPULATION_SIZE):
        chromosome = []
        for j in range (0, chromosome_len):
            chromosome.append(choice(ALPHABET))
        encoded = encode_string(chromosome)
        population.append(encoded)
    return population

def calculate_cost(chromosome):
    """
    calculate_cost - function calculates cost for the chromosome
    by counting how many bits are different than in the target string

    :param chromosome: chromosome (list of bits) to calculate cost from
    :return: number of wrong bits (cost of the chromosome)
    """ 
    global ENCODED_TARGET
    flat_target = [item for sublist in ENCODED_TARGET for item in sublist]
    flat_chromosome = [item for sublist in chromosome for item in sublist]
    cost = 0
    for gs, gt in zip(flat_chromosome, flat_target):
        if gs != gt: cost+=1
    return cost

from random import random, choice, randint
from math import ceil
from copy import deepcopy


POPULATION_SIZE = 200
SELECTION_RATE = 0.1
MUTATION_RATE = 0.01
ALPHABET = "abcdefghijklmnopqrstuvwxyz0123456789.,;:?!_+-*/ " + 16*"a"
TARGET = "marcelina suszczyk"
ENCODED_TARGET = encode_string(TARGET)
N_BITS = len(ENCODED_TARGET)*len(ENCODED_TARGET[0])

NUMBER_OF_ITERATIONS = 300

def main():
    generation = 1
    found = False
    population = create_population()
    
    while not found or generation<NUMBER_OF_ITERATIONS:
        population = sorted(population, key = lambda x:calculate_cost(x))
        
        if calculate_cost(population[0])<=0:
            found = True
            break
        
        new_generation = selection(population)

        children =[]
        for _ in range(POPULATION_SIZE - len(new_generation)): 
            parent1 = choice(new_generation) 
            parent2 = choice(new_generation) 
            child = mate(parent1, parent2)
            children.append(child) 
        
        children = mutate(children)
        new_generation.extend(children)
        population = new_generation
        
        print("Gen: {}\tSolution: {}\tFitness Score: {}".
              format(generation, 
              "".join(decode_chromosome(population[0])), 
              calculate_cost(population[0]))) 
  
        generation += 1

    print("Gen: {}\tSolution: {}\tFitness Score: {}".
          format(generation, 
          "".join(decode_chromosome(population[0])), 
          calculate_cost(population[0]))) 
    if found:
        print("Target was found in", generation, "iterations.")
    else:
        print("Target was not found.")

if __name__ == '__main__': 
    main()
