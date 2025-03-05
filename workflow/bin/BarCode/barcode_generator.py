#!/usr/bin/python
# encoding: utf-8

"""
@version: v1.0
@author: Zhongping Xu
@license: Apache Licence
@contact: 1099808298@qq.com
@site: http://tiramisutes.github.io/
@software: python2.7
@time: 2020/8/5 18:46
@description：Filter fasta for run BE
"""
import argparse

print '\n\nBARCODE GENERATOR by Zhongping Xu\n\t---***---\nGroup of Cotton Genetic Improvement\nHuazhong Agricultural University\nhttp://tiramisutes.github.io/\n'
print '''This program generates barcodes of a desired length, distance, and GC content.
First, your need the base primer paired end primers.\n'''


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--project", type = str, help = "1. What is the project name? Usually is gene id")
    parser.add_argument("--primerF", type = str, help = "2. The base primer of Forward from 5' to 3': ")
    parser.add_argument("--primerR", type = str, help = "3. The base primer of Reverse from 5' to 3': ")
    parser.add_argument("--length", type = int, help = "4. Barcode length; Enter LENGTH as an integer (i.e. 4)")
    parser.add_argument("--total", type = int, help = "5. Total number of barcodes (Primer Pairs) (default is LENGTH x 5)")
    parser.add_argument("--minimum_dif", type = int, help = "6. The minimum number of different bases between barcodes (default is LENGTH/2, i.e. 7->3, 4->2)")
    parser.add_argument("--minimum_GC_per", type = int, help = "7. Desired GC content minimum range in percentages (i.e. 50 ->50%%) (default is 0)", default = 0)
    parser.add_argument("--maximum_GC_per", type = int, help = "8. Desired GC content maximum range in percentages (i.e. 50 ->50%%) (default is 100)", default = 100)
    parser.add_argument("--random", type = int, help = "9. How many attempts? The default number of random codes to test is 10000. Do not enter more than a million", default = 10000)
    parser.add_argument("--output", type = str, help = "The output file name, can be add path (default is: project + '_barcode.txt')")
    
    return parser.parse_args()


comp_barcode = ''

def main():

    args = parse_args()

    #___________________________________________________raw input
    ProjectName = args.project
    primer_name = ProjectName
    Base_primer_F = args.primerF
    Base_primer_R = args.primerR
    length = args.length
    number = args.total
    if number == '':
        number = length*5
    else:
        number = number
    diffs = args.minimum_dif
    if diffs == '':
        diffs = length//2
    else:
        diffs = diffs
    mingc = args.minimum_GC_per
    if mingc == '':
        mingc = 0
    else:
        mingc = float(mingc) / 100
    maxgc = args.maximum_GC_per
    if maxgc == '':
        maxgc = 100
    else:
        maxgc = float(maxgc) / 100
    attempts = args.random
    if attempts == '':
        attempts = 10000
    else:
        attempts = attempts
    outputname = args.output
    if outputname == '':
        outputname = primer_name + '_barcode.txt'
    else:
        outputname = outputname

    #___________________________________________________process

    # make list of the four bases
    l1 = ['a', 'c', 'g', 't']

    # initialze the barcode list
    barcode_list = []

    # initialize the first barcode
    first_barcode = []

    # prime the tested list, for future counting
    tested = []

    # function to determine GC content
    def gc_cont (bar_code):
        gc = 0.0
        for base in range(length):
            if bar_code[base] == 'c' or bar_code[base] == 'g':
                gc += 1
            else:
                gc += 0
        cont = gc / length
        return cont

    # import random module
    import random

    # make the first barcode
    # add first barcode to barcode list. This is needed for the
    # first comparison of "compare_barcode" function
    while barcode_list == []:
        for i in range(length):
            first_barcode.append(random.choice(l1))
        if gc_cont(first_barcode) <= maxgc and gc_cont(first_barcode) >= mingc:
            barcode_list.append(first_barcode)
        else:
            first_barcode = []
            

    # the barcode "cradle": a place where each barcode will sit
    barcode = []

    #___________________________________________________define functions

    # function makes the barcode
    def make_barcode(length):
        global barcode
        # empties the barcode cradle
        barcode = []
        for i in range(length):
            barcode.append(random.choice(l1))

    # barcode is tested vs the previously generated barcodes
    def compare_barcode(length, barcode_l):
        count = 0
        global barcode
        # run barcode creator
        make_barcode(length)
        # keep track of it
        tested.append(barcode)
        # testing of barcode
        if barcode not in barcode_list:
            global count_list
            count_list = []
            # compare to barcodes in list
            for bc in barcode_l:
                # matches to existing barcodes
                # are scored as points
                count = 0
                for pos in range(length):
                    if barcode[pos] == bc[pos]:
                        count += 1
                    else:
                        count += 0
                # for each barcode a list of scores is made
                count_list.append(count)
            # if the barcode has enough unique bases
            # and the proper GC content, it is added
            # to the list of good barcodes
            if max(count_list) > length-diffs:
                count_list = []
            elif gc_cont(barcode) <= maxgc and gc_cont(barcode) >= mingc:
                barcode_list.append(barcode)
                count_list = []
            else:
                count_list = []
        else:
            pass
                
    #___________________________________________________run functions
        
    # initialize count

    count_list = []

    # program stalls if too many attempts are allowed
    # and few barcodes remain to be discovered
    # this loop keeps the attempts within the range allowed
    while len(tested) < attempts:
        if len(barcode_list) < number:
            compare_barcode(length, barcode_list)
        else:
            break

    barcode_list.sort()


    print "\n\nRESULTS\n\ngood barcodes and GC content:"

    for i in barcode_list:
        print i, int(gc_cont(i)*100), '%'

    print '\nnumber of tested barcodes:'

    print len(tested)

    print '\nnumber of good barcodes:'

    print len(barcode_list)

    #count base composition in each of the barcode position
    from collections import defaultdict

    print '\nbase compositions by position'

    for pos in range(length):
        list_l = [i[pos] for i in barcode_list]
        base_count = defaultdict(int)
        for base in list_l:
            base_count[base] += 1
        print base_count.keys(), base_count.values()

    #___________________________________________________manipulate barcodes and print results
    # make a file for the barcoded primer sequence
    # open file 
    barfile = open(outputname, 'wb')

    print >> barfile, 'These are the barcodes of length '+str(length)+' with a distance of '+str(diffs)+' bases.\n' + 'All the orientation of primer in this files is 5\' to 3\'!!!\n' + 'The total creat Primer Pairs: ' + str(number) + ', which ' + str(len(barcode_list)) + ' pairs are good barcodes.' + '\n\n'

    #___________________________________________________primer sequence

    # modify the primer as desired. This seq is the
    # illumina PE primer
    # Forward prime 5' to 3'
    primerA = Base_primer_F

    # note that the adapter is originally as below
    # GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG
    # however, an A is added to allow the sequencing
    # primer to anneal

    # this is the complementary primer from illumina
    # same for regular or PE 
    # Reverse prime 5' to 3'
    primerB = Base_primer_R

    # the output will be:
    # >adA2_cccca
    # ccccaAGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG
    # >adB2_cccca
    # ACACTCTTTCCCTACACGACGCTCTTCCGATCTtggggT
    # note that barcode is in lower case

    # define adapter names
    name_rootA = primer_name + '-F_'

    name_rootB = primer_name + '-R_'

    # initialize a holder name for the complement of barcode
    #comp_barcode = ''

    # define function to derive complement of any seq
    # try/except/finally is to make sure that the 'maketrans'
    # has been imported
    def reverse_comp(seq):
        try:
            maketrans  
        except NameError:
            from string import maketrans
        finally:
            comp_table = maketrans('actg','tgac')
            global comp_barcode
            comp_barcode = seq[::-1].translate(comp_table)

    for i in barcode_list:
        j = ''.join(i)
        reverse_comp(j)
        print >> barfile, '%s%s\t%s%s\n%s%s\t%s%s\n' %  (name_rootA, j, j, primerA, name_rootB, j, comp_barcode, primerB)

    # close file or it does not get updated
    barfile.close()

#___________________________________________________log

# changes from 2.7:
# rephrase raw input queries
# default barcodes to test to 10,000
# clean up shell results

if __name__ == '__main__':
    main()
