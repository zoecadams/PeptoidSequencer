from openpyxl import Workbook, load_workbook
import numpy as np
#from datascience import *

wb = load_workbook('DowSeq.xlsx', data_only=True)
seq = wb['Sequences(edit)']
mol_ion_mass = seq['B9':'B218']
frag = wb['Fragments(edit)']
etc = dict(Proton=1.01, Linker=326.21, Acetyl=43.02, Fluoro=388.08, MetAhx=213.12, Ahx=113.08)
side_chains = dict(Gly=115.03, Ala=129.04, Amb=143.06, EDA=100.06, Pip=191.06, But=113.09, Ben=147.07)
#[115.03, 129.04, ]

def mol_ion(mass, acetyl):
    mol_ion_mass_mod = list(mol_ion_mass)
    for obj in mol_ion_mass_mod:
        mz = obj[0]
        if int(mass) == int(mz.value + etc['Acetyl']):
            mass -= etc['Acetyl'] + etc['Proton'] + etc['Linker']#Just use int since our mass spec isn't super accurate so estimate/round to allow for error in measurement
            linear_combination(int(mass))
        if int(mass) == int(mz.value + etc['Fluoro']):
            mass -= etc['Fluoro'] + etc['Proton'] + etc['Linker']
            linear_combination(int(mass))
        if int(mass + 1) == int(mz.value + etc['Acetyl']):
            mass -= etc['Acetyl'] + etc['Proton'] + etc['Linker']
            linear_combination(int(mass))
        if int(mass + 1) == int(mz.value + etc['Fluoro']):
            mass -= etc['Fluoro'] + etc['Proton'] + etc['Linker']
            linear_combination(int(mass))
        if int(mass - 1) == int(mz.value + etc['Acetyl']):
            mass -= etc['Acetyl'] + etc['Proton'] + etc['Linker']
            linear_combination(int(mass))
        if int(mass - 1) == int(mz.value + etc['Fluoro']):
            mass -= etc['Fluoro'] + etc['Proton'] + etc['Linker']
            linear_combination(int(mass))


def factors_set():
    factors_set = ( (i, j , k ,l, m, n, p) for i in [0,1,2,3,4]
                          for j in [0,1,2,3,4]
                          for k in [0,1,2,3,4]
                          for l in [0,1,2,3,4]
                          for m in [0,1,2,3,4]
                          for n in [0,1,2,3,4]
                          for p in [0,1,2,3,4] )
    for factor in factors_set:
        yield factor
def memoize(f):
    results = {}
    def helper(x):
        if x not in results:
            results[x] = f(x)
        return results[x]
    return helper

def linear_combination(mass):
    """ returns the tuple (i,j,k,l) satisfying
        n = i*1 + j*3 + k*9 + l*27      """
    weighs = tuple(side_chains.values())
    for factors in factors_set():
       total = 0
       sequences = []
       sidechains = []
       list_side_chains = [list(side_chains.keys())]
       list_side_chains2 = np.transpose(list_side_chains)
       for i in range(len(factors)):
          total += factors[i] * weighs[i]
       if int(total) == mass or mass <= int(total) + 2 and mass >= int(total) - 2 :
          if sum(list(factors)) == 4:
            list_factors = [list(factors)]
            list_factors2 = np.transpose(list_factors)
            array1 = np.column_stack((list_side_chains2, list_factors2))

            for i in range(0,int(array1[0,1])):
              sidechains.append(array1[0,0])
            for i in range(0,int(array1[1,1])):
              sidechains.append(array1[1,0])
            for i in range(0,int(array1[2,1])):
              sidechains.append(array1[2,0])
            for i in range(0,int(array1[3,1])):
              sidechains.append(array1[3,0])
            for i in range(0,int(array1[4,1])):
              sidechains.append(array1[4,0])
            for i in range(0,int(array1[5,1])):
              sidechains.append(array1[5,0])
            for i in range(0,int(array1[6,1])):
              sidechains.append(array1[6,0])
          if end_res in sidechains:
              perms = all_perms(sidechains)
              listperms = list(perms)
              filtered_listperms = filter1(listperms,end_res)
              # print(filtered_listperms)
              ions = ion_finder(filtered_listperms)


def all_perms(elements):
    if len(elements) <=1:
        yield elements
    else:
        for perm in all_perms(elements[1:]):
            for i in range(len(elements)):
                # nb elements[0:1] works in both string and list contexts
                yield perm[:i] + elements[0:1] + perm[i:]

def filter1(lst, criteria):
    filterperms = []
    for i in range(0,len(lst)):
        if lst[i][3] == criteria:
          filterperms.append(lst[i])
    return filterperms
def ion_finder(filtered_listperms1):
  #determine possible y type ions
  ys = []
  for  q in range(0,len(filtered_listperms1)):
    y1 = round(etc['Proton']*2 + etc['MetAhx'],2)
    y2 = round(etc['Proton']*2 + etc['MetAhx']+etc['Ahx'],2)
    y3 = round(etc['Proton']*2 + etc['MetAhx']+etc['Ahx']+side_chains[filtered_listperms1[q][0]],2)
    y4 = round(etc['Proton']*2 + etc['MetAhx']+etc['Ahx']+side_chains[filtered_listperms1[q][0]]+side_chains[filtered_listperms1[q][1]],2)
    y5 = round(etc['Proton']*2 + etc['MetAhx']+etc['Ahx']+side_chains[filtered_listperms1[q][0]]+side_chains[filtered_listperms1[q][1]]+side_chains[filtered_listperms1[q][2]],2)
    y6 = round(etc['Proton']*2 + etc['MetAhx']+etc['Ahx']+side_chains[filtered_listperms1[q][0]]+side_chains[filtered_listperms1[q][1]]+side_chains[filtered_listperms1[q][2]]+side_chains[filtered_listperms1[q][3]],2)
  
    #print('y:', q,  y3, y4, y5)
    ys.append([y3,y4,y5])
  print ('y:', ys)

  #determine possible b type ions
  bs = []
  for  r in range(0,len(filtered_listperms1)):
    b1 = round(etc['Acetyl'],2)
    b2 = round(etc['Acetyl']+side_chains[filtered_listperms1[r][3]],2)
    b3 = round(etc['Acetyl']+side_chains[filtered_listperms1[r][3]]+side_chains[filtered_listperms1[r][2]],2)
    b4 = round(etc['Acetyl']+side_chains[filtered_listperms1[r][3]]+side_chains[filtered_listperms1[r][2]]+side_chains[filtered_listperms1[r][1]],2)
    b5 = round(etc['Acetyl']+side_chains[filtered_listperms1[r][3]]+side_chains[filtered_listperms1[r][2]]+side_chains[filtered_listperms1[r][1]]+side_chains[filtered_listperms1[r][0]],2)
    b6 = round(etc['Acetyl']+side_chains[filtered_listperms1[r][3]]+side_chains[filtered_listperms1[r][2]]+side_chains[filtered_listperms1[r][1]]+side_chains[filtered_listperms1[r][0]]+etc['Ahx'],2)
    bs.append([b2,b3,b4])
    #print('b:',r, b2, b3, b4)
  print('b:', bs)

mass = 827#input('Input Mass:')
acetyl = True#input('Acetylated? y or n     ')
y=True
n=False
end_res = 'Ala'#input('What is the last residue?   ')






### DON'T DEFINE ANY FUNCTIONS BELOW THIS POINT ###
mol_ion(int(mass),acetyl)


              









'''

EDA == side_chains.keys(EDA)
Gly == side_chains.keys(Gly)
Ala == side_chains.keys(Ala)
Pip == side_chains.keys(Pip)
Amb == side_chains.keys(Amb)
But == side_chains.keys(But)
Ben == side_chains.keys(Ben)

def sequence(mass):
    def sequence_helper(mass, monomer, monomer_list, poss_seq):
        if mass >= -1 and mass <= 1:
            return poss_seq
        if mass < -1:
            return 0
        if mass > 1:
            for i in range(0, len(monomer_list)):
                without_monomer = sequencer_helper(mass, monomer_list[i+1], monomer_list, poss_seq)
                with_monomer = sequence_helper(mass-monomer_list[i], monomer_list[i], monomer_list, poss_seq)
    return sequencer_helper(mass, side_chains.values()[0], side_chains.values(), [0, 0, 0, 0, 0, 0, 0])
'''
