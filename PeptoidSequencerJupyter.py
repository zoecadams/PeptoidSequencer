from openpyxl import Workbook, load_workbook
import numpy as np
from datascience import *
#%matplotlib inline
import matplotlib.pyplot as plots
plots.style.use('fivethirtyeight')

wb = load_workbook('DowSeq.xlsx', data_only=True)
seq = wb['Sequences(edit)']
mol_ion_mass = seq['B9':'B218']
frag = wb['Fragments(edit)']
etc = dict(Proton=1.01, Linker=326.21, Acetyl=43.02, Fluoro=388.08, MetAhx=213.12, Ahx=113.08)
side_chains = dict(Gly=115.03, Ala=129.04, Amb=143.06, EDA=100.06, Pip=191.06, But=113.09, Ben=147.07)

def mol_ion(mass, acetyl):
    end_res = str(input('What is the last residue?   '))#'Ala'#
    mol_ion_mass_mod = list(mol_ion_mass)
    for obj in mol_ion_mass_mod:
        mz = obj[0]
        if int(mass) == int(mz.value + etc['Acetyl']):
            #print('A')
            mass -= etc['Acetyl'] + etc['Proton'] + etc['Linker']#Just use int since our mass spec isn't super accurate so estimate/round to allow for error in measurement
            return linear_combination(int(mass), end_res)
        if int(mass) == int(mz.value + etc['Fluoro']):
            #print('B')
            mass -= etc['Fluoro'] + etc['Proton'] + etc['Linker']
            return linear_combination(int(mass), end_res)
        if int(mass + 1) == int(mz.value + etc['Acetyl']):
            #print('C')
            mass -= etc['Acetyl'] + etc['Proton'] + etc['Linker']
            return linear_combination(int(mass),end_res)
        if int(mass + 1) == int(mz.value + etc['Fluoro']):
            #print('D')
            mass -= etc['Fluoro'] + etc['Proton'] + etc['Linker']
            return linear_combination(int(mass),end_res)
        if int(mass - 1) == int(mz.value + etc['Acetyl']):
            #print('E')
            mass -= etc['Acetyl'] + etc['Proton'] + etc['Linker']
            return linear_combination(int(mass),end_res)
        if int(mass - 1) == int(mz.value + etc['Fluoro']):
            #print('F')
            mass -= etc['Fluoro'] + etc['Proton'] + etc['Linker']
            return linear_combination(int(mass),end_res)
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
def linear_combination(mass,end_res):
    """ returns the tuple (i,j,k,l) satisfying
        n = i*1 + j*3 + k*9 + l*27      """
    weighs = tuple(side_chains.values())
    ions = Table().with_columns('list','','q','','y3','','y4','','y5','','b2','','b3','','b4','')
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
            
              ions = ions.append(ion_finder(filtered_listperms))
    #print(ions.exclude(0))
    return ions.exclude(0)

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
  #determine possible b type ions
    qs = make_array()
    y3s = make_array()
    y4s = make_array()
    y5s = make_array()
    b2s = make_array()
    b3s = make_array()
    b4s = make_array()
    for  q in range(0,len(filtered_listperms1)):
        y1 = round(etc['Proton']*2 + etc['MetAhx'],2)
        y2 = round(etc['Proton']*2 + etc['MetAhx']+etc['Ahx'],2)
        y3 = round(etc['Proton']*2 + etc['MetAhx']+etc['Ahx']+side_chains[filtered_listperms1[q][0]],2)
        y4 = round(etc['Proton']*2 + etc['MetAhx']+etc['Ahx']+side_chains[filtered_listperms1[q][0]]+side_chains[filtered_listperms1[q][1]],2)
        y5 = round(etc['Proton']*2 + etc['MetAhx']+etc['Ahx']+side_chains[filtered_listperms1[q][0]]+side_chains[filtered_listperms1[q][1]]+side_chains[filtered_listperms1[q][2]],2)
        y6 = round(etc['Proton']*2 + etc['MetAhx']+etc['Ahx']+side_chains[filtered_listperms1[q][0]]+side_chains[filtered_listperms1[q][1]]+side_chains[filtered_listperms1[q][2]]+side_chains[filtered_listperms1[q][3]],2)
        b1 = round(etc['Acetyl'],2)
        b2 = round(etc['Acetyl']+side_chains[filtered_listperms1[q][3]],2)
        b3 = round(etc['Acetyl']+side_chains[filtered_listperms1[q][3]]+side_chains[filtered_listperms1[q][2]],2)
        b4 = round(etc['Acetyl']+side_chains[filtered_listperms1[q][3]]+side_chains[filtered_listperms1[q][2]]+side_chains[filtered_listperms1[q][1]],2)
        b5 = round(etc['Acetyl']+side_chains[filtered_listperms1[q][3]]+side_chains[filtered_listperms1[q][2]]+side_chains[filtered_listperms1[q][1]]+side_chains[filtered_listperms1[q][0]],2)
        b6 = round(etc['Acetyl']+side_chains[filtered_listperms1[q][3]]+side_chains[filtered_listperms1[q][2]]+side_chains[filtered_listperms1[q][1]]+side_chains[filtered_listperms1[q][0]]+etc['Ahx'],2)

        qs = np.append(qs, q)
        y3s =np.append(y3s, y3)
        y4s = np.append(y4s, y4)
        y5s = np.append(y5s, y5)
        b2s = np.append(b2s, b2)
        b3s = np.append(b3s, b3)
        b4s = np.append(b4s, b4)
    tbl_filtered_listperms1 = Table().with_columns('listperms1',filtered_listperms1)
    str_filtered_listperms1 = tbl_filtered_listperms1.apply(str,'listperms1')
    y_b_table = Table().with_columns('list', str_filtered_listperms1,'q', qs, 'y3', y3s, 'y4', y4s, 'y5', y5s,'b2', b2s, 'b3', b3s, 'b4', b4s)
    return(y_b_table)

def select_masses_to_show(table):
    y3_vals = make_array()
    y4_vals = make_array()
    y5_vals = make_array()
    b2_vals = make_array()
    b3_vals = make_array()
    b4_vals = make_array()
    for i in np.arange(len(table.column('q'))):
        if table.column('y3').item(i) not in y3_vals:
            y3_vals = np.append(y3_vals, table.column('y3').item(i))
        if table.column('y4').item(i) not in y4_vals:
            y4_vals = np.append(y4_vals, table.column('y4').item(i))
        if table.column('y5').item(i) not in y5_vals:
            y5_vals = np.append(y5_vals, table.column('y5').item(i))
        if table.column('b2').item(i) not in b2_vals:
            b2_vals = np.append(b2_vals, table.column('b2').item(i))
        if table.column('b3').item(i) not in b3_vals:
            b3_vals = np.append(b3_vals, table.column('b3').item(i))
        if table.column('b4').item(i) not in b4_vals:
            b4_vals = np.append(b4_vals, table.column('b4').item(i))
    masses_to_show = np.append(np.append(np.append(np.append(np.append(y3_vals, y4_vals), y5_vals), b2_vals), b3_vals), b4_vals)
    return masses_to_show
def interact_with_you(show_masses1, table):
    reps = float(input('How many peaks do you want to test?'))
    print(show_masses1)
    seen_vals = make_array()
    for i in np.arange(reps):
        seen_val = float(input('For which of the masses do you see a similar peak in your spectrum? Please type the exact number.   '))
        seen_vals = np.append(seen_vals, seen_val)
    numbers = make_array()
    for j in np.arange(len(table.column('q'))):
        true_count = make_array()
        for i in np.arange(2,7):
            true = np.count_nonzero(float(table.column(i).item(j)) in seen_vals)
            true_count = np.append(true_count, true)
            number = sum(true_count)
        numbers = np.append(numbers,number)
    table_new = table.with_columns('Number of matches', numbers)
    maximum = max(table_new.column('Number of matches'))
    most_likely = table_new.where('Number of matches',maximum).column('list')
    most_likely_concise = make_array()
    for q in np.arange(len(most_likely)):
        if most_likely.item(q) not in most_likely_concise:
            most_likely_concise = np.append(most_likely_concise, most_likely.item(q))
    print('The most likely sequences are:', most_likely_concise )
def run_all():
    mass = input('Input Mass:')#827#
    acetyl = input('Acetylated? y or n     ')#True#
    y=True
    n=False
    table = mol_ion(int(mass),acetyl)
    masses_to_show = select_masses_to_show(table)
    return interact_with_you(masses_to_show, table)
run_all()