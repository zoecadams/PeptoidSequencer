from openpyxl import Workbook, load_workbook

wb = load_workbook('DowSeq.xlsx', data_only=True)
seq = wb['Sequences(edit)']
mol_ion_mass = seq['B9':'B218']
frag = wb['Fragments(edit)']
etc = dict(Proton=1.01, Linker=326.21, Acetyl=43.02, Fluoro=388.08)
side_chains = dict(Gly=115.03, Ala=129.04, Amb=143.05, EDA=100.06, Pip=191.06, But=113.09, Ben=147.07)
[115.03, 129.04, ]
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
       for i in range(len(factors)):
          total += factors[i] * weighs[i]
       if int(total) == mass or int(total) + 1 == mass or int(total) - 1 == mass:
          if sum(list(factors)) == 4:
              print(factors)
'''
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
