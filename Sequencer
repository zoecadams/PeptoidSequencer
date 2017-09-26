from openpyxl import Workbook, load_workbook

wb = load_workbook('DowSeq.xlsx', data_only=True)
seq = wb['Sequences(edit)']
mol_ion_mass = seq['B9':'B218']
frag = wb['Fragments(edit)']
etc = dict(Proton=1.01, Linker=326.21, Acetyl=43.02, Fluoro=388.08)
side_chains = dict(Gly=115.03, Ala=129.04, Amb=143.05, EDA=100.06, Pip=191.06, But=113.09, Ben=147.07)

def mol_ion(mass, acetyl):
    mol_ion_mass_mod = list(mol_ion_mass)
    for obj in mol_ion_mass_mod:
        mz = obj[0]
        if mass == int(mz.value + etc['Acetyl']) or mass == int(mz.value + etc['Fluoro']): #Just use int since our mass spec isn't super accurate so estimate/round to allow for error in measurement
            print('YAY')
        elif mass + 1 == int(mz.value + etc['Acetyl']) or mass == int(mz.value + etc['Fluoro']):
            print('NAY')
        elif mass - 1 == int(mz.value + etc['Acetyl']) or mass == int(mz.value + etc['Fluoro']):
            print('OOPS')
