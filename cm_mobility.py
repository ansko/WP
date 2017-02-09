#!/usr/bin/env python3

import math
import copy
import pprint
pprint = pprint.PrettyPrinter(indent=1).pprint

#folder = '/home/anton/Desktop/2016-11-16/pa_mob_data2/' # 160 molecules
#folder = '/home/anton/Desktop/2016-10-12/2016-11-16/very_big_poly/' # 1280 molecules
#folder = '/media/anton/Seagate Expansion Drive/Article-MMT/Cluster calculations for article/BiggerSystems/Polymer/2 - 300K relaxation/1748460/' # 343 molecules
folder = '/media/anton/Seagate Expansion Drive/Article-MMT/Cluster calculations for article/BiggerSystems/Polymer/2 - 300K relaxation/1748460/'
#folder = 'comp_mob_data2/'
masses = [ 
    0, # there is no 0 type
    1.00797,
    12.0112,
    1.00797,
    12.0112,
    15.9994,
    14.0067,
    14.0067,
    12.0112
]
masses2 = [
0,
26.9815,
24.305,
28.0855,
15.9994,
15.9994,
15.9994,
15.9994,
1.00797,
14.0067,
12.0112,
12.0112,
1.00797,
1.00797,
12.0112,
15.9994,
14.0067,
14.0067,
]
masses = masses

# прочитать датафайл и вернуть атомы и границы в списке
def read_datafile(filename):
    f = open(filename, 'r')
    flag = 0
    atoms = []
    bonds = []
    bounds = []
    atom = 0
    bond = 0
    for line in f:
        if line.endswith('xlo xhi\n'):
            line_s = line.split()
            bounds.append(float(line_s[0]))
            bounds.append(float(line_s[1]))
        if line.endswith('ylo yhi\n'):
            line_s = line.split()
            bounds.append(float(line_s[0]))
            bounds.append(float(line_s[1]))
        if line.endswith('zlo zhi\n'):
            line_s = line.split()
            bounds.append(float(line_s[0]))
            bounds.append(float(line_s[1]))
        if line.startswith('Atoms # full'):
            flag = 1
            continue
        if line.startswith('Velocities'):
            flag = 2
        if line.startswith('Bonds'):
            flag = 3
            continue
        if line.startswith('Angles'):
            break
        if flag == 1:
            line = line.split()
            atoms.append([])
            for word in line:
                try:
                    word = int(word)
                except:
                    word = float(word)
                atoms[atom].append(word)
            atom += 1
        if flag == 3:
            line = line.split()
            bonds.append([])
            for word in line:
                try:
                    word = int(word)
                except:
                    word = float(word)
                bonds[bond].append(word)
            bond += 1
    atoms.pop(0)
    atoms.pop(len(atoms) - 1)
    atoms.sort()
    bonds.pop(0)
    bonds.pop(len(bonds) - 1)
    bonds.sort()

    return (atoms, bonds, bounds)

# рассчитать смещение с учётом возможного перехода границ задаются старые
# и новые координаты чего-то одного (атома или цм); границы ячейки
def calculate_disp(old, new, bounds=None):
    displacement=[0, 0, 0]
    bounds_def = [
        3.6574177352488846e-01, 9.3767139801735425e+01,
       -1.7173796716130170e+00, 7.8533787176849529e+01,
        4.2580775568494040e+00, 4.6615800720898385e+01
    ]
    if bounds is None:
        bounds = bounds_def
    # x
    delta = new[0] - old[0]
    delta2 = bounds[1] - bounds[0] - abs(delta)
    if abs(delta) < abs(delta2):
        displacement[0] += delta
    else:
        displacement[0] = -math.copysign(delta2, delta)
    # y
    delta = new[1] - old[1]
    delta2 = bounds[3] - bounds[2] - abs(delta)
    if abs(delta) < abs(delta2):
        displacement[1] += delta
    else:
        displacement[1] = -math.copysign(delta2, delta)
    # z
    delta = new[2] - old[2]
    delta2 = bounds[5] - bounds[4] - abs(delta)
    if abs(delta) < abs(delta2):
        displacement[2] += delta
    else:
        displacement[2] = -math.copysign(delta2, delta)

    return (displacement[0], displacement[1], displacement[2])

# формируется имя файла для чтения
def make_file_name(j):
    fname = folder + 'co.' + str(j * 50000 + 0) + '.data'
    return fname

# формируется список, состоящий из списков, в которых
# лежат номера атомов, относящихся к каждой молекуле
def make_molecules():
    molecules = []
    for mol_number in range(343):
        molecule = [mol_number * 192 + j + 1 for j in range(192)] #poly
        molecules.append(molecule)
    return molecules

# рассчитать координаты центра масс группы атомов
def calculate_com(atoms, bounds):
    mass = 0
    com = [0, 0, 0] # center-of-mass coordinates
    head_atom = atoms[0]
    for i in range(1, len(atoms)):
        atom = atoms[i]
        lift_set = [-1, 0, 1]
        distance = 1000
        for x_lift in lift_set:
            for y_lift in lift_set:
                for z_lift in lift_set:
                    atom_x = atom[4] + x_lift * 90
                    atom_y = atom[5]
                    atom_z = atom[6]
        atom_type = atom[2]
        atom_mass = masses[atom_type]
        mass += atom_mass
        atom_x = atom[4]
        atom_y = atom[5]
        atom_z = atom[6]
        com[0] += atom_mass * atom_x
        com[1] += atom_mass * atom_y
        com[2] += atom_mass * atom_z
    com[0] /= mass
    com[1] /= mass
    com[2] /= mass

    return com

def find_neighbours(atoms, bonds, bounds, head_atom, molecule, depth, multips):

    old_xi = multips[0]
    old_yi = multips[1]
    old_zi = multips[2]
    lx = bounds[1] - bounds[0]
    ly = bounds[3] - bounds[2]
    lz = bounds[5] - bounds[4]

    molecule.append([head_atom[0], 
                     head_atom[1],
                     head_atom[2],
                     head_atom[3],
                     head_atom[4] + old_xi * lx,
                     head_atom[5] + old_yi * ly,
                     head_atom[6] + old_zi * lz, 
                     ])
    depth += 1
    head_atom_number = head_atom[0]
    head_x = head_atom[4]
    head_y = head_atom[5]
    head_z = head_atom[6]

    multips = [-1, 0, 1]
    min_r = 10000
    min_xi = None
    min_yi = None
    min_zi = None

    for i in range(len(bonds)):
        length = len(bonds)
        if i>=length:
            break
        bond = bonds[i]
        if bond[2] == head_atom_number:
            bonds.pop(i)
            head_atom_new = atoms[bond[3] - 1]
            for xi in multips:
                for yi in multips:
                    for zi in multips:
                        new_head_x = head_atom_new[4] + xi * lx
                        new_head_y = head_atom_new[4] + xi * ly
                        new_head_z = head_atom_new[4] + xi * lz
                        dx = new_head_x - head_x
                        dy = new_head_y - head_y
                        dz = new_head_z - head_z
                        r = dx**2 + dy**2 + dz**2
                        if r < min_r:
                            min_xi = xi
                            min_yi = yi
                            min_zi = zi
                            min_r = r
            new_multips = [multips[0] + min_xi,
                           multips[1] + min_yi,
                           multips[2] + min_zi]
            if head_atom_new is None:
                return molecule
            molecule = find_neighbours(atoms, 
                                       bonds, 
                                       bounds, 
                                       head_atom_new,
                                       molecule,
                                       depth,
                                       new_multips)
        if bond[3] == head_atom_number:
            bonds.pop(i)
            head_atom_new = atoms[bond[2] - 1]
            for xi in multips:
                for yi in multips:
                    for zi in multips:
                        new_head_x = head_atom_new[4] + xi * lx
                        new_head_y = head_atom_new[4] + xi * ly
                        new_head_z = head_atom_new[4] + xi * lz
                        dx = new_head_x - head_x
                        dy = new_head_y - head_y
                        dz = new_head_z - head_z
                        r = dx**2 + dy**2 + dz**2
                        if r < min_r:
                            min_xi = xi
                            min_yi = yi
                            min_zi = zi
                            min_r = r
            new_multips = [multips[0] + min_xi,
                           multips[1] + min_yi,
                           multips[2] + min_zi]
            if head_atom_new is None:
                return molecule
            molecule = find_neighbours(atoms, 
                                       bonds, 
                                       bounds, 
                                       head_atom_new,
                                       molecule,
                                       depth,
                                       new_multips)
    return molecule

# получает некие координаты и путём параллельного переноса
# вписывает их внутрь ячейки
def to_the_box(coords, bounds):
    lx = bounds[1] - bounds[0]
    ly = bounds[3] - bounds[2]
    lz = bounds[5] - bounds[4]
    multips = [-1, 0, 1]
    xi_act = None
    yi_act = None
    zi_act = None
    for xi in multips:
        if bounds[0] < coords[0] + xi * lx < bounds[1]:
            xi_act = xi
    for yi in multips:
        if bounds[2] < coords[1] + yi * ly < bounds[3]:
            yi_act = yi
    for zi in multips:
        if bounds[4] < coords[2] + zi * lz < bounds[5]:
            zi_act = zi

    return [coords[0] + xi_act * lx,
            coords[1] + yi_act * ly,
            coords[2] + zi_act * lz]

# получаем молекулу, в которой координаты атомов могут вылезать
# за границы ячейки
def make_molecule3(mol_num, mol_len, atoms, bounds):
    lx = bounds[1] - bounds[0]
    ly = bounds[3] - bounds[2]
    lz = bounds[5] - bounds[4]
    multips = [-1, 0, 1]
    min_r = 1000
    molecule = [atoms[mol_len * mol_num + i] for i in range(mol_len)] # for polymers
    #molecule = [atoms[1560 + mol_len * mol_num + i] for i in range(mol_len)] # for comp
    molecule1 = [molecule[0]]
    for i in range(1, len(molecule)):
        head_atom = molecule1[len(molecule1) - 1]
        head_x = head_atom[4]
        head_y = head_atom[5]
        head_z = head_atom[6]
        next_atom = molecule[i]
        next_x = next_atom[4]
        next_y = next_atom[5]
        next_z = next_atom[6]
        min_xi = None
        min_yi = None
        min_zi = None
        for xi in multips:
            for yi in multips:
                for zi in multips:
                    x = next_x + xi * lx
                    y = next_y + yi * ly
                    z = next_z + zi * lz
                    r = (head_x - x)**2 + (head_y - y)**2 + (head_z - z)**2
                    if r <= min_r:
                        min_r = r
                        min_xi = xi
                        min_yi = yi
                        min_zi = zi
        molecule1.append([
            molecule[i][0],
            molecule[i][1],
            molecule[i][2],
            molecule[i][3],
            molecule[i][4] + min_xi * lx,
            molecule[i][5] + min_yi * ly,
            molecule[i][6] + min_zi * lz
        ])
        min_r = 1000

    return molecule1

# просто подсчёт центра масс, независимо от того, какие границы
#
def com_of_molecule(molecule):
    com = [0, 0, 0]
    mass = 0

    for i in range(len(molecule)):
        atom = molecule[i]
        atom_mass = masses[atom[2]]
        x = atom[4]
        y = atom[5]
        z = atom[6]
        mass += atom_mass
        com[0] += atom_mass * x
        com[1] += atom_mass * y
        com[2] += atom_mass * z

    com[0] /= mass
    com[1] /= mass
    com[2] /= mass
    
    return com

#----------------------------------------------------------------------------
def main2(mol_number=343, mol_len=192):
#1280, 160 - poly
    coms = [[0, 0, 0] for i in range(mol_number)]
    dr2 = [[0, 0, 0] for i in range(mol_number)]
    dr2ave = [0, 0, 0]

    fname = make_file_name(1)
    (atoms, bonds, bounds) = read_datafile(fname)
    for i in range(mol_number):
        molecule = make_molecule3(i, mol_len, atoms, bounds)
        com = com_of_molecule(molecule)
        com = to_the_box(com, bounds)
        coms[i] = com

    for filenum in range(1, 51):
        fname = make_file_name(filenum)
        (atoms, bonds, bounds) = read_datafile(fname)
        for i in range(mol_number):
            molecule = make_molecule3(i, mol_len, atoms, bounds)
            com = com_of_molecule(molecule)
            com = to_the_box(com, bounds)#это не надо (надо)
            #x
            delta = coms[i][0] - com[0]
            delta2 = bounds[1] - bounds[0] - abs(delta)
            if abs(delta) < abs(delta2):
                dr2[i][0] += delta
            else:
                dr2[i][0] -= math.copysign(delta2, delta)
            #y
            delta = coms[i][1] - com[1]
            delta2 = bounds[3] - bounds[2] - abs(delta)
            if abs(delta) < abs(delta2):
                dr2[i][1] += delta
            else:
                dr2[i][1] -= math.copysign(delta2, delta)
            #z
            delta = coms[i][2] - com[2]
            delta2 = bounds[1] - bounds[0] - abs(delta)
            if abs(delta) < abs(delta2):
                dr2[i][2] += delta
            else:
                dr2[i][2] -= math.copysign(delta2, delta)
            ##
            coms[i] = com
            dr2ave[0] += dr2[i][0]
            dr2ave[1] += dr2[i][1]
            dr2ave[2] += dr2[i][2]
        print((dr2ave[0]**2/mol_number + 
               dr2ave[1]**2/mol_number +
               dr2ave[2]**2/mol_number))
        dr2ave = [0, 0, 0]
    return None

main2()
