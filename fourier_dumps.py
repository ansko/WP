#!/usr/bin/env python3

import math

#######################################################
#
# reads dump file, return list stress_dump[stress_x][z]
#

def read_dump(datafile, mode, stress_dump):
    f=open(datafile, 'r')
    i=0

    if mode not in [2, 7]:
        print("Incorrect mode")
        return None

    for string in f:
        i += 1
        if i < 10:
            continue
        string_splitted=string.split()
        (stress_dump[int(string_splitted[0])])[0] = float(string_splitted[1]) # z
        (stress_dump[int(string_splitted[0])])[1] = float(string_splitted[2]) # stress_x

    return stress_dump

#####################################
#
# converts atomic stresses to pressure
# for comp + 2 phases
#

def pressure_two_phases(stress_dump):
    cellvolume = 310321
    softvolume = 243321

    pres_comp = pres_mmt = pres_soft = 0

    for i in range(1, 31321):
        pres_comp += stress_dump[i][1]
        if (i < 721 or
           (i > 3480 and i < 4201) or
           (i > 6960 and i < 7681) or 
           (i > 10440 and i < 11160) or
           (i > 13920 and i < 14641) or
           (i > 17420 and i < 18121) or
           (i > 20880 and i < 21601) or
           (i > 24360 and i < 25081) or (i > 27840 and i < 28561)):
            pres_mmt += stress_dump[i][1]
        else:
            pres_soft += stress_dump[i][1]

    return [pres_comp / cellvolume, 
            pres_mmt / (cellvolume - softvolume),
            pres_soft / softvolume]

def appro(pressures):
    i = a0 = a1 = a2 = b1 = b2 = 0
    suma0 = suma1 = suma2 = sumb1 = sumb2 = magnitude = 0

    period = len(pressures)
    print(period)

    for i in range(period):
        suma0 += float(pressures[i])
        suma1 += float(pressures[i]) * math.cos(2 * math.pi / period * (i + 1))
        sumb1 += float(pressures[i]) * math.sin(2 * math.pi / period * (i + 1))
        suma2 += float(pressures[i]) * math.cos(4 * math.pi / period * (i + 1))
        sumb2 += float(pressures[i]) * math.sin(4 * math.pi / period * (i + 1))

    a0 = suma0 / period
    a1 = 2 * suma1 / period
    a2 = 2 * suma2 / period
    b1 = 2 * sumb1 / period
    b2 = 2 * sumb2 / period
    magnitude = 5 * math.sqrt(a1**2 + b1**2)

    print(magnitude)
    return (a0, a1, b1)

def error(pressures, a0, a1, b1):
    error = 0
    period = len(pressures)
    for i in range(period):
        error += abs((a0 +
                      a1 * math.cos(2 * math.pi / period * (i+1)) +
                      b1 * math.sin(2 * math.pi / period * (i+1))) - pressures[i])**2
    print(math.sqrt(error/len(pressures)))
    return None

folder = '500k/'
stress_dump = [[0, 0] for i in range(31321)]
pressures_comp = []
pressures_mmt = []
pressures_soft = []
for i in range(25001):
    fname = folder + 'ALLstress.' + str(i * 100)
    stresses = read_dump(fname, 2, stress_dump)
    (pres_comp, pres_mmt, pres_soft) = pressure_two_phases(stresses)
    pressures_comp.append(pres_comp)
    pressures_mmt.append(pres_mmt)
    pressures_soft.append(pres_soft)
    print(i, pres_comp, pres_mmt, pres_soft)

(a0, a1, b1) = appro(pressures_comp)
error(pressures_comp, a0, a1, b1)
(a0, a1, b1) = appro(pressures_mmt)
error(pressures_mmt, a0, a1, b1)
(a0, a1, b1) = appro(pressures_soft)
error(pressures_soft, a0, a1, b1)
