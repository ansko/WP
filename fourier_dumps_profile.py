#!/usr/bin/env python3

import math

#systemsize = 31230
systemsize = 48420

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
        (stress_dump[int(string_splitted[0])])[1] = float(string_splitted[2]) # c_all

    return stress_dump

#####################################
#
# converts atomic stresses to pressure
#

def pressure_layers(stress_dump):
    lx = 93.67633450252626
    ly = 80.20140987336211
    lz = 40.92597694615771
    cellvolume = lx * ly * lz
    softvolume = lx * ly * (lz - 8.1) # 8.1 is thickness of mmt

    pres_comp = pres_mmt = pres_soft = 0
    press1 = press2 = press3 = press4 = press5 = press6 = press7 = 0

    for i in range(1, systemsize + 1):
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

        # 5x20
        #if(stress_dump[i][0] < -67.5 or stress_dump[i][0] > -31):
        #    press2 += stress_dump[i][1]
        #elif (stress_dump[i][0] < -63.5):
        #    press3 += stress_dump[i][1]
        #elif (stress_dump[i][0] < -59.5):
        #    press4 += stress_dump[i][1]
        #elif (stress_dump[i][0] < -54.5):
        #    press5 += stress_dump[i][1]
        #elif (stress_dump[i][0] < -50.5):
        #    press6 += stress_dump[i][1]
        #elif (stress_dump[i][0] < -45):
        #    press7 += stress_dump[i][1]
        #elif (stress_dump[i][0] > -36.5):
        #    press1 += stress_dump[i][1]        

    #return [pres_comp / cellvolume, 
    #        pres_mmt / (cellvolume - softvolume),
    #        pres_soft / softvolume,
    #        press1 / (lx * ly * 5),
    #        press2 / (lx * ly * 4),
    #        press3 / (lx * ly * 4),
    #        press4 / (lx * ly * 4.5),
    #        press5 / (lx * ly * 4),
    #        press6 / (lx * ly * 4),
    #        press7 / (lx * ly * 5)]

        # 10x10
        if(stress_dump[i][0] < 1.5 or stress_dump[i][0] > 38):
            press4 += stress_dump[i][1]
        elif (stress_dump[i][0] < 6):
            press5 += stress_dump[i][1]
        elif (stress_dump[i][0] < 10.5):
            press6 += stress_dump[i][1]
        elif (stress_dump[i][0] < 15.5):
            press7 += stress_dump[i][1]
        elif (23.5 < stress_dump[i][0] < 29):
            press1 += stress_dump[i][1]
        elif (29 < stress_dump[i][0] < 33.5):
            press2 += stress_dump[i][1]
        elif (38 > stress_dump[i][0] > 33.5):
            press3 += stress_dump[i][1]

    return [pres_comp / cellvolume, 
            pres_mmt / (cellvolume - softvolume),
            pres_soft / softvolume,
            press1 / (lx * ly * 5.5),
            press2 / (lx * ly * 4.5),
            press3 / (lx * ly * 4.5),
            press4 / (lx * ly * 26.8),
            press5 / (lx * ly * 4.5),
            press6 / (lx * ly * 4.5),
            press7 / (lx * ly * 5)]

def appro(pressures):
    i = a0 = a1 = a2 = b1 = b2 = 0
    suma0 = suma1 = suma2 = sumb1 = sumb2 = magnitude = 0

    period = len(pressures)
    #print(period)

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

def main(atomsnum=systemsize, step=100):
    folder = '/home/anton/Cluster/10chains/3 - wiggle/1 - 1766011/Dumps/'
    stress_dump = [[0, 0] for i in range(atomsnum + 1)] # z, 
    pressures_comp = []
    pressures_mmt = []
    pressures_soft = []
    pressures1 = []
    pressures2 = []
    pressures3 = []
    pressures4 = []
    pressures5 = []
    pressures6 = []
    pressures7 = []
    for i in range(int(2500000 / step) + 1):
        fname = folder + 'ALLstress.' + str(i * step)
        stresses = read_dump(fname, 7, stress_dump)
        (pres_comp, pres_mmt, pres_soft, press1, press2, press3, press4, press5, press6, press7) = pressure_layers(stresses)
        pressures_comp.append(pres_comp)
        pressures_mmt.append(pres_mmt)
        pressures_soft.append(pres_soft)
        pressures1.append(press1)
        pressures2.append(press2)
        pressures3.append(press3)
        pressures4.append(press4)
        pressures5.append(press5)
        pressures6.append(press6)
        pressures7.append(press7)
        print(i, pres_comp, pres_mmt, pres_soft)

    for k in range(15):
        period = len(pressures_comp)
        print("Period: ", period)

        print("Comp")
        (a0, a1, b1) = appro(pressures_comp)
        error(pressures_comp, a0, a1, b1)
        print("MMT")
        (a0, a1, b1) = appro(pressures_mmt)
        error(pressures_mmt, a0, a1, b1)
        print("Soft")
        (a0, a1, b1) = appro(pressures_soft)
        error(pressures_soft, a0, a1, b1)
        print("L1")
        (a0, a1, b1) = appro(pressures1)
        error(pressures1, a0, a1, b1)
        print("L2")
        (a0, a1, b1) = appro(pressures2)
        error(pressures2, a0, a1, b1)
        print("L3")
        (a0, a1, b1) = appro(pressures3)
        error(pressures3, a0, a1, b1)
        print("L4")
        (a0, a1, b1) = appro(pressures4)
        error(pressures4, a0, a1, b1)
        print("L5")
        (a0, a1, b1) = appro(pressures5)
        error(pressures5, a0, a1, b1)
        print("L6")
        (a0, a1, b1) = appro(pressures6)
        error(pressures6, a0, a1, b1)
        print("L7")
        (a0, a1, b1) = appro(pressures7)
        error(pressures7, a0, a1, b1)

        for j in range(period - 1):
            pressures_comp[j] += pressures_comp[j + 1]
            pressures_comp[j] /= 2
        if len(pressures_comp) % 2 == 0:
            pressures_comp = pressures_comp[::2]
        else:
            endelement = pressures_comp[-1]
            pressures_comp = pressures_comp[::2]
            pressures_comp.append(endelement)

        for j in range(period - 1):
            pressures_mmt[j] += pressures_mmt[j + 1]
            pressures_mmt[j] /= 2
        if len(pressures_mmt) % 2 == 0:
            pressures_mmt = pressures_mmt[::2]
        else:
            endelement = pressures_mmt[-1]
            pressures_mmt = pressures_mmt[::2]
            pressures_mmt.append(endelement)

        for j in range(period - 1):
            pressures_soft[j] += pressures_soft[j + 1]
            pressures_soft[j] /= 2
        if len(pressures_soft) % 2 == 0:
            pressures_soft = pressures_soft[::2]
        else:
            endelement = pressures_soft[-1]
            pressures_soft = pressures_soft[::2]
            pressures_soft.append(endelement)

        for j in range(period - 1):
            pressures1[j] += pressures1[j + 1]
            pressures1[j] /= 2
        if len(pressures1) % 2 == 0:
            pressures1 = pressures1[::2]
        else:
            endelement = pressures1[-1]
            pressures1 = pressures1[::2]
            pressures1.append(endelement)

        for j in range(period - 1):
            pressures2[j] += pressures2[j + 1]
            pressures2[j] /= 2
        if len(pressures2) % 2 == 0:
            pressures2 = pressures2[::2]
        else:
            endelement = pressures2[-1]
            pressures2 = pressures2[::2]
            pressures2.append(endelement)

        for j in range(period - 1):
            pressures3[j] += pressures3[j + 1]
            pressures3[j] /= 2
        if len(pressures3) % 2 == 0:
            pressures3 = pressures3[::2]
        else:
            endelement = pressures3[-1]
            pressures3 = pressures3[::2]
            pressures3.append(endelement)

        for j in range(period - 1):
            pressures4[j] += pressures4[j + 1]
            pressures4[j] /= 2
        if len(pressures4) % 2 == 0:
            pressures4 = pressures4[::2]
        else:
            endelement = pressures4[-1]
            pressures4 = pressures4[::2]
            pressures4.append(endelement)

        for j in range(period - 1):
            pressures5[j] += pressures5[j + 1]
            pressures5[j] /= 2
        if len(pressures5) % 2 == 0:
            pressures5 = pressures5[::2]
        else:
            endelement = pressures5[-1]
            pressures5 = pressures5[::2]
            pressures5.append(endelement)

        for j in range(period - 1):
            pressures6[j] += pressures6[j + 1]
            pressures6[j] /= 2
        if len(pressures6) % 2 == 0:
            pressures6 = pressures6[::2]
        else:
            endelement = pressures6[-1]
            pressures6 = pressures6[::2]
            pressures6.append(endelement)

        for j in range(period - 1):
            pressures7[j] += pressures7[j + 1]
            pressures7[j] /= 2
        if len(pressures7) % 2 == 0:
            pressures7 = pressures7[::2]
        else:
            endelement = pressures7[-1]
            pressures7 = pressures7[::2]
            pressures7.append(endelement)
        k += 1

main()
