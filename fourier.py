#!/usr/bin/env python3

import math
import pprint
pprint = pprint.PrettyPrinter(indent=4).pprint

fname = '343polywiggle.clean'
#fname = '/media/anton/Seagate Expansion Drive/Article-MMT/Cluster calculations for article/BiggerSystems/Polymer/1.1 - 500K dynamics/'

def appro_fourier(pressure, period):
    meanings = [[]]
    i = a0 = a1 = a2 = b1 = b2 = 0
    suma0 = suma1 = suma2 = sumb1 = sumb2 = magnitude = 0

    for i in range(1, period):
        meanings.append(float(pressure[i][7])) # 7777
        suma0 += float(meanings[i])
        suma1 += float(meanings[i]) * math.cos(2 * math.pi / period * (i + 1))
        sumb1 += float(meanings[i]) * math.sin(2 * math.pi / period * (i + 1))
        suma2 += float(meanings[i]) * math.cos(4 * math.pi / period * (i + 1))
        sumb2 += float(meanings[i]) * math.sin(4 * math.pi / period * (i + 1))

    a0 = suma0 / period
    a1 = 2 * suma1 / period
    a2 = 2 * suma2 / period
    b1 = 2 * sumb1 / period
    b2 = 2 * sumb2 / period
    magnitude = 5 * math.sqrt(a1**2 + b1**2)
    print(magnitude)
    return (a0, a1, b1)
   # return None

def error(pressures, a0, a1, b1):
    error = 0
    period = len(pressures)
    for i in range(period):
        error += abs((a0 +
                      a1 * math.cos(2 * math.pi / period * (i+1)) +
                      b1 * math.sin(2 * math.pi / period * (i+1))) - pressures[i])**2
    print(math.sqrt(error/len(pressures)))
    return None

def appro():
    f = open(fname, 'r')
    pressures = []
    i = a0 = a1 = a2 = b1 = b2 = 0
    suma0 = suma1 = suma2 = sumb1 = sumb2 = magnitude = 0
    for line in f:
        i += 1
        if (i - 1) % 5000 == 0:
            i += 1
            continue
        line_s = line.split()
        if len(line_s) < 1:
            continue
        value  = float(line_s[8])
        pressures.append(value)

    #pressures = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

    for k in range(15):
        period = len(pressures)
        print("Period: ", period)

        suma0 = 0
        suma1 = 0
        sumb1 = 0
        sumb2 = 0

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

        print("Mag: ", magnitude)

        error = 0
        for i in range(1, period):
            error += abs((a0 +
                          a1 * math.cos(2 * math.pi / period * (i+1)) +
                          b1 * math.sin(2 * math.pi / period * (i+1))) - pressures[i])**2
            #print(a0 +
            #      a1 * math.cos(2 * math.pi / period * (i+1)) +
            #      b1 * math.sin(2 * math.pi / period * (i+1)), pressures[i])
        print("Dev: ", math.sqrt(error/period))

        #pprint(pressures)
        for j in range(period - 1):
            #pprint(pressures)
            pressures[j] += pressures[j + 1]
            pressures[j] /= 2
        if len(pressures) % 2 == 0:
            pressures = pressures[::2]
        else:
            endelement = pressures[-1]
            pressures = pressures[::2]
            pressures.append(endelement)
        k += 1

    #pprint(pressures)

    return None
appro()
