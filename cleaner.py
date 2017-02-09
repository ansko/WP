#!/usr/bin/env python3
# coding utf-8

f = open('log.lammps', 'r')

for line in f:
    i = 0
    line_s = line.split()
    for word in line_s:
        try:
            float(word)
        except:
            i = 1
    if i == 0:
        print(line)
