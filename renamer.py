#!/usr/bin/env python3
# coding utf-8

import shutil

def rename():
    foldernames = ['1507014 after 1481676 no wiggle/',
                   '1508763 after 1507014 no wiggle/',
                   '1515448 after 1508763 no wiggle/',
                   '1518790 after 1515448 no wiggle/']
    j = 1
    for foname in foldernames:
        for i in range(1, 51):
            old_fname = 'pa500mob/' + foname + "co." + str(i * 50000) + ".data"
            new_fname = 'pa500mob/' + "co." + str(j * 50000) + ".data"
            shutil.copyfile(old_fname, new_fname)
            j += 1

rename()
