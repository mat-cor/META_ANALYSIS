#!/usr/bin/python

import sys, re

POP = ['EUR', 'AFR', 'FIN', 'ARAB', 'SAS', 'EAS', 'OTH', 'CSA', 'AMR', 'NON-EUR', 'HIS', 'ALL']
ANA = ['A1', 'A2', 'B1', 'B2', 'B3', 'C1', 'C2', 'D1']
SEX = ['ALL', 'MALE', 'FEMALE', 'M', 'F']
AGE = ['ALL', 'GT_60', 'LT_60', 'LE_60']

with open(sys.argv[1], 'rt') as f:
    for line in f:
        ok = True
        line = line.strip()
        filename = line.split('/').pop()
        filename = re.sub('.txt|.gz|.bgz', '', filename)
        fields = filename.split('.')
        if (len(fields) != 11):
            print(str(len(fields)) + ' fields\t' + line)
            continue
        try:
            freeze = int(fields[3])
            n_cases = int(fields[7])
            n_controls = int(fields[8])
            date = int(fields[10])
        except ValueError:
            print('Number error\t{}'.format(line))
            ok = False
        if fields[6] not in POP:
            print('Unexpected pop {}\t{}'.format(fields[6], line))
            ok = False
        ana_ok = False
        for ana in ANA:
            if ana in fields[2]:
                ana_ok = True
                break
        if not ana_ok:
            print('Unexpected ANA {}\t{}'.format(fields[2], line))
            ok = False
        if fields[4].upper() not in AGE:
            print('Unexpected age {}\t{}'.format(fields[4], line))
            ok = False
        if fields[5].upper() not in SEX:
            print('Unexpected sex {}\t{}'.format(fields[5], line))
            ok = False
        if ok:
            print('ALLGOOD\t{}'.format(line))
