#!/usr/bin/env python
#_*_coding:utf-8_*_

import re
from collections import Counter

def NAC(fastas, **kw):
    NA = kw['order'] if kw['order'] != None else 'ACGT'
    encodings = []
    header = ['#', 'label']
    for i in NA:
        header.append(i)
    encodings.append(header)

    for i in fastas:
        name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
        count = Counter(sequence)
        for key in count:
            count[key] = count[key]/len(sequence)
        code = [name, label]
        for na in NA:
            code.append(count[na])
        encodings.append(code)
    return encodings