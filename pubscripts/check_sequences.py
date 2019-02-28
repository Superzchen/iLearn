#!/usr/bin/env python
#_*_coding:utf-8_*_

import re

def check_fasta_with_equal_length(fastas):
    status = True
    lenList = set()
    for i in fastas:
        lenList.add(len(i[1]))
    if len(lenList) == 1:
        return True
    else:
        return False

def get_min_sequence_length(fastas):
    minLen = 10000
    for i in fastas:
        if minLen > len(i[1]):
            minLen = len(i[1])
    return minLen

def get_min_sequence_length_1(fastas):
    minLen = 10000
    for i in fastas:
        if minLen > len(re.sub('-', '', i[1])):
            minLen = len(re.sub('-', '', i[1]))
    return minLen