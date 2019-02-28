#!/usr/bin/env python
# _*_ coding: utf-8 _*_

def read_config(file):
    parameters = {}
    with open(file) as f:
        records = f.readlines()
    for line in records:
        if line[0] == '#' or line.strip() == '':
            continue
        else:
            # print(line.strip())
            array = line.strip().split('=')
            parameters[array[0]] = array[1]

    for key in parameters:
        if parameters[key].isdigit():
            parameters[key] = int(parameters[key])

    if parameters['Weight_Value'] != '':
        parameters['Weight_Value'] = float(parameters['Weight_Value'])
    if parameters['Cost'] != '':
        parameters['Cost'] = float(parameters['Cost'])

    default = {
        'Sequence_Type': 'DNA',
        'Method': 'binary',
        'Kmer_Size': 3,
        'Sliding_Window': 5,
        'K_Space': 5,
        'Lag_Value': 2,
        'Weight_Value': 0.1,
        'Lamada_Value': 2,
        'Di-DNA-Phychem': 'Twist;Tilt;Roll;Shift;Slide;Rise',
        'Tri-DNA-Phychem': 'Dnase I;Bendability (DNAse)',
        'Di-DNA-Phychem-default6': 'Rise;Roll;Shift;Slide;Tilt;Twist',
        'Di-RNA-Phychem': 'Rise (RNA);Roll (RNA);Shift (RNA);Slide (RNA);Tilt (RNA);Twist (RNA)',
        'AAindex': 'ANDN920101;ARGP820101;ARGP820102;ARGP820103;BEGF750101;BEGF750102;BEGF750103;BHAR880101',
        'PseKRAAC_Model': 'g-gap',
        'Ktuple': 2,
        'GapLamada': 2,
        'RAACCluster1': 2,
        'RAACCluster2': 2,
        'RAACCluster3A': 2,
        'RAACCluster3B': 2,
        'RAACCluster4': 5,
        'RAACCluster5': 2,
        'RAACCluster6A': 5,
        'RAACCluster6B': 5,
        'RAACCluster6C': 5,
        'RAACCluster7': 2,
        'RAACCluster8': 2,
        'RAACCluster9': 2,
        'RAACCluster10': 2,
        'RAACCluster11': 2,
        'RAACCluster12': 2,
        'RAACCluster13': 4,
        'RAACCluster14': 2,
        'RAACCluster15': 2,
        'RAACCluster16': 2,
        'Output_Format': 'svm',
        'Clustering_Algorithm': '',
        'Kmean_Cluster_Number': 2,
        'Clustering_Type': 'sample',
        'Feature_Normalization_Algorithm': '',
        'Feature_Selection_Algorithm': '',
        'Dimension_Reduction_Algorithm': '',
        'Selected_Feature_Number': 3,
        'ML': '',
        'Kernel': 'rbf',
        'Cost': 1.0,
        'Tree_Number': 100,
        'K_Nearest_Neighbour': 3,
        'Hidden_Layer_Size': '32;32',
        'Validation': 5,
    }

    for key in default:
        if key in parameters:
            if parameters[key] == '':
                parameters[key] = default[key]
    return parameters