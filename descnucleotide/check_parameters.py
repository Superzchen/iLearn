#!/usr/bin/env python
# _*_coding:utf-8_*_

import sys, os, platform

pPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(pPath)
data_path = os.path.abspath(
    os.path.dirname(pPath) + os.path.sep + ".") + r'\data' if platform.system() == 'Windows' else os.path.abspath(
    os.path.dirname(pPath) + os.path.sep + ".") + r'/data'
import pickle

didna_list = ['Base stacking', 'Protein induced deformability', 'B-DNA twist', 'Dinucleotide GC Content', 'A-philicity',
              'Propeller twist', 'Duplex stability:(freeenergy)',
              'Duplex tability(disruptenergy)', 'DNA denaturation', 'Bending stiffness', 'Protein DNA twist',
              'Stabilising energy of Z-DNA', 'Aida_BA_transition', 'Breslauer_dG', 'Breslauer_dH',
              'Breslauer_dS', 'Electron_interaction', 'Hartman_trans_free_energy', 'Helix-Coil_transition',
              'Ivanov_BA_transition', 'Lisser_BZ_transition', 'Polar_interaction', 'SantaLucia_dG',
              'SantaLucia_dH', 'SantaLucia_dS', 'Sarai_flexibility', 'Stability', 'Stacking_energy',
              'Sugimoto_dG', 'Sugimoto_dH', 'Sugimoto_dS', 'Watson-Crick_interaction', 'Twist', 'Tilt', 'Roll',
              'Shift', 'Slide', 'Rise',
              'Clash Strength', 'Roll_roll', 'Twist stiffness', 'Tilt stiffness', 'Shift_rise',
              'Adenine content', 'Direction', 'Twist_shift', 'Enthalpy1', 'Twist_twist', 'Roll_shift',
              'Shift_slide', 'Shift2', 'Tilt3', 'Tilt1', 'Tilt4', 'Tilt2', 'Slide (DNA-protein complex)1',
              'Tilt_shift', 'Twist_tilt', 'Twist (DNA-protein complex)1', 'Tilt_rise', 'Roll_rise',
              'Stacking energy', 'Stacking energy1', 'Stacking energy2', 'Stacking energy3', 'Propeller Twist',
              'Roll11', 'Rise (DNA-protein complex)', 'Tilt_tilt', 'Roll4', 'Roll2', 'Roll3', 'Roll1',
              'Minor Groove Size', 'GC content', 'Slide_slide', 'Enthalpy', 'Shift_shift', 'Slide stiffness',
              'Melting Temperature1', 'Flexibility_slide', 'Minor Groove Distance',
              'Rise (DNA-protein complex)1', 'Tilt (DNA-protein complex)', 'Guanine content',
              'Roll (DNA-protein complex)1', 'Entropy', 'Cytosine content', 'Major Groove Size', 'Twist_rise',
              'Major Groove Distance', 'Twist (DNA-protein complex)', 'Purine (AG) content',
              'Melting Temperature', 'Free energy', 'Tilt_slide', 'Major Groove Width', 'Major Groove Depth',
              'Wedge', 'Free energy8', 'Free energy6', 'Free energy7', 'Free energy4', 'Free energy5',
              'Free energy2', 'Free energy3', 'Free energy1', 'Twist_roll', 'Shift (DNA-protein complex)',
              'Rise_rise', 'Flexibility_shift', 'Shift (DNA-protein complex)1', 'Thymine content', 'Slide_rise',
              'Tilt_roll', 'Tip', 'Keto (GT) content', 'Roll stiffness', 'Minor Groove Width', 'Inclination',
              'Entropy1', 'Roll_slide', 'Slide (DNA-protein complex)', 'Twist1', 'Twist3', 'Twist2', 'Twist5',
              'Twist4', 'Twist7', 'Twist6', 'Tilt (DNA-protein complex)1', 'Twist_slide', 'Minor Groove Depth',
              'Roll (DNA-protein complex)', 'Rise2', 'Persistance Length', 'Rise3', 'Shift stiffness',
              'Probability contacting nucleosome core', 'Mobility to bend towards major groove', 'Slide3',
              'Slide2', 'Slide1', 'Shift1', 'Bend', 'Rise1', 'Rise stiffness',
              'Mobility to bend towards minor groove']

tridna_list = ['Dnase I', 'Bendability (DNAse)', 'Bendability (consensus)', 'Trinucleotide GC Content',
               'Nucleosome positioning', 'Consensus_roll', 'Consensus-Rigid', 'Dnase I-Rigid', 'MW-Daltons',
               'MW-kg', 'Nucleosome', 'Nucleosome-Rigid']

dirna_list = ['Slide (RNA)', 'Adenine content', 'Hydrophilicity (RNA)', 'Tilt (RNA)', 'Stacking energy (RNA)',
              'Twist (RNA)', 'Entropy (RNA)', 'Roll (RNA)', 'Purine (AG) content', 'Hydrophilicity (RNA)1',
              'Enthalpy (RNA)1', 'GC content', 'Entropy (RNA)1', 'Rise (RNA)', 'Free energy (RNA)',
              'Keto (GT) content', 'Free energy (RNA)1', 'Enthalpy (RNA)', 'Guanine content', 'Shift (RNA)',
              'Cytosine content', 'Thymine content']

myDict = {
    'DAC': {'DNA': didna_list, 'RNA': dirna_list},
    'DCC': {'DNA': didna_list, 'RNA': dirna_list},
    'DACC': {'DNA': didna_list, 'RNA': dirna_list},
    'TAC': {'DNA': tridna_list, 'RNA': []},
    'TCC': {'DNA': tridna_list, 'RNA': []},
    'TACC': {'DNA': tridna_list, 'RNA': []},
    'PseDNC': {'DNA': didna_list, 'RNA': dirna_list},
    'PseKNC': {'DNA': didna_list, 'RNA': dirna_list},
    'PCPseDNC': {'DNA': didna_list, 'RNA': dirna_list},
    'PCPseTNC': {'DNA': tridna_list, 'RNA': []},
    'SCPseDNC': {'DNA': didna_list, 'RNA': dirna_list},
    'SCPseTNC': {'DNA': tridna_list, 'RNA': []},
}

myDictDefault = {
    'DAC': {'DNA': ['Rise', 'Roll', 'Shift', 'Slide', 'Tilt', 'Twist'],
            'RNA': ['Rise (RNA)', 'Roll (RNA)', 'Shift (RNA)', 'Slide (RNA)', 'Tilt (RNA)', 'Twist (RNA)']},
    'DCC': {'DNA': ['Rise', 'Roll', 'Shift', 'Slide', 'Tilt', 'Twist'],
            'RNA': ['Rise (RNA)', 'Roll (RNA)', 'Shift (RNA)', 'Slide (RNA)', 'Tilt (RNA)', 'Twist (RNA)']},
    'DACC': {'DNA': ['Rise', 'Roll', 'Shift', 'Slide', 'Tilt', 'Twist'],
             'RNA': ['Rise (RNA)', 'Roll (RNA)', 'Shift (RNA)', 'Slide (RNA)', 'Tilt (RNA)', 'Twist (RNA)']},
    'TAC': {'DNA': ['Dnase I', 'Bendability (DNAse)'], 'RNA': []},
    'TCC': {'DNA': ['Dnase I', 'Bendability (DNAse)'], 'RNA': []},
    'TACC': {'DNA': ['Dnase I', 'Bendability (DNAse)'], 'RNA': []},
    'PseDNC': {'DNA': ['Rise', 'Roll', 'Shift', 'Slide', 'Tilt', 'Twist'],
               'RNA': ['Rise (RNA)', 'Roll (RNA)', 'Shift (RNA)', 'Slide (RNA)', 'Tilt (RNA)', 'Twist (RNA)']},
    'PseKNC': {'DNA': ['Rise', 'Roll', 'Shift', 'Slide', 'Tilt', 'Twist'],
               'RNA': ['Rise (RNA)', 'Roll (RNA)', 'Shift (RNA)', 'Slide (RNA)', 'Tilt (RNA)', 'Twist (RNA)']},
    'PCPseDNC': {
        'DNA': ['Base stacking', 'Protein induced deformability', 'B-DNA twist', 'A-philicity', 'Propeller twist',
                'Duplex stability:(freeenergy)', 'DNA denaturation', 'Bending stiffness', 'Protein DNA twist',
                'Aida_BA_transition', 'Breslauer_dG', 'Breslauer_dH', 'Electron_interaction',
                'Hartman_trans_free_energy', 'Helix-Coil_transition', 'Lisser_BZ_transition', 'Polar_interaction',
                'SantaLucia_dG', 'SantaLucia_dS', 'Sarai_flexibility', 'Stability', 'Sugimoto_dG', 'Sugimoto_dH',
                'Sugimoto_dS', 'Duplex tability(disruptenergy)', 'Stabilising energy of Z-DNA', 'Breslauer_dS',
                'Ivanov_BA_transition', 'SantaLucia_dH', 'Stacking_energy', 'Watson-Crick_interaction',
                'Dinucleotide GC Content', 'Twist', 'Tilt', 'Roll', 'Shift', 'Slide', 'Rise'],
        'RNA': ['Rise (RNA)', 'Roll (RNA)', 'Shift (RNA)', 'Slide (RNA)', 'Tilt (RNA)', 'Twist (RNA)']},
    'PCPseTNC': {'DNA': ['Dnase I', 'Bendability (DNAse)'], 'RNA': []},
    'SCPseDNC': {'DNA': ['Rise', 'Roll', 'Shift', 'Slide', 'Tilt', 'Twist'],
                 'RNA': ['Rise (RNA)', 'Roll (RNA)', 'Shift (RNA)', 'Slide (RNA)', 'Tilt (RNA)', 'Twist (RNA)']},
    'SCPseTNC': {'DNA': ['Dnase I', 'Bendability (DNAse)'], 'RNA': []},
}

myKmer = {
    'DAC': 2, 'DCC': 2, 'DACC': 2,
    'TAC': 3, 'TCC': 3, 'TACC': 3
}

myDataFile = {
    'DAC': {'DNA': 'didnaPhyche.data', 'RNA': 'dirnaPhyche.data'},
    'DCC': {'DNA': 'didnaPhyche.data', 'RNA': 'dirnaPhyche.data'},
    'DACC': {'DNA': 'didnaPhyche.data', 'RNA': 'dirnaPhyche.data'},
    'TAC': {'DNA': 'tridnaPhyche.data', 'RNA': ''},
    'TCC': {'DNA': 'tridnaPhyche.data', 'RNA': ''},
    'TACC': {'DNA': 'tridnaPhyche.data', 'RNA': ''},
    'PseDNC': {'DNA': 'didnaPhyche.data', 'RNA': 'dirnaPhyche.data'},
    'PseKNC': {'DNA': 'didnaPhyche.data', 'RNA': 'dirnaPhyche.data'},
    'PCPseDNC': {'DNA': 'didnaPhyche.data', 'RNA': 'dirnaPhyche.data'},
    'PCPseTNC': {'DNA': 'tridnaPhyche.data', 'RNA': ''},
    'SCPseDNC': {'DNA': 'didnaPhyche.data', 'RNA': 'dirnaPhyche.data'},
    'SCPseTNC': {'DNA': 'tridnaPhyche.data', 'RNA': ''},
}


def check_acc_arguments(args):
    if not os.path.exists(args.file):
        print('Error: the input file does not exist.')
        sys.exit(1)
    if args.lag <= 0:
        print('Error: the lag value should be larger than 0.')
        sys.exit(1)

    myNum = 0
    if args.index:
        myNum = myNum + 1
    if args.udi:
        myNum = myNum + 1
    if args.all_index:
        myNum = myNum + 1

    if myNum > 1:
        print(
            'Error: argument is incorrect, "--index", "--udi" and "--all_index" can not be assigned at the same time.')
        sys.exit(1)

    kmer = myKmer[args.method]
    myIndex = []
    myProperty = {}
    dataFile = ''

    if myNum == 0:
        myIndex = myDictDefault[args.method][args.type]
        dataFile = myDataFile[args.method][args.type]
    else:
        if args.index:
            with open(args.index) as f:
                records = f.read().strip().split('\n')
                for i in records:
                    if i in myDict[args.method][args.type]:
                        myIndex.append(i)
                    else:
                        print('Error: there is no "%s" in the index list.' % i)
                        sys.exit(1)
        if args.all_index:
            myIndex = myDict[args.method][args.type]

        if args.udi:
            if not os.path.exists(args.udi):
                print('Error: The user-defined indices file does not exist.')
                sys.exit(1)
            else:
                with open(args.udi) as f:
                    record = f.read().strip().split('\n')[1:]
                for line in record:
                    array = line.strip().split()
                    myProperty[array[0]] = array[1:]
                    myIndex.append(array[0])
        else:
            dataFile = myDataFile[args.method][args.type]

    if dataFile != '':
        with open(data_path + '/' + dataFile, 'rb') as f:
            myProperty = pickle.load(f)

    if len(myIndex) == 0 or len(myProperty) == 0:
        print('Error: arguments is incorrect.')
        sys.exit(1)
    return myIndex, myProperty, kmer


def check_Pse_arguments(args, fastas):
    if not os.path.exists(args.file):
        print('Error: the input file does not exist.')
        sys.exit(1)
    if not 0 < args.weight < 1:
        print('Error: the weight factor ranged from 0 ~ 1.')
        sys.exit(1)
    if not 0 < args.kmer < 10:
        print('Error: the kmer value ranged from 1 - 10')
        sys.exit(1)

    fastaMinLength = 100000000
    for i in fastas:
        if len(i[1]) < fastaMinLength:
            fastaMinLength = len(i[1])
    if not 0 <= args.lamadaValue <= (fastaMinLength - 2):
        print('Error: lamada value error, please see the manual for details.')
        sys.exit(1)

    myNum = 0
    if args.index:
        myNum = myNum + 1
    if args.udi:
        myNum = myNum + 1
    if args.all_index:
        myNum = myNum + 1

    if myNum > 1:
        print(
            'Error: argument is incorrect, "--index", "--udi" and "--all_index" can not be assigned at the same time.')
        sys.exit(1)

    myIndex = []
    myProperty = {}
    dataFile = ''

    if myNum == 0:
        myIndex = myDictDefault[args.method][args.type]
        dataFile = myDataFile[args.method][args.type]
    else:
        if args.index:
            with open(args.index) as f:
                records = f.read().strip().split('\n')
                for i in records:
                    if i in myDict[args.method][args.type]:
                        myIndex.append(i)
                    else:
                        print('Error: there is no "%s" in the index list.' % i)
                        sys.exit(1)
        if args.all_index:
            myIndex = myDict[args.method][args.type]

        if args.udi:
            if not os.path.exists(args.udi):
                print('Error: The user-defined indices file does not exist.')
                sys.exit(1)
            else:
                with open(args.udi) as f:
                    record = f.read().strip().split('\n')[1:]
                for line in record:
                    array = line.strip().split()
                    myProperty[array[0]] = array[1:]
                    myIndex.append(array[0])
        else:
            dataFile = myDataFile[args.method][args.type]

    if dataFile != '':
        with open(data_path + '/' + dataFile, 'rb') as f:
            myProperty = pickle.load(f)

    if len(myIndex) == 0 or len(myProperty) == 0:
        print('Error: arguments is incorrect.')
        sys.exit(1)

    return myIndex, myProperty, args.lamadaValue, args.weight, args.kmer

def check_acc_arguments_pipeline(parameters, method):
    if int(parameters['Lag_Value']) <= 0:
        print('Error: the lag value should be larger than 0.')
        sys.exit(1)

    type = parameters['Sequence_Type']
    kmer = myKmer[method]
    myIndex = []
    myProperty = {}
    dataFile = ''

    if parameters['All_Property'] == 'True':
        myIndex = myDict[method][type]
    else:
        if type == 'DNA':
            if method in ('DAC', 'DCC', 'DACC'):
                myIndex = parameters['Di-DNA-Phychem'].strip().split(';')
            else:
                myIndex = parameters['Tri-DNA-Phychem'].strip().split(';')
        if type == 'RNA':
            if method in ('DAC', 'DCC', 'DACC'):
                myIndex = parameters['Di-RNA-Phychem'].strip().split(';')

    dataFile = myDataFile[method][type]
    if dataFile != '':
        with open(data_path + '/' + dataFile, 'rb') as f:
            myProperty = pickle.load(f)

    if len(myIndex) == 0 or len(myProperty) == 0:
        print('Error: arguments is incorrect.')
        sys.exit(1)
    return myIndex, myProperty, kmer


def check_Pse_arguments_pipeline(parameters, method, fastas):
    if not 0 < float(parameters['Weight_Value']) < 1:
        print('Error: the weight factor ranged from 0 ~ 1.')
        sys.exit(1)
    if not 0 < int(parameters['Kmer_Size']) < 10:
        print('Error: the kmer value ranged from 1 - 10')
        sys.exit(1)

    fastaMinLength = 100000000
    for i in fastas:
        if len(i[1]) < fastaMinLength:
            fastaMinLength = len(i[1])
    if not 0 <= int(parameters['Lamada_Value']) <= (fastaMinLength - 2):
        print('Error: lamada value error, please see the manual for details.')
        sys.exit(1)

    myIndex = []
    myProperty = {}
    dataFile = ''

    type = parameters['Sequence_Type']
    if parameters['All_Property'] == 'True':
        myIndex = myDict[method][type]
    else:
        if type == 'DNA':
            if method in ('PseDNC', 'PseKNC'):
                myIndex = parameters['Di-DNA-Phychem-default6'].strip().split(';')
            if method in ('PCPseDNC', 'SCPseDNC'):
                myIndex = parameters['Di-DNA-Phychem'].strip().split(';')
            if method in ('PCPseTNC', 'SCPseTNC'):
                myIndex = parameters['Tri-DNA-Phychem'].strip().split(';')
        if type == 'RNA':
            if method in ('PseDNC', 'PseKNC'):
                myIndex = parameters['Di-RNA-Phychem'].strip().split(';')
            if method in ('PCPseDNC', 'SCPseDNC'):
                myIndex = parameters['Di-RNA-Phychem'].strip().split(';')

    dataFile = myDataFile[method][type]
    if dataFile != '':
        with open(data_path + '/' + dataFile, 'rb') as f:
            myProperty = pickle.load(f)

    if len(myIndex) == 0 or len(myProperty) == 0:
        print('Error: arguments is incorrect.')
        sys.exit(1)

    return myIndex, myProperty, int(parameters['Lamada_Value']), float(parameters['Weight_Value'])
