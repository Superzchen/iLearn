#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse
import sys, os, re
import numpy as np
from pubscripts import *
from clusters import *
from dimreduction import lda as dimlda
from dimreduction import pca as dimpca
from dimreduction import tsne as dimtsne
from featurenormalization import *
from featureselection import *
from machinelearning import *
from itertools import combinations

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="it's usage tip.",
                                     description="running the iLearn pipeline")
    parser.add_argument("--config", default='config.txt', help="the config file")
    args = parser.parse_args()

    # read config
    parameters = read_config.read_config(args.config)
    if parameters['Sequence_Type'] == 'Protein':
        from descproteins import *
    else:
        from descnucleotide import *

    # commands for encoding
    dna_cmd_coding = {
        'Kmer': ['Kmer.Kmer(training_data, k=%s, **kw)' % parameters['Kmer_Size'], 'Kmer.Kmer(testing_data, k=%s, **kw)' % parameters['Kmer_Size']],
        'RCKmer': ['RCKmer.RCKmer(training_data, k=%s, **kw)' % parameters['Kmer_Size'], 'RCKmer.RCKmer(testing_data, k=%s, **kw)' % parameters['Kmer_Size']],
        'NAC': ['NAC.NAC(training_data, **kw)', 'NAC.NAC(testing_data, **kw)'],
        'DNC': ['DNC.DNC(training_data, **kw)', 'DNC.DNC(testing_data, **kw)'],
        'TNC': ['TNC.TNC(training_data, **kw)', 'TNC.TNC(testing_data, **kw)'],
        'ANF': ['ANF.ANF(training_data, **kw)', 'ANF.ANF(testing_data, **kw)'],
        'ENAC': ['ENAC.ENAC(training_data, window=%s, **kw)' % parameters['Sliding_Window'], 'ENAC.ENAC(testing_data, window=%s, **kw)' % parameters['Sliding_Window']],
        'binary': ['binary.binary(training_data, **kw)', 'binary.binary(testing_data, **kw)'],
        'CKSNAP': ['CKSNAP.CKSNAP(training_data, gap=%s, **kw)' % parameters['K_Space'], 'CKSNAP.CKSNAP(testing_data, gap=%s, **kw)' % parameters['K_Space']],
        'NCP': ['NCP.NCP(training_data, **kw)', 'NCP.NCP(testing_data, **kw)'],
        'PSTNPss': ['PSTNPss.PSTNPss(PSTNP_training_data, **kw)', 'PSTNPss.PSTNPss(PSTNP_testing_data, **kw)'],
        'PSTNPds': ['PSTNPds.PSTNPds(PSTNP_training_data, **kw)', 'PSTNPds.PSTNPds(PSTNP_testing_data, **kw)'],
        'EIIP': ['EIIP.EIIP(training_data, **kw)', 'EIIP.EIIP(testing_data, **kw)'],
        'PseEIIP': ['PseEIIP.PseEIIP(training_data, **kw)', 'PseEIIP.PseEIIP(testing_data, **kw)'],
        'DAC': ['ACC.make_ac_vector(training_data, my_property_name, my_property_value, %s, my_kmer)' % parameters['Lag_Value'], 'ACC.make_ac_vector(testing_data, my_property_name, my_property_value, %s, my_kmer)' % parameters['Lag_Value']],
        'DCC': ['ACC.make_cc_vector(training_data, my_property_name, my_property_value, %s, my_kmer)' % parameters['Lag_Value'], 'ACC.make_cc_vector(testing_data, my_property_name, my_property_value, %s, my_kmer)' % parameters['Lag_Value']],
        'DACC': ['ACC.make_acc_vector(training_data, my_property_name, my_property_value, %s, my_kmer)' % parameters['Lag_Value'], 'ACC.make_acc_vector(testing_data, my_property_name, my_property_value, %s, my_kmer)' % parameters['Lag_Value']],
        'TAC': ['ACC.make_ac_vector(training_data, my_property_name, my_property_value, %s, my_kmer)' % parameters['Lag_Value'], 'ACC.make_ac_vector(testing_data, my_property_name, my_property_value, %s, my_kmer)' % parameters['Lag_Value']],
        'TCC': ['ACC.make_cc_vector(training_data, my_property_name, my_property_value, %s, my_kmer)' % parameters['Lag_Value'], 'ACC.make_cc_vector(testing_data, my_property_name, my_property_value, %s, my_kmer)' % parameters['Lag_Value']],
        'TACC': ['ACC.make_acc_vector(training_data, my_property_name, my_property_value, %s, my_kmer)' % parameters['Lag_Value'], 'ACC.make_acc_vector(testing_data, my_property_name, my_property_value, %s, my_kmer)' % parameters['Lag_Value']],
        'PseDNC': ['Pse.make_PseDNC_vector(training_data, my_property_name, my_property_value, my_lamada, my_weight)','Pse.make_PseDNC_vector(testing_data, my_property_name, my_property_value, my_lamada, my_weight)'],
        'PseKNC': ['Pse.make_PseKNC_vector(training_data, my_property_name, my_property_value, my_lamada, my_weight, %d)' % int(parameters['Kmer_Size']), 'Pse.make_PseKNC_vector(testing_data, my_property_name, my_property_value, my_lamada, my_weight, %d)' % int(parameters['Kmer_Size'])],
        'PCPseDNC': ['Pse.make_PseDNC_vector(training_data, my_property_name, my_property_value, my_lamada, my_weight)','Pse.make_PseDNC_vector(testing_data, my_property_name, my_property_value, my_lamada, my_weight)'],
        'PCPseTNC': ['Pse.make_PCPseTNC_vector(training_data, my_property_name, my_property_value, my_lamada, my_weight)','Pse.make_PCPseTNC_vector(testing_data, my_property_name, my_property_value, my_lamada, my_weight)'],
        'SCPseDNC': ['Pse.make_SCPseDNC_vector(training_data, my_property_name, my_property_value, my_lamada, my_weight)','Pse.make_SCPseDNC_vector(testing_data, my_property_name, my_property_value, my_lamada, my_weight)'],
        'SCPseTNC': ['Pse.make_SCPseTNC_vector(training_data, my_property_name, my_property_value, my_lamada, my_weight)', 'Pse.make_SCPseTNC_vector(testing_data, my_property_name, my_property_value, my_lamada, my_weight)'],
    }
    protein_cmd_coding = {
        'AAC': ['AAC.AAC(training_data, **kw)', 'AAC.AAC(testing_data, **kw)'],
        'EAAC': ['EAAC.EAAC(training_data, window=%d, **kw)' % int(parameters['Sliding_Window']), 'EAAC.EAAC(testing_data, window=%d, **kw)' % int(parameters['Sliding_Window'])],
        'CKSAAP': ['CKSAAP.CKSAAP(training_data, gap=%d, **kw)' % int(parameters['K_Space']), 'CKSAAP.CKSAAP(testing_data, gap=%d, **kw)' % int(parameters['K_Space'])],
        'DPC': ['DPC.DPC(training_data, **kw)', 'DPC.DPC(testing_data, **kw)'],
        'DDE': ['DDE.DDE(training_data, **kw)', 'DDE.DDE(testing_data, **kw)'],
        'TPC': ['TPC.TPC(training_data, **kw)', 'TPC.TPC(testing_data, **kw)'],
        'binary': ['binary.binary(training_data, **kw)', 'binary.binary(testing_data, **kw)'],
        'GAAC': ['GAAC.GAAC(training_data, **kw)', 'GAAC.GAAC(testing_data, **kw)'],
        'EGAAC': ['EGAAC.EGAAC(training_data, window=%d, **kw)' % int(parameters['Sliding_Window']), 'EGAAC.EGAAC(testing_data, window=%d, **kw)' % int(parameters['Sliding_Window'])],
        'CKSAAGP': ['CKSAAGP.CKSAAGP(training_data, gap=%d, **kw)' % int(parameters['K_Space']), 'CKSAAGP.CKSAAGP(testing_data, gap=%d, **kw)' % int(parameters['K_Space'])],
        'GDPC': ['GDPC.GDPC(training_data, **kw)', 'GDPC.GDPC(testing_data, **kw)'],
        'GTPC': ['GTPC.GTPC(training_data, **kw)', 'GTPC.GTPC(testing_data, **kw)'],
        'AAINDEX': ['AAINDEX.AAINDEX(training_data, props=props, **kw)', 'AAINDEX.AAINDEX(testing_data, props=props, **kw)'],
        'ZSCALE': ['ZSCALE.ZSCALE(training_data, **kw)', 'ZSCALE.ZSCALE(testing_data, **kw)'],
        'BLOSUM62': ['BLOSUM62.BLOSUM62(training_data, **kw)', 'BLOSUM62.BLOSUM62(testing_data, **kw)'],
        'Moran': ['Moran.Moran(training_data, props=props, nlag=%d, **kw)' % int(parameters['Lag_Value']), 'Moran.Moran(testing_data, props=props, nlag=%d, **kw)' % int(parameters['Lag_Value'])],
        'Geary': ['Geary.Geary(training_data, props=props, nlag=%d, **kw)' % int(parameters['Lag_Value']), 'Geary.Geary(testing_data, props=props, nlag=%d, **kw)' % int(parameters['Lag_Value'])],
        'NMBroto': ['NMBroto.NMBroto(training_data, props=props, nlag=%d, **kw)' % int(parameters['Lag_Value']), 'NMBroto.NMBroto(testing_data, props=props, nlag=%d, **kw)' % int(parameters['Lag_Value'])],
        'CTDC': ['CTDC.CTDC(training_data, **kw)', 'CTDC.CTDC(testing_data, **kw)'],
        'CTDT': ['CTDT.CTDT(training_data, **kw)', 'CTDT.CTDT(testing_data, **kw)'],
        'CTDD': ['CTDD.CTDD(training_data, **kw)', 'CTDD.CTDD(testing_data, **kw)'],
        'CTriad': ['CTriad.CTriad(training_data, gap=0, **kw)', 'CTriad.CTriad(testing_data, gap=0, **kw)'],
        'KSCTriad': ['KSCTriad.KSCTriad(training_data, gap=%d, **kw)' % int(parameters['K_Space']), 'KSCTriad.KSCTriad(testing_data, gap=%d, **kw)' % int(parameters['K_Space'])],
        'SOCNumber': ['SOCNumber.SOCNumber(training_data, nlag=%d, **kw)' % int(parameters['Lag_Value']), 'SOCNumber.SOCNumber(testing_data, nlag=%d, **kw)' % int(parameters['Lag_Value'])],
        'QSOrder': ['QSOrder.QSOrder(training_data, nlag=%d, w=%f, **kw)' % (int(parameters['Lag_Value']), float(parameters['Weight_Value'])), 'QSOrder.QSOrder(testing_data, nlag=%d, w=%f, **kw)' % (int(parameters['Lag_Value']), float(parameters['Weight_Value']))],
        'PAAC': ['PAAC.PAAC(training_data, lambdaValue=%d, w=%f, **kw)' % (int(parameters['Lamada_Value']), float(parameters['Weight_Value'])), 'PAAC.PAAC(testing_data, lambdaValue=%d, w=%f, **kw)' % (int(parameters['Lamada_Value']), float(parameters['Weight_Value']))],
        'APAAC': ['APAAC.APAAC(training_data, lambdaValue=%d, w=%f, **kw)' % (int(parameters['Lamada_Value']), float(parameters['Weight_Value'])), 'APAAC.APAAC(testing_data, lambdaValue=%d, w=%f, **kw)' % (int(parameters['Lamada_Value']), float(parameters['Weight_Value']))],
        'KNNprotein': ['KNNprotein.KNNprotein(PSTNP_training_data, **kw)', 'KNNprotein.KNNprotein(PSTNP_testing_data, **kw)'],
        'KNNpeptide': ['KNNpeptide.KNNpeptide(PSTNP_training_data, **kw)', 'KNNpeptide.KNNpeptide(PSTNP_testing_data, **kw)'],
        'Kmer': ['Kmer.Kmer(training_data, k=%d, type="Protein", **kw)' % int(parameters['Kmer_Size']), 'Kmer.Kmer(testing_data, k=%d, type="Protein", **kw)' % int(parameters['Kmer_Size'])],
        'type1': ['type1.type1(training_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster1']), int(parameters['Ktuple']), int(parameters['GapLamada'])), 'type1.type1(testing_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster1']), int(parameters['Ktuple']), int(parameters['GapLamada']))],
        'type2': ['type2.type1(training_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster2']), int(parameters['Ktuple']), int(parameters['GapLamada'])), 'type2.type1(testing_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster2']), int(parameters['Ktuple']), int(parameters['GapLamada']))],
        'type3A': ['type3A.type1(training_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster3A']), int(parameters['Ktuple']), int(parameters['GapLamada'])), 'type3A.type1(testing_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster3A']), int(parameters['Ktuple']), int(parameters['GapLamada']))],
        'type3B': ['type3B.type1(training_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster3B']), int(parameters['Ktuple']), int(parameters['GapLamada'])), 'type3B.type1(testing_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster3B']), int(parameters['Ktuple']), int(parameters['GapLamada']))],
        'type4': ['type4.type1(training_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster4']), int(parameters['Ktuple']), int(parameters['GapLamada'])), 'type4.type1(testing_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster4']), int(parameters['Ktuple']), int(parameters['GapLamada']))],
        'type5': ['type5.type1(training_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster5']), int(parameters['Ktuple']), int(parameters['GapLamada'])), 'type5.type1(testing_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster5']), int(parameters['Ktuple']), int(parameters['GapLamada']))],
        'type6A': ['type6A.type1(training_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster6A']), int(parameters['Ktuple']), int(parameters['GapLamada'])), 'type6A.type1(testing_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster6A']), int(parameters['Ktuple']), int(parameters['GapLamada']))],
        'type6B': ['type6B.type1(training_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster6B']), int(parameters['Ktuple']), int(parameters['GapLamada'])), 'type6B.type1(testing_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster6B']), int(parameters['Ktuple']), int(parameters['GapLamada']))],
        'type6C': ['type6C.type1(training_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster6C']), int(parameters['Ktuple']), int(parameters['GapLamada'])), 'type6C.type1(testing_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster6C']), int(parameters['Ktuple']), int(parameters['GapLamada']))],
        'type7': ['type7.type1(training_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster7']), int(parameters['Ktuple']), int(parameters['GapLamada'])), 'type7.type1(testing_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster7']), int(parameters['Ktuple']), int(parameters['GapLamada']))],
        'type8': ['type8.type1(training_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster8']), int(parameters['Ktuple']), int(parameters['GapLamada'])), 'type8.type1(testing_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster8']), int(parameters['Ktuple']), int(parameters['GapLamada']))],
        'type9': ['type9.type1(training_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster9']), int(parameters['Ktuple']), int(parameters['GapLamada'])), 'type9.type1(testing_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster9']), int(parameters['Ktuple']), int(parameters['GapLamada']))],
        'type10': ['type10.type1(training_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster10']), int(parameters['Ktuple']), int(parameters['GapLamada'])), 'type10.type1(testing_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster10']), int(parameters['Ktuple']), int(parameters['GapLamada']))],
        'type11': ['type11.type1(training_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster11']), int(parameters['Ktuple']), int(parameters['GapLamada'])), 'type11.type1(testing_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster11']), int(parameters['Ktuple']), int(parameters['GapLamada']))],
        'type12': ['type12.type1(training_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster12']), int(parameters['Ktuple']), int(parameters['GapLamada'])), 'type12.type1(testing_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster12']), int(parameters['Ktuple']), int(parameters['GapLamada']))],
        'type13': ['type13.type1(training_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster13']), int(parameters['Ktuple']), int(parameters['GapLamada'])), 'type13.type1(testing_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster13']), int(parameters['Ktuple']), int(parameters['GapLamada']))],
        'type14': ['type14.type1(training_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster14']), int(parameters['Ktuple']), int(parameters['GapLamada'])), 'type14.type1(testing_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster14']), int(parameters['Ktuple']), int(parameters['GapLamada']))],
        'type15': ['type15.type1(training_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster15']), int(parameters['Ktuple']), int(parameters['GapLamada'])), 'type15.type1(testing_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster15']), int(parameters['Ktuple']), int(parameters['GapLamada']))],
        'type16': ['type16.type1(training_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster16']), int(parameters['Ktuple']), int(parameters['GapLamada'])), 'type16.type1(testing_data, "%s", %d, %d, %d)' %(parameters['PseKRAAC_Model'], int(parameters['RAACCluster16']), int(parameters['Ktuple']), int(parameters['GapLamada']))],
    }

    # Error information
    error_array = []

    # read fasta sequence and specify cmd
    fastas = []
    cmd_coding = {}
    if parameters['Sequence_Type'] in ('DNA', 'RNA'):
        fastas = read_fasta_sequences.read_nucleotide_sequences(parameters['Sequence_File'])
        cmd_coding = dna_cmd_coding
    elif parameters['Sequence_Type'] == 'Protein':
        fastas = read_fasta_sequences.read_protein_sequences(parameters['Sequence_File'])
        cmd_coding = protein_cmd_coding
    else:
        error_array.append('Sequence type can only be selected in "DNA", "RNA" or "Protein".')

    kw = {'nclusters': 3, 'sof': 'sample', 'order': ''}
    kw['order'] = 'ACGT' if parameters['Sequence_Type'] == 'DNA' or parameters['Sequence_Type'] == 'RNA' else 'ACDEFGHIKLMNPQRSTVWY'

    # divide training and testing data
    training_data = []
    testing_data = []
    PSTNP_training_data = []
    PSTNP_testing_data = []
    for sequence in fastas:
        if sequence[3] == 'training':
            training_data.append(sequence)
            PSTNP_training_data.append(sequence)
            PSTNP_training_data.append([sequence[0], sequence[1], sequence[2], 'testing'])
            PSTNP_testing_data.append(sequence)
        else:
            testing_data.append(sequence)
            PSTNP_testing_data.append(sequence)

    # get property for AAindex, NMBroto, Geary, Moran
    props = parameters['AAindex'].split(';') if parameters['AAindex'] != '' else ['CIDH920105', 'BHAR880101',
                                                                                  'CHOC760101', 'BIGC670101',
                                                                                  'CHAM810101', 'DAYM780201']


    # get property for ACC descriptors and Pse descriptors
    my_property_name, my_property_value, my_kmer, my_lamada, my_weight = 0, 0, 0, 0, 0

    # calculate descriptor for training data
    training_code_dict = {}
    testing_code_dict = {}
    method_array = parameters['Method'].split(';')
    for method in method_array:
        # before extract these 12 descriptor, the property need be extracted
        if method in ('DAC', 'DCC', 'DACC', 'TAC', 'TCC', 'TACC'):
            my_property_name, my_property_value, my_kmer = check_parameters.check_acc_arguments_pipeline(parameters,
                                                                                                         method)
        if method in ('PseDNC', 'PseKNC', 'PCPseDNC', 'PCPseTNC', 'SCPseDNC', 'SCPseTNC'):
            my_property_name, my_property_value, my_lamada, my_weight = check_parameters.check_Pse_arguments_pipeline(
                parameters, method, fastas)
        # calculate descriptors
        training_code_dict[method] = eval(cmd_coding[method][0])
        if len(testing_data) > 0:
            testing_code_dict[method] = eval(cmd_coding[method][1])

    training_code = np.array(training_code_dict[method_array[0]])
    testing_code = []
    if len(testing_data) > 0:
        testing_code = np.array(testing_code_dict[method_array[0]])

    for i in range(1, len(method_array)):
        if training_code_dict[method_array[i]] != 0:
            training_code = np.concatenate((training_code, np.array(training_code_dict[method_array[i]])[:, 2:]), axis=1)
            if len(testing_data) > 0:
                if testing_code_dict[method_array[i]] != 0:
                    testing_code = np.concatenate((testing_code, np.array(testing_code_dict[method_array[i]])[:, 2:]), axis=1)

    if len(testing_data) != 0 and training_code.shape[1] != testing_code.shape[1]:
        error_array.append('Descriptor(s) for testing data calculating failed.')
        testing_data = []

    training_code = training_code.tolist()
    save_file.save_file(training_code, format=parameters['Output_Format'], file='training_code.txt')
    save_file.save_file(training_code, format='tsv_1', file='training_code_1.tsv')

    if len(testing_data) > 0:
        testing_code = testing_code.tolist()
        save_file.save_file(testing_code, format=parameters['Output_Format'], file='testing_code.txt')
        save_file.save_file(testing_code, format='tsv_1', file='testing_code_1.tsv')

    # clustering for training data
    kw['sof'] = parameters['Clustering_Type'] if parameters['Clustering_Type'] != '' else 'sample'
    kw['nclusters'] = parameters['Kmean_Cluster_Number'] if parameters['Kmean_Cluster_Number'] else 2
    if parameters['Clustering_Algorithm'] in ('kmeans', 'hcluster', 'apc', 'meanshift', 'dbscan'):
        cluster_method = parameters['Clustering_Algorithm'].strip()
        training_clustering_data = read_code.read_tsv_1('training_code_1.tsv')
        cmd = cluster_method + '.' + cluster_method + '(training_clustering_data, **kw)'
        clusters_res, e = eval(cmd)
        save_file.save_cluster_result(clusters_res, e, cluster_method + '.txt')
        draw_plot.plot_clustering_2d(training_clustering_data, clusters_res, cluster_method, **kw)
        if e != '':
            error_array.append(e)

    # prepare data for feature normalization, selection and dimension reduction
    training_code_file = 'training_code.txt'
    testing_code_file = 'testing_code.txt'
    training_code_used, training_labels = read_code.read_code(training_code_file, format=parameters['Output_Format'])
    testing_code_used, testing_labels = [], []
    if len(testing_data) > 0:
        testing_code_used, testing_labels = read_code.read_code(testing_code_file, format=parameters['Output_Format'])

    # feature normalization
    training_code_normalized = []
    testing_code_normalized = []
    if parameters['Feature_Normalization_Algorithm'] != '' and parameters['Feature_Normalization_Algorithm'] in ('ZScore, MinMax'):
        cmd = parameters['Feature_Normalization_Algorithm'] + '.' + parameters[
            'Feature_Normalization_Algorithm'] + '(training_code_used, training_labels)'
        training_code_normalized, e = eval(cmd)
        save_file.save_file(training_code_normalized, format=parameters['Output_Format'],
                            file='training_code_normalized.txt')
        training_code_file = 'training_code_normalized.txt'

        if len(testing_data) > 0:
            cmd = parameters['Feature_Normalization_Algorithm'] + '.' + parameters[
                'Feature_Normalization_Algorithm'] + '(testing_code_used, testing_labels)'
            testing_code_normalized, e = eval(cmd)
            save_file.save_file(testing_code_normalized, format=parameters['Output_Format'],
                                file='testing_code_normalized.txt')
            testing_code_file = 'testing_code_normalized.txt'

    # feature selection
    training_code_selected = []
    testing_code_selected = []
    if parameters['Feature_Selection_Algorithm'] != '' and parameters['Feature_Selection_Algorithm'] in (
    'CHI2', 'IG', 'MIC', 'pearsonr', 'Fscore'):
        training_code_used, training_labels = read_code.read_code(training_code_file,
                                                                  format=parameters['Output_Format'])
        if len(testing_data) > 0:
            testing_code_used, testing_labels = read_code.read_code(testing_code_file,
                                                                    format=parameters['Output_Format'])

        cmd = parameters['Feature_Selection_Algorithm'] + '.' + parameters[
            'Feature_Selection_Algorithm'] + '(training_code_used, training_labels)'
        selected_features, e = eval(cmd)
        save_file.save_FS_result(selected_features, e, parameters['Feature_Selection_Algorithm'], 'feaure_rank.txt')
        training_code_selected = select_features.select_features(training_code_used, training_labels,
                                                                 selected_features,
                                                                 parameters['Selected_Feature_Number']).tolist()
        save_file.save_file(training_code_selected, parameters['Output_Format'], 'training_code_selected.txt')
        training_code_file = 'training_code_selected.txt'

        if len(testing_data) > 0:
            testing_code_selected = select_features.select_features(testing_code_used, testing_labels,
                                                                    selected_features,
                                                                    parameters['Selected_Feature_Number']).tolist()
            save_file.save_file(testing_code_selected, parameters['Output_Format'], 'testing_code_selected.txt')
            testing_code_file = 'testing_code_selected.txt'

    # dimension reduction
    if parameters['Dimension_Reduction_Algorithm'] != '' and parameters['Dimension_Reduction_Algorithm'] in ('pca', 'lda', 'tsne'):
        training_code_used, training_labels = read_code.read_code(training_code_file,
                                                                  format=parameters['Output_Format'])
        if len(testing_data) > 0:
            testing_code_used, testing_labels = read_code.read_code(testing_code_file,
                                                                    format=parameters['Output_Format'])
        n_components = parameters['Dimension_Reduction_Number'] if parameters['Dimension_Reduction_Number'] != '' else 3
        training_code_reduced = []
        if parameters['Dimension_Reduction_Algorithm'] == 'pca':
            training_code_reduced = dimpca.pca(training_code_used, n_components=n_components)
        if parameters['Dimension_Reduction_Algorithm'] == 'lda':
            training_code_reduced = dimlda.lda(training_code_used, training_labels, n_components=n_components)
        if parameters['Dimension_Reduction_Algorithm'] == 'tsne':
            training_code_reduced = dimtsne.tsne(training_code_used[1:, 1:].astype(float), no_dims=n_components)
            new_data = np.zeros((training_code_reduced.shape[0], training_code_reduced.shape[1] + 1))
            new_data[:, 1:] = training_code_reduced
            new_data = new_data.astype(str)
            new_data[:, 0] = np.array(training_code_used[1:])[:, 0]
            training_code_reduced = new_data.tolist()
        save_file.save_reduction_result(training_code_reduced, file='training_code_dimension_reduction.txt')
        if n_components >= 2:
            draw_plot.plot_2d(training_code_reduced, training_labels, file='training_dimension_reduction_2d.png')
        if n_components >= 3:
            draw_plot.plot_3d(training_code_reduced, training_labels, file='training_dimension_reduction_3d.png')

    # machine learning
    ML_array = parameters['ML'].split(';')
    if parameters['ML'] != '' and parameters['Validation'] != '':
        X, y, independent = 0, 0, np.array([])
        X, y = read_code_ml.read_code(training_code_file, format=parameters['Output_Format'])
        classes = sorted(list(set(y)))
        if len(testing_data) > 0:
            ind_X, ind_y = read_code_ml.read_code(testing_code_file, format=parameters['Output_Format'])
            independent = np.zeros((ind_X.shape[0], ind_X.shape[1] + 1))
            independent[:, 0], independent[:, 1:] = ind_y, ind_X

        fold = parameters['Validation'] if parameters['Validation'] != '' else 5
        # svm
        kernel = parameters['Kernel'] if parameters['Kernel'] != '' else 'rbf'
        auto = False
        cost = 1.0
        gamma = 1 / len(training_data)
        if parameters['Auto_Opterimize'] == 'True':
            auto = True
            cost = 1.0
            gamma = 1 / len(training_data)
        else:
            cost = float(parameters['Cost']) if parameters['Cost'] != '' else 1.0
            gamma = float(parameters['Gamma']) if parameters['Gamma'] != '' else 1 / len(training_data)

        # RF
        n_trees = parameters['Tree_Number'] if parameters['Tree_Number'] != '' else 100
        # KNN
        top_K = int(parameters['K_Nearest_Neighbour']) if parameters['K_Nearest_Neighbour'] != '' else 3
        # ANN
        hidden_layer_size = eval('(' + re.sub(';', ', ', parameters['Hidden_Layer_Size']) + ')') if parameters['Hidden_Layer_Size'] != '' else (32, 32)

        ensemble_train_data = {}  # for ensemble learning
        ensemble_test_data = {}
        # training model
        if len(classes) == 2:
            ML_AUCs_dict = {}
            para_info_dict = {}
            for ML in ML_array:
                para_info, cv_res, ind_res = 0, 0, 0
                if ML == 'SVM':
                    para_info, cv_res, ind_res = SVM.SVM_Classifier(X, y, indep=independent, fold=fold, batch=0.8,
                                                                    auto=auto, kernel=kernel, gamma=gamma, C=cost)
                    con = np.concatenate(tuple((f for f in cv_res)), axis=0)
                    ensemble_train_data[ML] = con
                    ensemble_test_data[ML] = ind_res
                if ML == 'RF':
                    para_info, cv_res, ind_res = RF.RF_Classifier(X, y, indep=independent, fold=fold, n_trees=n_trees)
                    con = np.concatenate(tuple((f for f in cv_res)), axis=0)
                    ensemble_train_data[ML] = con
                    ensemble_test_data[ML] = ind_res
                if ML == 'LR':
                    para_info, cv_res, ind_res = LR.LR_Classifier_binary(X, y, indep=independent, fold=fold)
                    con = np.concatenate(tuple((f for f in cv_res)), axis=0)
                    ensemble_train_data[ML] = con
                    ensemble_test_data[ML] = ind_res
                if ML == 'KNN':
                    para_info, cv_res, ind_res = KNN.KNN_Classifier_binary(X, y, indep=independent, fold=fold,
                                                                           n_neighbors=top_K)
                    con = np.concatenate(tuple((f for f in cv_res)), axis=0)
                    ensemble_train_data[ML] = con
                    ensemble_test_data[ML] = ind_res
                if ML == 'MLP':
                    para_info, cv_res, ind_res = MLP.MLP_Classifier(X, y, indep=independent, fold=fold,
                                                                    hidden_layer_size=hidden_layer_size)
                    con = np.concatenate(tuple((f for f in cv_res)), axis=0)
                    ensemble_train_data[ML] = con
                    ensemble_test_data[ML] = ind_res
                para_info_dict[ML] = para_info

                save_file.save_CV_result_binary(cv_res, '%s_CV.txt' % ML, para_info)
                ML_AUCs_dict[ML] = draw_plot.plot_roc_cv(cv_res, '%s_ROC_CV.png' % ML, label_column=0, score_column=2)
                mean_auprc = draw_plot.plot_prc_CV(cv_res, '%s_PRC_CV.png' % ML, label_column=0, score_column=2)
                cv_metrics = calculate_prediction_metrics.calculate_metrics_cv(cv_res, label_column=0, score_column=2)
                save_file.save_prediction_metrics_cv(cv_metrics, '%s_metrics_CV.txt' % ML)

                if len(testing_data) > 0:
                    save_file.save_IND_result_binary(ind_res, '%s_IND.txt' % ML, para_info)
                    ind_auc = draw_plot.plot_roc_ind(ind_res, '%s_ROC_IND.png' % ML, label_column=0, score_column=2)
                    ind_auprc = draw_plot.plot_prc_ind(ind_res, '%s_PRC_IND.png' % ML, label_column=0, score_column=2)
                    ind_metrics = calculate_prediction_metrics.calculate_metrics(ind_res[:, 0], ind_res[:, 2])
                    save_file.save_prediction_metrics_ind(ind_metrics, '%s_metrics_IND.txt' % ML)

            best_ML = ''
            best_AUC = 0
            for ml in ML_array:
                if best_AUC < ML_AUCs_dict[ml]:
                    best_AUC = ML_AUCs_dict[ml]
                    best_ML = ml
            print('The model with best performance is : %s' %best_ML)

        if len(classes) >= 3:
            acc_dict = {}
            para_info_dict = {}
            for ML in ML_array:
                para_info, cv_res, ind_res = 0, 0, 0
                if ML == 'SVM':
                    para_info, cv_res, ind_res = SVM.SVM_Classifier(X, y, indep=independent, fold=fold, batch=0.8,
                                                                    auto=auto, kernel=kernel, gamma=gamma, C=cost)
                    con = np.concatenate(tuple((f for f in cv_res)), axis=0)
                    ensemble_train_data[ML] = con
                    ensemble_test_data[ML] = ind_res
                if ML == 'RF':
                    para_info, cv_res, ind_res = RF.RF_Classifier(X, y, indep=independent, fold=fold, n_trees=n_trees)
                    con = np.concatenate(tuple((f for f in cv_res)), axis=0)
                    ensemble_train_data[ML] = con
                    ensemble_test_data[ML] = ind_res
                if ML == 'LR':
                    para_info, cv_res, ind_res = LR.LR_Classifier_binary(X, y, indep=independent, fold=fold)
                    con = np.concatenate(tuple((f for f in cv_res)), axis=0)
                    ensemble_train_data[ML] = con
                    ensemble_test_data[ML] = ind_res
                if ML == 'KNN':
                    para_info, cv_res, ind_res = KNN.KNN_Classifier_binary(X, y, indep=independent, fold=fold,
                                                                           n_neighbors=top_K)
                    con = np.concatenate(tuple((f for f in cv_res)), axis=0)
                    ensemble_train_data[ML] = con
                    ensemble_test_data[ML] = ind_res
                if ML == 'MLP':
                    para_info, cv_res, ind_res = MLP.MLP_Classifier(X, y, indep=independent, fold=fold,
                                                                    hidden_layer_size=hidden_layer_size)
                    con = np.concatenate(tuple((f for f in cv_res)), axis=0)
                    ensemble_train_data[ML] = con
                    ensemble_test_data[ML] = ind_res
                para_info_dict[ML] = para_info

                save_file.save_CV_result(cv_res, classes, '%s_CV.txt' % ML, para_info)
                cv_metrics = calculate_prediction_metrics.calculate_metrics_cv_muti(cv_res, classes, label_column=0)
                save_file.save_prediction_metrics_cv_muti(cv_metrics, classes, '%s_metrics_CV.txt' % ML)

                # calculate mean acc
                mean_acc = 0
                for res in cv_metrics:
                    for c in classes:
                        mean_acc += res[c]
                acc_dict[ML] = mean_acc

                if len(testing_data) > 0:
                    save_file.save_IND_result(ind_res, classes, '%s_IND.txt' % ML, para_info)
                    ind_metrics = calculate_prediction_metrics.calculate_metrics_ind_muti(ind_res, classes,
                                                                                          label_column=0)
                    save_file.save_prediction_metrics_ind_muti(ind_metrics, classes, '%s_metrics_IND.txt' % ML)

            best_ML = ''
            best_acc = 0
            for ML in ML_array:
                if best_acc < acc_dict[ML]:
                    best_acc = acc_dict[ML]
                    best_ML = ML
            print('The model with best performance is : %s' %best_ML)

    # ensemble learning
    if len(ML_array) >= 2 and parameters['Ensemble'].upper() == 'YES':
        fold = parameters['Validation'] if parameters['Validation'] != '' else 5
        ML_combinations = []
        for i in range(2, len(ML_array) + 1):
            ML_combinations += list(combinations(ML_array, i))

        if len(classes) == 2:
            Cb_AUCs_dict = {}
            for c in ML_combinations:
                ensemble_data_X, ensemble_data_y = np.concatenate(
                    tuple(ensemble_train_data[m][:, 2].reshape((-1, 1)) for m in c), axis=1), ensemble_train_data[c[0]][
                                                                                              :, 0].astype(int)
                ensemble_data_ind = np.array([])
                if len(testing_data) > 0:
                    ensemble_data_ind = np.concatenate(tuple(ensemble_test_data[m][:, 2].reshape((-1, 1)) for m in c),
                                                       axis=1)
                    ensemble_data_ind = np.concatenate(
                        (ensemble_test_data[c[0]][:, 0].reshape((-1, 1)), ensemble_data_ind), axis=1)
                para_info, cv_res, ind_res = LR.LR_Classifier_binary(ensemble_data_X, ensemble_data_y,
                                                                     indep=ensemble_data_ind, fold=fold)

                save_file.save_CV_result_binary(cv_res, '%s_CV.txt' % '-'.join(c), '')
                Cb_AUCs_dict['-'.join(c)] = draw_plot.plot_roc_cv(cv_res, '%s_ROC_CV.png' % '-'.join(c), label_column=0,
                                                                  score_column=2)
                mean_auprc = draw_plot.plot_prc_CV(cv_res, '%s_PRC_CV.png' % '-'.join(c), label_column=0,
                                                   score_column=2)
                cv_metrics = calculate_prediction_metrics.calculate_metrics_cv(cv_res, label_column=0, score_column=2)
                save_file.save_prediction_metrics_cv(cv_metrics, '%s_metrics_CV.txt' % '-'.join(c))

                if len(testing_data) > 0:
                    save_file.save_IND_result_binary(ind_res, '%s_IND.txt' % '-'.join(c), para_info)
                    ind_auc = draw_plot.plot_roc_ind(ind_res, '%s_ROC_IND.png' % '-'.join(c), label_column=0,
                                                     score_column=2)
                    ind_auprc = draw_plot.plot_prc_ind(ind_res, '%s_PRC_IND.png' % '-'.join(c), label_column=0,
                                                       score_column=2)
                    ind_metrics = calculate_prediction_metrics.calculate_metrics(ind_res[:, 0], ind_res[:, 2])
                    save_file.save_prediction_metrics_ind(ind_metrics, '%s_metrics_IND.txt' % '-'.join(c))

            best_ML_comb = ''
            best_AUC_comb = 0
            for cb in Cb_AUCs_dict:
                if best_AUC_comb < Cb_AUCs_dict[cb]:
                    best_AUC_comb = Cb_AUCs_dict[cb]
                    best_ML_comb = cb
            print('The best model combination is: %s' %cb)

        if len(classes) >= 3:
            Cb_acc_dict = {}
            for c in ML_combinations:
                ensemble_data_X, ensemble_data_y = np.concatenate(tuple(ensemble_train_data[m][:, 1:] for m in c),
                                                                  axis=1), ensemble_train_data[c[0]][:, 0].astype(int)
                ensemble_data_ind = np.array([])
                if len(testing_data) > 0:
                    ensemble_data_ind = np.concatenate(tuple(ensemble_test_data[m][:, 1:] for m in c), axis=1)
                    ensemble_data_ind = np.concatenate(
                        (ensemble_test_data[c[0]][:, 0].reshape((-1, 1)), ensemble_data_ind), axis=1)
                para_info, cv_res, ind_res = LR.LR_Classifier_binary(ensemble_data_X, ensemble_data_y,
                                                                     indep=ensemble_data_ind, fold=fold)
                save_file.save_CV_result(cv_res, classes, '%s_CV.txt' % '-'.join(c), '')
                cv_metrics = calculate_prediction_metrics.calculate_metrics_cv_muti(cv_res, classes, label_column=0)
                save_file.save_prediction_metrics_cv_muti(cv_metrics, classes, '%s_metrics_CV.txt' % '-'.join(c))

                ## calculate mean acc
                mean_acc = 0
                for res in cv_metrics:
                    for cl in classes:
                        mean_acc += res[cl]
                Cb_acc_dict['-'.join(c)] = mean_acc

                if len(testing_data) > 0:
                    save_file.save_IND_result(ind_res, classes, '%s_IND.txt' % '-'.join(c), para_info)
                    ind_metrics = calculate_prediction_metrics.calculate_metrics_ind_muti(ind_res, classes,
                                                                                          label_column=0)
                    save_file.save_prediction_metrics_ind_muti(ind_metrics, classes, '%s_metrics_IND.txt' % '-'.join(c))

            best_ML_comb = ''
            best_acc_comb = 0
            for cb in Cb_acc_dict:
                if best_acc_comb < Cb_acc_dict[cb]:
                    best_acc_comb = Cb_acc_dict[cb]
                    best_ML_comb = cb
            print('The best model combination is: %s' %cb)

