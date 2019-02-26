#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse
import sys, os, re
import numpy as np
from pubscripts import *
from machinelearning import *


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
        'Kmer': ['Kmer.Kmer(training_data, k=%s, **kw)' % parameters['Kmer_Size'],
                 'Kmer.Kmer(testing_data, k=%s, **kw)' % parameters['Kmer_Size']],
        'RCKmer': ['RCKmer.RCKmer(training_data, k=%s, **kw)' % parameters['Kmer_Size'],
                   'RCKmer.RCKmer(testing_data, k=%s, **kw)' % parameters['Kmer_Size']],
        'NAC': ['NAC.NAC(training_data, **kw)', 'NAC.NAC(testing_data, **kw)'],
        'DNC': ['DNC.DNC(training_data, **kw)', 'DNC.DNC(testing_data, **kw)'],
        'TNC': ['TNC.TNC(training_data, **kw)', 'TNC.TNC(testing_data, **kw)'],
        'ANF': ['ANF.ANF(training_data, **kw)', 'ANF.ANF(testing_data, **kw)'],
        'ENAC': ['ENAC.ENAC(training_data, window=%s, **kw)' % parameters['Sliding_Window'],
                 'ENAC.ENAC(testing_data, window=%s, **kw)' % parameters['Sliding_Window']],
        'binary': ['binary.binary(training_data, **kw)', 'binary.binary(testing_data, **kw)'],
        'CKSNAP': ['CKSNAP.CKSNAP(training_data, gap=%s, **kw)' % parameters['K_Space'],
                   'CKSNAP.CKSNAP(testing_data, gap=%s, **kw)' % parameters['K_Space']],
        'NCP': ['NCP.NCP(training_data, **kw)', 'NCP.NCP(testing_data, **kw)'],
        'PSTNPss': ['PSTNPss.PSTNPss(PSTNP_training_data, **kw)', 'PSTNPss.PSTNPss(PSTNP_testing_data, **kw)'],
        'PSTNPds': ['PSTNPds.PSTNPds(PSTNP_training_data, **kw)', 'PSTNPds.PSTNPds(PSTNP_testing_data, **kw)'],
        'EIIP': ['EIIP.EIIP(training_data, **kw)', 'EIIP.EIIP(testing_data, **kw)'],
        'PseEIIP': ['PseEIIP.PseEIIP(training_data, **kw)', 'PseEIIP.PseEIIP(testing_data, **kw)'],
        'DAC': ['ACC.make_ac_vector(training_data, my_property_name, my_property_value, %s, my_kmer)' % parameters[
            'Lag_Value'],
                'ACC.make_ac_vector(testing_data, my_property_name, my_property_value, %s, my_kmer)' % parameters[
                    'Lag_Value']],
        'DCC': ['ACC.make_cc_vector(training_data, my_property_name, my_property_value, %s, my_kmer)' % parameters[
            'Lag_Value'],
                'ACC.make_cc_vector(testing_data, my_property_name, my_property_value, %s, my_kmer)' % parameters[
                    'Lag_Value']],
        'DACC': [
            'ACC.make_acc_vector(training_data, my_property_name, my_property_value, %s, my_kmer)' % parameters[
                'Lag_Value'],
            'ACC.make_acc_vector(testing_data, my_property_name, my_property_value, %s, my_kmer)' % parameters[
                'Lag_Value']],
        'TAC': ['ACC.make_ac_vector(training_data, my_property_name, my_property_value, %s, my_kmer)' % parameters[
            'Lag_Value'],
                'ACC.make_ac_vector(testing_data, my_property_name, my_property_value, %s, my_kmer)' % parameters[
                    'Lag_Value']],
        'TCC': ['ACC.make_cc_vector(training_data, my_property_name, my_property_value, %s, my_kmer)' % parameters[
            'Lag_Value'],
                'ACC.make_cc_vector(testing_data, my_property_name, my_property_value, %s, my_kmer)' % parameters[
                    'Lag_Value']],
        'TACC': [
            'ACC.make_acc_vector(training_data, my_property_name, my_property_value, %s, my_kmer)' % parameters[
                'Lag_Value'],
            'ACC.make_acc_vector(testing_data, my_property_name, my_property_value, %s, my_kmer)' % parameters[
                'Lag_Value']],
        'PseDNC': [
            'Pse.make_PseDNC_vector(training_data, my_property_name, my_property_value, my_lamada, my_weight)',
            'Pse.make_PseDNC_vector(testing_data, my_property_name, my_property_value, my_lamada, my_weight)'],
        'PseKNC': [
            'Pse.make_PseKNC_vector(training_data, my_property_name, my_property_value, my_lamada, my_weight, %d)' % int(
                parameters['Kmer_Size']),
            'Pse.make_PseKNC_vector(testing_data, my_property_name, my_property_value, my_lamada, my_weight, %d)' % int(
                parameters['Kmer_Size'])],
        'PCPseDNC': [
            'Pse.make_PseDNC_vector(training_data, my_property_name, my_property_value, my_lamada, my_weight)',
            'Pse.make_PseDNC_vector(testing_data, my_property_name, my_property_value, my_lamada, my_weight)'],
        'PCPseTNC': [
            'Pse.make_PCPseTNC_vector(training_data, my_property_name, my_property_value, my_lamada, my_weight)',
            'Pse.make_PCPseTNC_vector(testing_data, my_property_name, my_property_value, my_lamada, my_weight)'],
        'SCPseDNC': [
            'Pse.make_SCPseDNC_vector(training_data, my_property_name, my_property_value, my_lamada, my_weight)',
            'Pse.make_SCPseDNC_vector(testing_data, my_property_name, my_property_value, my_lamada, my_weight)'],
        'SCPseTNC': [
            'Pse.make_SCPseTNC_vector(training_data, my_property_name, my_property_value, my_lamada, my_weight)',
            'Pse.make_SCPseTNC_vector(testing_data, my_property_name, my_property_value, my_lamada, my_weight)'],
    }
    protein_cmd_coding = {
        'AAC': ['AAC.AAC(training_data, **kw)', 'AAC.AAC(testing_data, **kw)'],
        'EAAC': ['EAAC.EAAC(training_data, window=%d, **kw)' % int(parameters['Sliding_Window']),
                 'EAAC.EAAC(testing_data, window=%d, **kw)' % int(parameters['Sliding_Window'])],
        'CKSAAP': ['CKSAAP.CKSAAP(training_data, gap=%d, **kw)' % int(parameters['K_Space']),
                   'CKSAAP.CKSAAP(testing_data, gap=%d, **kw)' % int(parameters['K_Space'])],
        'DPC': ['DPC.DPC(training_data, **kw)', 'DPC.DPC(testing_data, **kw)'],
        'DDE': ['DDE.DDE(training_data, **kw)', 'DDE.DDE(testing_data, **kw)'],
        'TPC': ['TPC.TPC(training_data, **kw)', 'TPC.TPC(testing_data, **kw)'],
        'binary': ['binary.binary(training_data, **kw)', 'binary.binary(testing_data, **kw)'],
        'GAAC': ['GAAC.GAAC(training_data, **kw)', 'GAAC.GAAC(testing_data, **kw)'],
        'EGAAC': ['EGAAC.EGAAC(training_data, window=%d, **kw)' % int(parameters['Sliding_Window']),
                  'EGAAC.EGAAC(testing_data, window=%d, **kw)' % int(parameters['Sliding_Window'])],
        'CKSAAGP': ['CKSAAGP.CKSAAGP(training_data, gap=%d, **kw)' % int(parameters['K_Space']),
                    'CKSAAGP.CKSAAGP(testing_data, gap=%d, **kw)' % int(parameters['K_Space'])],
        'GDPC': ['GDPC.GDPC(training_data, **kw)', 'GDPC.GDPC(testing_data, **kw)'],
        'GTPC': ['GTPC.GTPC(training_data, **kw)', 'GTPC.GTPC(testing_data, **kw)'],
        'AAINDEX': ['AAINDEX.AAINDEX(training_data, props=props, **kw)',
                    'AAINDEX.AAINDEX(testing_data, props=props, **kw)'],
        'ZSCALE': ['ZSCALE.ZSCALE(training_data, **kw)', 'ZSCALE.ZSCALE(testing_data, **kw)'],
        'BLOSUM62': ['BLOSUM62.BLOSUM62(training_data, **kw)', 'BLOSUM62.BLOSUM62(testing_data, **kw)'],
        'Moran': ['Moran.Moran(training_data, props=props, nlag=%d, **kw)' % int(parameters['Lag_Value']),
                  'Moran.Moran(testing_data, props=props, nlag=%d, **kw)' % int(parameters['Lag_Value'])],
        'Geary': ['Geary.Geary(training_data, props=props, nlag=%d, **kw)' % int(parameters['Lag_Value']),
                  'Geary.Geary(testing_data, props=props, nlag=%d, **kw)' % int(parameters['Lag_Value'])],
        'NMBroto': ['NMBroto.NMBroto(training_data, props=props, nlag=%d, **kw)' % int(parameters['Lag_Value']),
                    'NMBroto.NMBroto(testing_data, props=props, nlag=%d, **kw)' % int(parameters['Lag_Value'])],
        'CTDC': ['CTDC.CTDC(training_data, **kw)', 'CTDC.CTDC(testing_data, **kw)'],
        'CTDT': ['CTDT.CTDT(training_data, **kw)', 'CTDT.CTDT(testing_data, **kw)'],
        'CTDD': ['CTDD.CTDD(training_data, **kw)', 'CTDD.CTDD(testing_data, **kw)'],
        'CTriad': ['CTriad.CTriad(training_data, gap=0, **kw)', 'CTriad.CTriad(testing_data, gap=0, **kw)'],
        'KSCTriad': ['KSCTriad.KSCTriad(training_data, gap=%d, **kw)' % int(parameters['K_Space']),
                     'KSCTriad.KSCTriad(testing_data, gap=%d, **kw)' % int(parameters['K_Space'])],
        'SOCNumber': ['SOCNumber.SOCNumber(training_data, nlag=%d, **kw)' % int(parameters['Lag_Value']),
                      'SOCNumber.SOCNumber(testing_data, nlag=%d, **kw)' % int(parameters['Lag_Value'])],
        'QSOrder': ['QSOrder.QSOrder(training_data, nlag=%d, w=%f, **kw)' % (
            int(parameters['Lag_Value']), float(parameters['Weight_Value'])),
                    'QSOrder.QSOrder(testing_data, nlag=%d, w=%f, **kw)' % (
                        int(parameters['Lag_Value']), float(parameters['Weight_Value']))],
        'PAAC': ['PAAC.PAAC(training_data, lambdaValue=%d, w=%f, **kw)' % (
            int(parameters['Lamada_Value']), float(parameters['Weight_Value'])),
                 'PAAC.PAAC(testing_data, lambdaValue=%d, w=%f, **kw)' % (
                     int(parameters['Lamada_Value']), float(parameters['Weight_Value']))],
        'APAAC': ['APAAC.APAAC(training_data, lambdaValue=%d, w=%f, **kw)' % (
            int(parameters['Lamada_Value']), float(parameters['Weight_Value'])),
                  'APAAC.APAAC(testing_data, lambdaValue=%d, w=%f, **kw)' % (
                      int(parameters['Lamada_Value']), float(parameters['Weight_Value']))],
        'KNNprotein': ['KNNprotein.KNNprotein(PSTNP_training_data, **kw)', 'KNNprotein.KNNprotein(PSTNP_testing_data, **kw)'],
        'KNNpeptide': ['KNNpeptide.KNNpeptide(PSTNP_training_data, **kw)', 'KNNpeptide.KNNpeptide(PSTNP_testing_data, **kw)'],
        'Kmer': ['Kmer.Kmer(training_data, k=%d, type="Protein", **kw)' % int(parameters['Kmer_Size']),
                 'Kmer.Kmer(testing_data, k=%d, type="Protein", **kw)' % int(parameters['Kmer_Size'])],
        'type1': ['type1.type1(training_data, "%s", %d, %d, %d)' % (
            parameters['PseKRAAC_Model'], int(parameters['RAACCluster1']), int(parameters['Ktuple']),
            int(parameters['GapLamada'])), 'type1.type1(testing_data, "%s", %d, %d, %d)' % (
                      parameters['PseKRAAC_Model'], int(parameters['RAACCluster1']), int(parameters['Ktuple']),
                      int(parameters['GapLamada']))],
        'type2': ['type2.type1(training_data, "%s", %d, %d, %d)' % (
            parameters['PseKRAAC_Model'], int(parameters['RAACCluster2']), int(parameters['Ktuple']),
            int(parameters['GapLamada'])), 'type2.type1(testing_data, "%s", %d, %d, %d)' % (
                      parameters['PseKRAAC_Model'], int(parameters['RAACCluster2']), int(parameters['Ktuple']),
                      int(parameters['GapLamada']))],
        'type3A': ['type3A.type1(training_data, "%s", %d, %d, %d)' % (
            parameters['PseKRAAC_Model'], int(parameters['RAACCluster3A']), int(parameters['Ktuple']),
            int(parameters['GapLamada'])), 'type3A.type1(testing_data, "%s", %d, %d, %d)' % (
                       parameters['PseKRAAC_Model'], int(parameters['RAACCluster3A']), int(parameters['Ktuple']),
                       int(parameters['GapLamada']))],
        'type3B': ['type3B.type1(training_data, "%s", %d, %d, %d)' % (
            parameters['PseKRAAC_Model'], int(parameters['RAACCluster3B']), int(parameters['Ktuple']),
            int(parameters['GapLamada'])), 'type3B.type1(testing_data, "%s", %d, %d, %d)' % (
                       parameters['PseKRAAC_Model'], int(parameters['RAACCluster3B']), int(parameters['Ktuple']),
                       int(parameters['GapLamada']))],
        'type4': ['type4.type1(training_data, "%s", %d, %d, %d)' % (
            parameters['PseKRAAC_Model'], int(parameters['RAACCluster4']), int(parameters['Ktuple']),
            int(parameters['GapLamada'])), 'type4.type1(testing_data, "%s", %d, %d, %d)' % (
                      parameters['PseKRAAC_Model'], int(parameters['RAACCluster4']), int(parameters['Ktuple']),
                      int(parameters['GapLamada']))],
        'type5': ['type5.type1(training_data, "%s", %d, %d, %d)' % (
            parameters['PseKRAAC_Model'], int(parameters['RAACCluster5']), int(parameters['Ktuple']),
            int(parameters['GapLamada'])), 'type5.type1(testing_data, "%s", %d, %d, %d)' % (
                      parameters['PseKRAAC_Model'], int(parameters['RAACCluster5']), int(parameters['Ktuple']),
                      int(parameters['GapLamada']))],
        'type6A': ['type6A.type1(training_data, "%s", %d, %d, %d)' % (
            parameters['PseKRAAC_Model'], int(parameters['RAACCluster6A']), int(parameters['Ktuple']),
            int(parameters['GapLamada'])), 'type6A.type1(testing_data, "%s", %d, %d, %d)' % (
                       parameters['PseKRAAC_Model'], int(parameters['RAACCluster6A']), int(parameters['Ktuple']),
                       int(parameters['GapLamada']))],
        'type6B': ['type6B.type1(training_data, "%s", %d, %d, %d)' % (
            parameters['PseKRAAC_Model'], int(parameters['RAACCluster6B']), int(parameters['Ktuple']),
            int(parameters['GapLamada'])), 'type6B.type1(testing_data, "%s", %d, %d, %d)' % (
                       parameters['PseKRAAC_Model'], int(parameters['RAACCluster6B']), int(parameters['Ktuple']),
                       int(parameters['GapLamada']))],
        'type6C': ['type6C.type1(training_data, "%s", %d, %d, %d)' % (
            parameters['PseKRAAC_Model'], int(parameters['RAACCluster6C']), int(parameters['Ktuple']),
            int(parameters['GapLamada'])), 'type6C.type1(testing_data, "%s", %d, %d, %d)' % (
                       parameters['PseKRAAC_Model'], int(parameters['RAACCluster6C']), int(parameters['Ktuple']),
                       int(parameters['GapLamada']))],
        'type7': ['type7.type1(training_data, "%s", %d, %d, %d)' % (
            parameters['PseKRAAC_Model'], int(parameters['RAACCluster7']), int(parameters['Ktuple']),
            int(parameters['GapLamada'])), 'type7.type1(testing_data, "%s", %d, %d, %d)' % (
                      parameters['PseKRAAC_Model'], int(parameters['RAACCluster7']), int(parameters['Ktuple']),
                      int(parameters['GapLamada']))],
        'type8': ['type8.type1(training_data, "%s", %d, %d, %d)' % (
            parameters['PseKRAAC_Model'], int(parameters['RAACCluster8']), int(parameters['Ktuple']),
            int(parameters['GapLamada'])), 'type8.type1(testing_data, "%s", %d, %d, %d)' % (
                      parameters['PseKRAAC_Model'], int(parameters['RAACCluster8']), int(parameters['Ktuple']),
                      int(parameters['GapLamada']))],
        'type9': ['type9.type1(training_data, "%s", %d, %d, %d)' % (
            parameters['PseKRAAC_Model'], int(parameters['RAACCluster9']), int(parameters['Ktuple']),
            int(parameters['GapLamada'])), 'type9.type1(testing_data, "%s", %d, %d, %d)' % (
                      parameters['PseKRAAC_Model'], int(parameters['RAACCluster9']), int(parameters['Ktuple']),
                      int(parameters['GapLamada']))],
        'type10': ['type10.type1(training_data, "%s", %d, %d, %d)' % (
            parameters['PseKRAAC_Model'], int(parameters['RAACCluster10']), int(parameters['Ktuple']),
            int(parameters['GapLamada'])), 'type10.type1(testing_data, "%s", %d, %d, %d)' % (
                       parameters['PseKRAAC_Model'], int(parameters['RAACCluster10']), int(parameters['Ktuple']),
                       int(parameters['GapLamada']))],
        'type11': ['type11.type1(training_data, "%s", %d, %d, %d)' % (
            parameters['PseKRAAC_Model'], int(parameters['RAACCluster11']), int(parameters['Ktuple']),
            int(parameters['GapLamada'])), 'type11.type1(testing_data, "%s", %d, %d, %d)' % (
                       parameters['PseKRAAC_Model'], int(parameters['RAACCluster11']), int(parameters['Ktuple']),
                       int(parameters['GapLamada']))],
        'type12': ['type12.type1(training_data, "%s", %d, %d, %d)' % (
            parameters['PseKRAAC_Model'], int(parameters['RAACCluster12']), int(parameters['Ktuple']),
            int(parameters['GapLamada'])), 'type12.type1(testing_data, "%s", %d, %d, %d)' % (
                       parameters['PseKRAAC_Model'], int(parameters['RAACCluster12']), int(parameters['Ktuple']),
                       int(parameters['GapLamada']))],
        'type13': ['type13.type1(training_data, "%s", %d, %d, %d)' % (
            parameters['PseKRAAC_Model'], int(parameters['RAACCluster13']), int(parameters['Ktuple']),
            int(parameters['GapLamada'])), 'type13.type1(testing_data, "%s", %d, %d, %d)' % (
                       parameters['PseKRAAC_Model'], int(parameters['RAACCluster13']), int(parameters['Ktuple']),
                       int(parameters['GapLamada']))],
        'type14': ['type14.type1(training_data, "%s", %d, %d, %d)' % (
            parameters['PseKRAAC_Model'], int(parameters['RAACCluster14']), int(parameters['Ktuple']),
            int(parameters['GapLamada'])), 'type14.type1(testing_data, "%s", %d, %d, %d)' % (
                       parameters['PseKRAAC_Model'], int(parameters['RAACCluster14']), int(parameters['Ktuple']),
                       int(parameters['GapLamada']))],
        'type15': ['type15.type1(training_data, "%s", %d, %d, %d)' % (
            parameters['PseKRAAC_Model'], int(parameters['RAACCluster15']), int(parameters['Ktuple']),
            int(parameters['GapLamada'])), 'type15.type1(testing_data, "%s", %d, %d, %d)' % (
                       parameters['PseKRAAC_Model'], int(parameters['RAACCluster15']), int(parameters['Ktuple']),
                       int(parameters['GapLamada']))],
        'type16': ['type16.type1(training_data, "%s", %d, %d, %d)' % (
            parameters['PseKRAAC_Model'], int(parameters['RAACCluster16']), int(parameters['Ktuple']),
            int(parameters['GapLamada'])), 'type16.type1(testing_data, "%s", %d, %d, %d)' % (
                       parameters['PseKRAAC_Model'], int(parameters['RAACCluster16']), int(parameters['Ktuple']),
                       int(parameters['GapLamada']))],
    }

    # read fasta sequence and specify cmd
    fastas = []
    cmd_coding = {}
    if parameters['Sequence_Type'] in ('DNA', 'RNA'):
        fastas = read_fasta_sequences.read_nucleotide_sequences(parameters['Sequence_File'])
        cmd_coding = dna_cmd_coding
    if parameters['Sequence_Type'] == 'Protein':
        fastas = read_fasta_sequences.read_protein_sequences(parameters['Sequence_File'])
        cmd_coding = protein_cmd_coding


    kw = {'nclusters': 3, 'sof': 'sample', 'order': ''}
    kw['order'] = 'ACGT' if parameters['Sequence_Type'] == 'DNA' or parameters['Sequence_Type'] == 'RNA' else 'ACDEFGHIKLMNPQRSTVWY'

    # divide training and testing data
    training_data = []
    testing_data = []
    PSTNP_training_data = []
    PSTNP_testing_data = []
    classes = set()
    for sequence in fastas:
        if sequence[3] == 'training':
            training_data.append(sequence)
            PSTNP_training_data.append(sequence)
            PSTNP_training_data.append([sequence[0], sequence[1], sequence[2], 'testing'])
            PSTNP_testing_data.append(sequence)
            classes.add(sequence[2])
        else:
            testing_data.append(sequence)
            PSTNP_testing_data.append(sequence)

    # get property for AAindex, NMBroto, Geary, Moran
    props = parameters['AAindex'].split(';') if parameters['AAindex'] != '' else ['CIDH920105', 'BHAR880101',
                                                                                  'CHOC760101', 'BIGC670101',
                                                                                  'CHAM810101', 'DAYM780201']
    # get property for ACC descriptors and Pse descriptors
    my_property_name, my_property_value, my_kmer, my_lamada, my_weight = 0, 0, 0, 0, 0

    #
    method_array = parameters['Method'].split(';')
    ML_array = parameters['ML'].split(';')

    for ML in ML_array:
        # parameters for binary classification task
        CV_records = {}
        IND_records = {}
        CV_record_metrics = {}
        IND_record_metrics = {}
        # parameters for muti-classification task
        muti_record_metrics = {}
        for method in method_array:
            print('Building model for %s - %s' %(ML, method))
            if method in ('DAC', 'DCC', 'DACC', 'TAC', 'TCC', 'TACC'):
                my_property_name, my_property_value, my_kmer = check_parameters.check_acc_arguments_pipeline(parameters,
                                                                                                             method)
            if method in ('PseDNC', 'PseKNC', 'PCPseDNC', 'PCPseTNC', 'SCPseDNC', 'SCPseTNC'):
                my_property_name, my_property_value, my_lamada, my_weight = check_parameters.check_Pse_arguments_pipeline(
                    parameters, method, fastas)

            training_code = eval(cmd_coding[method][0])
            testing_dode = []
            save_file.save_file(training_code, format=parameters['Output_Format'], file='training_code.txt')
            if len(testing_data) > 0:
                testing_dode = eval(cmd_coding[method][1])
                save_file.save_file(testing_dode, format=parameters['Output_Format'], file='testing_code.txt')

            X, y, independent = 0, 0, np.array([])
            X, y = read_code_ml.read_code('training_code.txt', format=parameters['Output_Format'])
            if len(testing_data) > 0:
                ind_X, ind_y = read_code_ml.read_code('testing_code.txt', format=parameters['Output_Format'])
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

            # training model
            para_info, cv_res, ind_res = 0, 0, 0
            if ML == 'SVM':
                para_info, cv_res, ind_res = SVM.SVM_Classifier(X, y, indep=independent, fold=fold, batch=0.8, auto=auto, kernel=kernel, gamma=gamma, C=cost)
            if ML == 'RF':
                para_info, cv_res, ind_res = RF.RF_Classifier(X, y, indep=independent, fold=fold, n_trees=n_trees)
            if ML == 'LR':
                para_info, cv_res, ind_res = LR.LR_Classifier_binary(X, y, indep=independent, fold=fold)
            if ML == 'KNN':
                para_info, cv_res, ind_res = KNN.KNN_Classifier_binary(X, y, indep=independent, fold=fold, n_neighbors=top_K)
            if ML == 'MLP':
                para_info, cv_res, ind_res = MLP.MLP_Classifier(X, y, indep=independent, fold=fold, hidden_layer_size=hidden_layer_size)

            # binary classification task
            if len(classes) == 2:
                CV_records['%s-%s' %(ML, method)] = cv_res
                IND_records['%s-%s' %(ML, method)] = ind_res
                cv_metrics = calculate_prediction_metrics.calculate_metrics_cv(cv_res, label_column=0, score_column=2)
                my_metrics = {'Sensitivity': 0, 'Specificity': 0, 'Accuracy': 0, 'MCC': 0, 'Recall': 0, 'Precision': 0, 'F1-score': 0}
                for i in range(len(cv_metrics)):
                    for key in my_metrics:
                        try:
                            my_metrics[key] += cv_metrics[i][key]
                        except TypeError as e:
                            pass
                for key in my_metrics:
                    my_metrics[key] /= len(cv_metrics)
                CV_record_metrics['%s-%s' %(ML, method)] = my_metrics

                if len(testing_data) > 0:
                    ind_metrics = calculate_prediction_metrics.calculate_metrics(ind_res[:, 0], ind_res[:, 2])
                    IND_record_metrics['%s-%s' %(ML, method)] = ind_metrics

            # muti-classification task
            if len(classes) >= 3:
                muti_record_metrics['%s-%s' %(ML, method)] = [calculate_prediction_metrics.calculate_metrics_cv_muti_V2(cv_res, classes, label_column=0)]
                if len(testing_data) > 0:
                    muti_record_metrics['%s-%s' % (ML, method)].append(calculate_prediction_metrics.calculate_metrics_ind_muti_V2(ind_res, classes, label_column=0))
                else:
                    muti_record_metrics['%s-%s' % (ML, method)].append('NA')

        if len(classes) == 2:
            draw_plot.plot_mean_roc_cv(CV_records, '%s_CV_ROC.png' % ML, label_column=0, score_column=2)
            draw_plot.plot_mean_prc_CV(CV_records, '%s_CV_PRC.png' % ML, label_column=0, score_column=2)
            save_file.save_prediction_metrics_2_classes(CV_record_metrics, 'Metrics_CV.txt')

            if len(testing_data) > 0:
                draw_plot.plot_roc_muti_ind(IND_records, '%s_IND_ROC.png' % ML, label_column=0, score_column=2)
                draw_plot.plot_prc_muti_ind(IND_records, '%s_IND_PRC.png' % ML, label_column=0, score_column=2)
                save_file.save_prediction_metrics_2_classes(IND_record_metrics, 'Metrics_IND.txt')

        if len(classes) >= 3:
            save_file.save_prediction_metrics_muti_V2(muti_record_metrics, '%s_Metrics_CV.txt' %ML)

