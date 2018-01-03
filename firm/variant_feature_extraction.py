from __future__ import absolute_import, division, print_function

import numpy as np
import pandas as pd

from Bio.SubsMat import MatrixInfo

from .util import log
from .aa_scales import aa_scales
from .apply_scale import apply_scale
from .pfam_scores import get_ref_and_alt_scores
from .ml.feature_engineering import CombinerFeatureExtractor, SimpleFeatureExtractor, OneHotEncodingFeatureExtractor, \
        CounterFeatureExtractor

class FeatureExtractionSetup(object):
    def __init__(self, geneffect_setup):
        self.common_pfam_clans = geneffect_setup.pfam_data.common_pfam_clans
        self.all_uniprot_boolean_tracks = geneffect_setup.all_uniprot_boolean_tracks
        self.all_uniprot_categorical_tracks = geneffect_setup.all_uniprot_categorical_tracks

def get_snp_effect_feature_extractor(setup):
    
    protein_context_feature_extractor = CombinerFeatureExtractor(
            [create_global_length_and_location_feature_extractor()] + \
            [create_context_aa_and_scales_feature_extractor(context_name, left_size, right_size) for context_name, left_size, \
                    right_size in NEIGHBORHOOD_OPTIONS] + \
            [create_context_track_counts_feature_extractor(context_name, left_size, right_size, setup.all_uniprot_boolean_tracks, \
                    setup.all_uniprot_categorical_tracks) for context_name, left_size, right_size in NEIGHBORHOOD_OPTIONS]
    )

    aa_substitution_feature_extractor = CombinerFeatureExtractor(
            [
                OneHotEncodingFeatureExtractor('ref_aa_', ALL_AAS, lambda snp_effect: snp_effect.ref_aa),
                OneHotEncodingFeatureExtractor('alt_aa_', ALL_AAS, lambda snp_effect: snp_effect.alt_aa),
            ] + \
            [
                create_blosum_scores_feature_extractor(),
                create_track_distances_and_hits_feature_extractor(setup.all_uniprot_boolean_tracks, setup.all_uniprot_categorical_tracks),
            ] + \
            [create_scale_diff_feature_extractor(scale_name, scale_values) for scale_name, scale_values in AA_SCALE_ITEMS]
    )

    return CombinerFeatureExtractor([
        protein_context_feature_extractor,
        aa_substitution_feature_extractor,
        create_pfam_domain_feature_extractor(setup.common_pfam_clans),
    ])
    
def create_pfam_domain_feature_extractor(common_pfam_clans):

    def extraction_function(snp_effect):
        
        domain_record = snp_effect.affected_gene.uniprot_record.get_closest_pfam_domain(snp_effect.protein_coordinate)
        found = (domain_record is not None)
        
        if found:
            distance = domain_record['distance']
            hit = (distance == 0)
            relative_position = calc_relative_position_to_pfam_domain(domain_record, snp_effect.protein_coordinate)
            length = domain_record['length']
            reported_ref_score = domain_record['bit_score']
            clan = domain_record['clan']
            clan_exists = pd.notnull(clan)
        else:
            distance = INFINITE_DISTANCE_VALUE
            hit = False
            relative_position = 0
            length = 0
            reported_ref_score = 0
            clan = None
            clan_exists = False
            
        ref_score = 0
        alt_score = 0
            
        if hit:
            try:
                ref_score, alt_score = get_ref_and_alt_scores(domain_record['hmm_name'], snp_effect.affected_gene.uniprot_record.seq, \
                        snp_effect.protein_coordinate, snp_effect.ref_aa, snp_effect.alt_aa)
            except:
                log('%s: failed calculating pfam scores.' % snp_effect)
            
        absolute_score_change = alt_score - ref_score
        relative_score_change = absolute_score_change / ref_score if ref_score > 0 else 0
        destroyed_hit = (ref_score > 0) and (alt_score == 0)
        
        return found, hit, distance, relative_position, length, reported_ref_score, ref_score, alt_score, \
                absolute_score_change, relative_score_change, destroyed_hit, clan_exists, clan
    
    features = ['pfam_domain_found', 'pfam_domain_hit', 'pfam_domain_distance', 'pfam_domain_relative_position', \
            'pfam_domain_length', 'pfam_domain_reported_ref_score', 'pfam_domain_ref_score', 'pfam_domain_alt_score', \
            'pfam_domain_absolute_score_change', 'pfam_domain_relative_score_change', 'pfam_domain_destroyed_hit', \
            'pfam_domain_clan_exists'] + \
            [OneHotEncodingFeatureExtractor('pfam_domain_clan_', common_pfam_clans)]
    return SimpleFeatureExtractor(features, extraction_function)
    
def create_track_distances_and_hits_feature_extractor(all_uniprot_boolean_tracks, all_uniprot_categorical_tracks):

    def extraction_function(snp_effect):
    
        uniprot_record = snp_effect.affected_gene.uniprot_record
        boolean_track_features = []
        categorical_track_features = []
    
        for track_name in all_uniprot_boolean_tracks:
            if track_name in uniprot_record.boolean_tracks:
                
                track = uniprot_record.boolean_tracks[track_name]
                distance = track.calc_distance(snp_effect.protein_coordinate, len(uniprot_record))
                
                if distance is None:
                    distance = INFINITE_DISTANCE_VALUE
                    
                hit = (distance == 0)
                miss = (not hit and not track.get_values() <= set([False]))
                boolean_track_features.extend([distance, hit, miss])
            else:
                boolean_track_features.extend([INFINITE_DISTANCE_VALUE, False, False])
        
        for track_name, values in all_uniprot_categorical_tracks:
            for value in (values + [None]):
                if track_name in uniprot_record.categorical_tracks:
                
                    track = uniprot_record.categorical_tracks[track_name]
                    distance = track.calc_distance(snp_effect.protein_coordinate, len(uniprot_record), value)

                    if distance is None:
                        distance = INFINITE_DISTANCE_VALUE

                    hit = (distance == 0)
                    miss = (not hit and not track.get_values() <= set([None]))
                    categorical_track_features.extend([distance, hit, miss])
                else:
                    if value is None:
                        categorical_track_features.extend([0, True, False])
                    else:
                        categorical_track_features.extend([INFINITE_DISTANCE_VALUE, False, True])
            
        return boolean_track_features + categorical_track_features
    
    feature_names = ['track_%s_%s' % (track_name, feature_type) for track_name in \
            all_uniprot_boolean_tracks for feature_type in ['distance', 'hit', 'miss']] + \
            ['track_%s_%s_%s' % (track_name, value, feature_type) for track_name, values in \
            all_uniprot_categorical_tracks for value in (values + [None]) for feature_type in ['distance', 'hit', 'miss']]
    return SimpleFeatureExtractor(feature_names, extraction_function)
    
def create_global_length_and_location_feature_extractor():
    
    def extraction_function(snp_effect):
        protein_length = len(snp_effect.affected_gene.uniprot_record.seq)
        variant_position = snp_effect.protein_coordinate
        variant_relative_position = variant_position / protein_length
        return protein_length, variant_position, variant_relative_position
    
    feature_names = ['protein_length', 'variant_position', 'variant_relative_position']
    return SimpleFeatureExtractor(feature_names, extraction_function)
    
def create_context_aa_and_scales_feature_extractor(context_name, left_size, right_size):

    AGGREGATION_METHODS = [np.average, np.std]
    AGGREGATION_METHOD_NAMES = ['avg', 'std']
    SCALE_NAMES, SCALE_VALUES = zip(*AA_SCALE_ITEMS)
    
    def extraction_function(snp_effect):
        context_seq = get_neighborhood(snp_effect.affected_gene.uniprot_record.seq, snp_effect.protein_coordinate, \
                left_size, right_size)
        scale_arrays = [apply_scale(context_seq, scale_values) for scale_values in SCALE_VALUES]
        return [context_seq] + [_replace_nan(aggregator(scale_array)) for scale_array in scale_arrays for \
                aggregator in AGGREGATION_METHODS]
    
    features = [CounterFeatureExtractor('context_%s_aa_' % context_name, ALL_AAS_AND_DUMMY)] + \
            ['context_%s_scale_%s_%s' % (context_name, scale_name, aggregator_name) for scale_name in SCALE_NAMES for \
                    aggregator_name in AGGREGATION_METHOD_NAMES]
    return SimpleFeatureExtractor(features, extraction_function)
    
def create_context_track_counts_feature_extractor(context_name, left_size, right_size, all_uniprot_boolean_tracks, all_uniprot_categorical_tracks):
    
    def extraction_function(snp_effect):
        
        uniprot_record = snp_effect.affected_gene.uniprot_record
        total_length = len(uniprot_record)
        
        if left_size is None:
            start_bound = 1
        else:
            start_bound = max(snp_effect.protein_coordinate - left_size, 1)
            
        if right_size is None:
            end_bound = total_length
        else:
            end_bound = min(snp_effect.protein_coordinate + right_size, total_length)
        
        context_length = end_bound - start_bound + 1
        boolean_track_counts = []
        categorical_track_counts = []
        
        for track_name in all_uniprot_boolean_tracks:
            if track_name in uniprot_record.boolean_tracks:
                track_count = uniprot_record.boolean_tracks[track_name].count_values(start_bound = start_bound, \
                        end_bound = end_bound).get(True, 0)
                boolean_track_counts.append(track_count)
            else:
                boolean_track_counts.append(0)
                
        for track_name, values in all_uniprot_categorical_tracks:
            if track_name in uniprot_record.categorical_tracks:
                track_counts = uniprot_record.categorical_tracks[track_name].count_values(total_length = total_length, \
                        start_bound = start_bound, end_bound = end_bound)
                categorical_track_counts.extend([track_counts.get(value, 0) for value in (values + [None])])
            else:
                categorical_track_counts.extend(len(values) * [0] + [context_length])
                
        return boolean_track_counts + categorical_track_counts
                
    feature_names = ['context_%s_track_%s_count' % (context_name, track_name) for track_name in all_uniprot_boolean_tracks] + \
            ['context_%s_track_%s_count_%s' % (context_name, track_name, value) for track_name, values in all_uniprot_categorical_tracks for \
                    value in (values + [None])]
    return SimpleFeatureExtractor(feature_names, extraction_function)
    
def create_scale_diff_feature_extractor(scale_name, scale_values):
    
    def extraction_function(snp_effect):
        ref_value = scale_values.get(snp_effect.ref_aa, 0)
        alt_value = scale_values.get(snp_effect.alt_aa, 0)
        absolute_diff = alt_value - ref_value
        relative_diff = 0 if ref_value == 0 else absolute_diff / ref_value
        return ref_value, alt_value, absolute_diff, relative_diff
    
    feature_name_templates = ['scale_%s_ref', 'scale_%s_alt', 'scale_%s_absolute_diff', 'scale_%s_relative_diff']
    feature_names = [feature_name_template % scale_name for feature_name_template in feature_name_templates]
    return SimpleFeatureExtractor(feature_names, extraction_function)
    
def create_blosum_scores_feature_extractor():
    
    blosum_names, blosum_matrices = zip(*BLOSUM_OPTIONS)
    
    def extraction_function(snp_effect):
        return [_get_blosum_score(blosum_matrix, snp_effect.ref_aa, snp_effect.alt_aa) for blosum_matrix in blosum_matrices]
    
    return SimpleFeatureExtractor(blosum_names, extraction_function)
    
def get_neighborhood(seq, position, left_size, right_size):
    
    # Convert from the usual 1-based indexing to Python 0-based one.
    position -= 1
    
    if left_size is None:
        left_seq = seq[:position]
    else:
        start = max(0, position - left_size)
        left_seq = seq[start:position]
        left_seq = (left_size - len(left_seq)) * '_' + left_seq
        
    if right_size is None:
        right_seq = seq[(position + 1):]
    else:
        right_seq = seq[(position + 1):(position + right_size + 1)]
        right_seq = right_seq + (right_size - len(right_seq)) * '_'

    return left_seq + seq[position] + right_seq
    
def calc_relative_position_to_pfam_domain(pfam_record, position):

    start = pfam_record['alignment_start']
    end = pfam_record['alignment_end']

    if start >= end:
        return 0
    else:
        return (position - start) / (end - start)
    
def _get_blosum_score(blosum_matrix, aa1, aa2):
    if (aa1, aa2) in blosum_matrix:
        return blosum_matrix[(aa1, aa2)]
    else:
        return blosum_matrix[(aa2, aa1)]
    
def _replace_nan(value, nan_predicate = np.isnan, replace_with = 0):
    if nan_predicate(value):
        return replace_with
    else:
        return value

    
INFINITE_DISTANCE_VALUE = 1000

ALL_AAS = list('ACDEFGHIKLMNPQSRTVWY')
ALL_AAS_AND_DUMMY = ALL_AAS + ['_']

BLOSUM_OPTIONS = [
    ('BLOSUM45', MatrixInfo.blosum45),
    ('BLOSUM62', MatrixInfo.blosum62),
    ('BLOSUM80', MatrixInfo.blosum80),
]

NEIGHBORHOOD_OPTIONS = [
    ('global', None, None),
    ('5', 5, 5),
    ('left_10', 10, 0),
    ('right_10', 0, 10),
    ('left', None, 0),
    ('right', 0, None),
]

AA_SCALE_ITEMS = aa_scales.items()
