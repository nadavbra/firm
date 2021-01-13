from collections import defaultdict

import geneffect

from .util import log

def get_uniprot_id(gene):
    return gene.uniprot_record.id

def determine_extended_gene_effects_and_scores(variants, firm_predict_batch_adjusted_proba_function, gene_indexer = get_uniprot_id):
    
    '''
    Determine the gene effects and effect scores of a given list of variants, using an extended policy that also assigns scores to non-missense
    variants (using siplistic, rule-based criteria). In short, loss-of-function variants (including nonsense, frameshift and canonical splicing
    variants) will be assigned a score of 0, synonymous variants will be assigned a score of 1, in-frame indels will be assigned a score according
    to get_complex_aa_change_effect_score, and missense variants will be assigned a score by the FIRM classifier (after reversing between 0 to 1,
    i.e. applying the transformation x -> 1 - x). For more details, see our full paper (which is linked from the GitHub repository of the project).
    Note that the effect scores produced by this tool are reversed to the scores produced by FIRM's programmatic API (here 0 denotes complete damage,
    and 1 no effect).
    @param variants (an iterable of geneffect.Variant objects): The variants to determine the effects for.
    @param firm_predict_batch_adjusted_proba_function (function): A function that receives a list of geneffect.SnpCDSGeneEffect objects representing
    missense variants, and returns their predicted effect scores as numbers between 0 to 1. Here the scores are expected to be consistent with those
    produced by FIRM's programmatic API (i.e. 0 indicates the least damage), so this function will flip the output (turning x into 1 - x). For example,
    firm_classifier.predict_adjusted_proba could be a valid firm_predict_batch_adjusted_proba_function.
    @param gene_indexer (function): A function that receives a geneffect.Gene object and returns a unique identifier (by default will use the UniProt
    ID).
    @return: A list of dictionaries mapping the unique gene IDs (as defined by gene_indexer) to tuples of (effect_str, effect_score), where effect_str
    is a string explaining a variant effect and effect_score is a variant's effect score. Each dictionary in the list should indicate the results of
    the corresponding variant in the input list. 
    '''
        
    variants_gene_effects_and_scores = [defaultdict(list) for _ in range(len(variants))]
    missense_effects = []
    missense_effects_meta_data = []
    
    for variant_index, (variant, variant_gene_effects_and_scores) in enumerate(zip(variants, variants_gene_effects_and_scores)):
        
        if variant is None:
            continue
    
        for splicing_gene_effect in variant.splicing_gene_effects:
            assert isinstance(splicing_gene_effect, geneffect.variant_processing.CanonicalSplicingGeneEffect), 'Unknown splicing gene effect %s' % \
                    type(splicing_gene_effect)
            variant_gene_effects_and_scores[gene_indexer(splicing_gene_effect.affected_gene)].append((str(splicing_gene_effect), 0.0))
        
        for cds_gene_effect in variant.cds_gene_effects:
        
            gene_index = gene_indexer(cds_gene_effect.affected_gene)
            
            if cds_gene_effect.destroys_start_codon() or cds_gene_effect.is_frameshift or cds_gene_effect.introduced_stop_codon:
                variant_gene_effects_and_scores[gene_index].append((str(cds_gene_effect), 0.0))
            elif cds_gene_effect.is_snp():
                if cds_gene_effect.is_synonymous():
                    variant_gene_effects_and_scores[gene_index].append((str(cds_gene_effect), 1.0))
                elif cds_gene_effect.is_missense():
                    missense_effects.append(cds_gene_effect)
                    missense_effects_meta_data.append((variant_index, gene_index))
                else:
                    raise ValueError('Unexpected CDS gene effect type: %s' % cds_gene_effect)
            else:
                variant_gene_effects_and_scores[gene_index].append(('inframe_indel', get_complex_aa_change_effect_score(cds_gene_effect.ref_protein_seq, \
                        cds_gene_effect.alt_protein_seq)))
                        
    log('Applying FIRM classifier on %d missense effects...' % len(missense_effects))
    missense_effect_scores = 1 - firm_predict_batch_adjusted_proba_function(missense_effects)        

    for missense_effect, effect_score, (variant_index, gene_index) in zip(missense_effects, missense_effect_scores, missense_effects_meta_data):
        variants_gene_effects_and_scores[variant_index][gene_index].append((str(missense_effect), effect_score))
    
    return list(map(dict, variants_gene_effects_and_scores))
    
def get_complex_aa_change_effect_score(ref_aa_seq, alt_aa_seq):
        
    '''
    Calculates a heuristic-based effect score for a complex (in-frame) amino-acid alteration variant. The score is based solely on the number
    of substituted, inserted and deleted amino-acids, according to ClinVar's statistics (the probability of such variants to be considered
    "pathogenic"). For more details, see our full paper (which is linked from the GitHub repository of the project).
    '''
    
    while len(ref_aa_seq) > 0 and len(alt_aa_seq) > 0 and ref_aa_seq[0] == alt_aa_seq[0]:
        ref_aa_seq = ref_aa_seq[1:]
        alt_aa_seq = alt_aa_seq[1:]
        
    while len(ref_aa_seq) > 0 and len(alt_aa_seq) > 0 and ref_aa_seq[-1] == alt_aa_seq[-1]:
        ref_aa_seq = ref_aa_seq[:-1]
        alt_aa_seq = alt_aa_seq[:-1]
        
    ref_len = len(ref_aa_seq)
    alt_len = len(alt_aa_seq)
    
    n_subs = min(ref_len, alt_len)
    n_inserts = max(alt_len - n_subs, 0)
    n_dels = max(ref_len - n_subs, 0)
    
    return ((1 - 0.71) ** n_subs) * ((1 - 0.44) ** n_inserts) * ((1 - 0.85) ** n_dels)
