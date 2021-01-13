import os
import re
import subprocess
import traceback

from .config import PFAM_HMM_DIR, PFAM_HMM_PROFILES_DIR, TMP_FASTA_FILE_PATH, TMP_DOMAIN_RESULTS_FILE_PATH
from .util import log

def create_missing_hmm_profiles(hmm_names):

    hmm_names = list(sorted(set(hmm_names)))
    missing_hmm_names = list(filter(_is_hmm_profile_missing, hmm_names))
    n_failures = 0
    
    if len(missing_hmm_names) > 0:
        
        log('%d of %d HMM profiles are missing. FIRM will create them before proceeding (this is a one-time operation)...' % \
                (len(missing_hmm_names), len(hmm_names)))
        
        for hmm_name in missing_hmm_names:
            try:
                _create_hmm_profile(hmm_name)
            except:
                n_failures += 1
                log('Failed creating HMM %s:\n%s' % (hmm_name, traceback.format_exc()))
                
        if n_failures == len(hmm_names) and len(hmm_names) >= 5:
            raise Exception('Failed creating all %d HMM profiles. You probably haven\'t installed HMMER3 properly.' % len(hmm_names))
        elif n_failures > 0:
            log('Failed creating %d of %d missing HMM profiles.' % (n_failures, len(missing_hmm_names)))
                
def does_hmm_profile_exist(hmm_name):
    hmm_profile_file_path = _get_hmm_profile_file_path(hmm_name)
    return os.path.exists(hmm_profile_file_path) and os.path.getsize(hmm_profile_file_path) > 0

def get_hmm_ref_and_alt_scores(hmm_name, ref_seq, position, ref_aa, alt_aa):
    position_zero_based = position - 1
    assert ref_seq[position_zero_based] == ref_aa
    alt_seq = ref_seq[:position_zero_based] + alt_aa + ref_seq[(position_zero_based + 1):]
    return get_hmm_score(hmm_name, ref_seq, position), get_hmm_score(hmm_name, alt_seq, position)

def get_hmm_score(hmm_name, seq, position):
    _override_temp_fasta_file(seq)
    _run_process('hmmscan --noali --notextw --cut_ga --domtblout %s %s %s' % (_get_tmp_domain_results_file_path(), _get_hmm_profile_file_path(hmm_name), \
            _get_tmp_fasta_file_path()))
    return _get_domain_result_score(position)
    
def _is_hmm_profile_missing(hmm_name):
    return not os.path.exists(_get_hmm_profile_file_path(hmm_name))
    
def _create_hmm_profile(hmm_name):
    log('Creating HMM profile: %s' % hmm_name)
    hmm_profile_file_path = _get_hmm_profile_file_path(hmm_name)
    _run_process('hmmfetch %s %s > %s' % (os.path.join(PFAM_HMM_DIR, 'Pfam-A.hmm'), hmm_name, hmm_profile_file_path))
    _run_process('hmmpress %s' % hmm_profile_file_path)

def _get_domain_result_score(position):

    relevant_scores = [score for start_index, end_index, score in _parse_domain_results() if position >= start_index and position <= end_index]

    if len(relevant_scores) == 0:
        return 0.0
    else:
        return max(relevant_scores)

def _parse_domain_results():
    with open(_get_tmp_domain_results_file_path(), 'r') as f:
        for line in f.readlines():
            if not line.startswith('#'):
                raw_fields = re.split(r'\s+', line.strip())
                score = float(raw_fields[7])
                start_index = int(raw_fields[17])
                end_index = int(raw_fields[18])
                yield start_index, end_index, score

def _override_temp_fasta_file(seq):
    with open(_get_tmp_fasta_file_path(), 'w') as f:
        f.write('>seq' + '\n')
        f.write(str(seq) + '\n')

def _run_process(command_line):
    
    process = subprocess.Popen(command_line, stdout = subprocess.PIPE, shell = True)
    returncode = process.wait()
    
    if returncode != 0:
        raise IOError('Process failed with return-code %d: %s' % (returncode, command_line))
        
def _get_hmm_profile_file_path(hmm_name):
    return os.path.join(PFAM_HMM_PROFILES_DIR, '%s.hmm' % hmm_name)
        
def _get_tmp_fasta_file_path():
    return TMP_FASTA_FILE_PATH % os.getpid()
    
def _get_tmp_domain_results_file_path():
    return TMP_DOMAIN_RESULTS_FILE_PATH % os.getpid()
