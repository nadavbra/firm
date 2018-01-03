from __future__ import absolute_import, division, print_function

import os

SRC_DIR = os.path.dirname(__file__)

# PFAM_HMM_DIR should point to a directory with the all the necessary pfam HMM profiles.
# It should containt the files: Pfam-A.hmm, Pfam-A.hmm.dat, Pfam-A.hmm.h3f, Pfam-A.hmm.h3i, Pfam-A.hmm.h3m, Pfam-A.hmm.h3p, active_site.dat.
# It should also contain a subdirectory named profiles/. In the meantime, this subdirectory can be left empty. As you run FIRM on new proteins
# and encounter new pfam domains, files will automatically be created within this directory. For every domain name *, the following files will
# automatically be created: *.hmm, *.hmm.h3f, *.hmm.h3i, *.hmm.h3m, *.hmm.h3p. 
# Altogether, you should have ~3.1 GB of available storage for this directory.
# 
# In order to prepare the files within this directory (and allow the module pfam_scores.py to use them later), you will need to have HMMER3
# installed in your environment. If it's not already installed, run the following (from within a temporary directory such as /tmp):
# >> wget http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz
# >> tar -zxvf hmmer-3.1b2-linux-intel-x86_64.tar.gz
# >> cd hmmer-3.1b2-linux-intel-x86_64
# >> ./configure
# >> make
# >> sudo make install
#
# In order to download the relevant data and prepare the files required within this directory, run the following (from within your
# PFAM_HMM_DIR directory):
# >> wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.hmm.gz
# >> gzip -d Pfam-A.hmm.gz
# >> wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.hmm.dat.gz
# >> gzip -d Pfam-A.hmm.dat.gz
# >> hmmpress Pfam-A.hmm
# >> wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/active_site.dat.gz
# >> gzip -d active_site.dat.gz
PFAM_HMM_DIR = os.path.join(os.path.expanduser('~'), 'data/pfam/hmm')
PFAM_HMM_PROFILES_DIR = os.path.join(PFAM_HMM_DIR, 'profiles')

# Paths for temporary files created by pfam_scores.py (%s indicates process ID, to allow running it with multiprocessing)
TMP_FASTA_FILE_PATH = '/tmp/temp_pfam_seq_pid_%s.fasta'
TMP_DOMAIN_RESULTS_FILE_PATH = '/tmp/temp_pfam_domain_results_pid_%s'

# Unless you switch the file, do not change that.
TRAINED_CLASSIFER_DUMP_FILE_PATH = os.path.join(SRC_DIR, 'data/classifier.pkl')
