What is FIRM?
==========

FIRM (Functional Impact Rating at the Molecular-level) is a machine-learning model for predicting the functional impact of genetic variants (currently focusing only on missense variants affecting protein sequences). FIRM assigns variants with scores between 0 (harmless) to 1 (harmful). The model can be easily installed and run as a standalone Python module, and it is sufficiently fast (and parallelizable) to analyze full genomes.

Unlike most tools and methods that assess the effect of genetic variants, which consider functional impact at the level of the whole organism (e.g. whether it is pathogenic or benign), FIRM predicts functional impacts at the molecular level, namely whether the molecular funtion of the affected protein is compromised as a result of the variant, regardless of whether or not this molecular-level impact has a critical downstream effect at the organism level. To achieve this goal, it relies only on biochemical and biophysical proteomic features, without making use of evolutionary data (e.g. conservation metrics). The model relies on features extracted from the high-dimensional proteomic context of the variant (1,109 features in total), which include rich annotations from UniProt and Pfam.
The available classifier has been trained on 37,008 annotated variants from the ClinVar dataset. More details about the model, its training, and the features included in it, can be found in our publication (https://doi.org/10.1093/nar/gkz546). 

Currently, FIRM supports only the human genome.

If you use FIRM in a work contributing to a scientific publication, we ask that you cite our publication:

Nadav Brandes, Nathan Linial, Michal Linial, Quantifying gene selection in cancer through protein functional alteration bias, Nucleic Acids Research, gkz546, https://doi.org/10.1093/nar/gkz546


Usage
==========

Python interface
----
    
    >>> import multiprocessing, geneffect, firm
    >>> n_threads = 4
    >>> thread_pool = multiprocessing.Pool(n_threads)
    >>> geneffect_setup = geneffect.Setup('GRCh38') # Must be a human reference genome
    >>> firm.setup_uniprot_tracks(geneffect_setup)
    >>> firm_classifier = firm.load_classifier(geneffect_setup)
    
    >>> snp = geneffect_setup.variant_interpreter.process_snp('17', 43082434, 'G', 'C')
    >>> snp_gene_effect, = snp.cds_gene_effects
    >>> print(snp_gene_effect)
    P38398:R1443G
    >>> print(firm_classifier.predict_adjusted_proba(snp_gene_effect))
    0.11450486370956418
    
    >>> async_firm_classifier = firm.AsyncClassifier(firm_classifier.predict_adjusted_proba, thread_pool = thread_pool, n_threads = n_threads)
    >>> async_firm_classifier.submit_samples(2000 * [snp_gene_effect], callback = lambda scores: print(len(scores), scores[:10]))
    >>> async_firm_classifier.process_remaining_samples()
    2000 [0.11450486370956418, 0.11450486370956418, 0.11450486370956418, 0.11450486370956418, 0.11450486370956418, 0.11450486370956418, 0.11450486370956418, 0.11450486370956418, 0.11450486370956418, 0.11450486370956418]
    >>> async_firm_classifier.close()    


Command-line interface
----

The command **firm_determine_extended_gene_effects_and_scores** can be used to interpret a list of variants and determine their effect scores using a generalized scheme that also assigns rule-based scores to non-missense variants (the trained machine-learning predictor is only used for the missense variants). For details, see our recent paper: https://doi.org/10.1101/812289. Note that the effect scores produced by the extended tool are reversed to the scores produced by FIRM's standard API (here 0 denotes complete damage, and 1 no effect). For more details, run:

.. code-block:: cshell
    
    firm_determine_extended_gene_effects_and_scores --help
    
Another useful command is **firm_list_all_possible_cds_snps**, allowing to list all possible single-nucleotide variants affecting the coding sequences (CDS) of protein-coding genes, and their effect scores. For more details, run:

.. code-block:: cshell
    
    firm_list_all_possible_cds_snps --help
    
Pre-calculated effect scores for all possible CDS-affecting single-nucleotide variants are available for versions GRCh38 and hg19 of the human reference genome via FTP at ftp://ftp.cs.huji.ac.il/users/nadavb/firm_data/. 

Installation
==========

Dependencies:

* numpy
* pandas
* biopython
* scikit-learn (recommended version: 0.24)
* cython (recommended)
* geneffect (https://github.com/nadavbra/geneffect)


Automatic installation (using the installation script)
----------

    >>> wget https://raw.githubusercontent.com/nadavbra/firm/master/install_firm.sh
    >>> chmod a+x install_firm.sh
    >>> ./install_firm.sh
    
The installation script will also install geneffect and all the other dependencies.


Manual installation
----------

Clone the project and run:

    python setup.py install
    
    
If you haven't installed geneffect before, make sure it is properly configured (see instructions at: https://github.com/nadavbra/geneffect).

After installtion, open config.py in your installation (where the "firm" module has been installed), and go over the instructions there. Specifically, you will need to prepare the directory ~/data/pfam/hmm with all of Pfam's HMM profiles for the human proteome. 
        
    
Replicating the model training
==========

(only for advanced users interested in in-depth understanding of the model)

To replicate the model training, read and follow the instructions in the "Classifier Training with ClinVar Data" Jupyter Notebook.
