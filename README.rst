What is FIRM?
--------

FIRM (Functional Impact Rating at the Molecular-level) is a machine-learning model for predicting the functional impact of genetic variants (currently focusing only on missense variants affecting protein sequences). FIRM assigns variants with scores between 0 (harmless) to 1 (harmful). The model can be easily installed and run as a standalone Python module, and it is sufficiently fast (and parallelizable) to analyze full genomes.
Unlike most tools and methods that assess the effect of genetic variants, which consider functional impact at the level of the whole organism (e.g. whether it is pathogenic or benign), FIRM predicts functional impacts at the molecular level, namely whether the molecular funtion of the affected protein is compromised as a result of the variant, regardless of whether or not this molecular-level impact has a critical downstream effect at the organism level. To achieve this goal, it relies only on biochemical and biophysical proteomic features, without making use of evolutionary data (e.g. conservation metrics). The model relies on features extracted from the high-dimensional proteomic context of the variant (1,109 features in total), which include rich annotations from UniProt and Pfam.
The available classifier has been trained on 37,008 annotated variants from the ClinVar dataset. 
Currently, FIRM supports only the human genome.


Usage
--------

    >>> from __future__ import print_function
    >>> import multiprocessing, geneffect, firm
    >>> n_threads = 4
    >>> thread_pool = multiprocessing.Pool(n_threads)
    >>> geneffect_setup = geneffect.Setup('GRCh38') # Must be a human reference genome
    >>> firm.setup_uniprot_tracks(geneffect_setup)
    >>> firm_classifier = firm.load_classifier(geneffect_setup)
    
    >>> snp = geneffect_setup.variant_interpreter.process_snp('17', 43082434, 'G', 'C')
    >>> snp_gene_effect, = snp.gene_effects
    >>> print(snp_gene_effect)
    P38398:R1443G
    >>> print(firm_classifier.predict_adjusted_proba(snp_gene_effect))
    0.11450486371
    
    >>> async_firm_classifier = firm.AsyncClassifier(firm_classifier.predict_adjusted_proba, thread_pool = thread_pool, n_threads = n_threads)
    >>> async_firm_classifier.submit_samples(2000 * [snp_gene_effect], callback = lambda scores: print(len(scores), scores[:10]))
    >>> async_firm_classifier.process_remaining_samples()
    2000 [0.11450486370956418, 0.11450486370956418, 0.11450486370956418, 0.11450486370956418, 0.11450486370956418, 0.11450486370956418, 0.11450486370956418, 0.11450486370956418, 0.11450486370956418, 0.11450486370956418]
    >>> async_firm_classifier.close()    


Installation
--------

Dependencies:

* numpy
* pandas
* biopython
* scikit-learn
* geneffect (https://github.com/nadavbra/geneffect)


To install, just run:

    python setup.py install
    
    
After installtion, open config.py in your installation (where the "firm" module has been installed), and go over the instructions there. Specifically, you will need to prepare the directory ~/data/pfam/hmm with all of Pfam's HMM profiles for the human proteome. 
        
    
Replicating the model training
--------

(only for advanced users interested in in-depth understanding of the model)

To replicate the model training, read and follow the instructions in the "Classifier Training with ClinVar Data" Jupyter Notebook.
