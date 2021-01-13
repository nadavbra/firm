import os

import pandas as pd

from firm.util import log

def create_gene_record(gene):
    
    cds_coordinates = [coordinate for cds_exon in gene.canonical_cds_isoform.cds_exons for coordinate in \
            [cds_exon.chromosome_start, cds_exon.chromosome_end]]

    return [
        gene.uniprot_record.id,
        gene.symbol,
        gene.name,
        list(gene.refseq_ids),
        gene.canonical_cds_isoform.chromosome,
        min(cds_coordinates),
        max(cds_coordinates),
    ]
    
def create_gene_dataset(geneffect_setup):
    return pd.DataFrame(list(map(create_gene_record, geneffect_setup.genes)), columns = ['uniprot_id', 'symbol', 'name', 'refseq_ids', \
            'chr', 'cds_start', 'cds_end'])
            
def get_or_create_gene_dataset(genes_dir, geneffect_setup):
    
    gene_dataset = create_gene_dataset(geneffect_setup)
    gene_dataset_csv_file_path = os.path.join(genes_dir, 'genes_%s.csv' % geneffect_setup._config_setup.ref_genome)
    
    if os.path.exists(gene_dataset_csv_file_path):
        assert (gene_dataset['uniprot_id'] == pd.read_csv(gene_dataset_csv_file_path, index_col = 0)['uniprot_id']).all()
        log('The %d genes in the created gene dataset are identical to the ones in the existing CSV file: %s' % (len(gene_dataset), \
                gene_dataset_csv_file_path))
    else:
        gene_dataset.to_csv(gene_dataset_csv_file_path)
        log('Saved the gene dataset (%d genes) into the CSV file: %s' % (len(gene_dataset), gene_dataset_csv_file_path))
        
    return gene_dataset
