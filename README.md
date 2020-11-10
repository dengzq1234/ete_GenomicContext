# ete_GenomicContext



Example:  

`ipython GeCo_graphication.py -- --cluster 001_756_949 --notation KEGG`

### MONGO DB config
change mongo config in `mongo.cnf` file to connect correct GMGCdb

## Usage
Program defautly output phylogenetic trees and corresponding neighbour gene annotation(upon your request). For other annotations, simply add the their flags.


    # query KEGG annotation of neighbour gene of GMGC cluster 001_754_949 with biome annotations

    python gmgfam_plot.py --cluster 001_754_949 --notation KEGG --nside 10 --biome

    ## for HPC cluster user
    module load ETE
    module load Python/3.7.2-GCCcore-8.2.0
    module load SciPy-bundle
    
    # to avoid qt issue on the cluster
    QT_QPA_PLATFORM=offscreen xvfb-run python gmgfam_plot.py --cluster 001_754_949 --notation KEGG --nside 10 --biome
