# ete_GenomicContext



Example:  

`ipython GeCo_graphication.py -- --cluster 001_756_949 --notation KEGG`

## Usage
Program defautly output phylogenetic trees and corresponding neighbour gene annotation(upon your request). For other annotations, simply add the their flags.


    # query KEGG annotation of neighbour gene of GMGC cluster 001_754_949 with biome annotations

    python gmgfam_ply.py --cluster 001_754_949 --notation KEGG --nside 10 --biome
