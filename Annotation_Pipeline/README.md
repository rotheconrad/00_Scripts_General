# Annotation Pipeline for TrEMBL, UniProt, and KEGG

This is a workflow for annotating amino acid sequences in fasta format using
the three databases. Includes instructions for downloading and converting the
databases.

Uses Blast+ (Blastp) and Kofamscan.

Kofamscan and hmm databases can be downloaded from the Download section of [this
page](https://www.genome.jp/tools/kofamkoala/).
The publication is [here](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btz859/5631907)


Kofamscan can be easily installed using a [conda environment](https://docs.conda.io/en/latest/miniconda.html):
    conda create -n kofamscan hmmer parallel
    conda activate kofamscan
    conda install -c conda-forge ruby


Details of Blast+ can be found [here](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download).


Blast+ can be easily installed using a [conda environment](https://docs.conda.io/en/latest/miniconda.html):
    conda create -n blastplus
    conda activate blastplus
    conda install bioconda::blast=2.7.1 conda-forge::gnutls conda-forge::nettle


Step 01: Download and parse databases.
    There are various ways to do this. A straighforward approach is simply with wget.
    wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
    wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
    wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
    wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.dat.gz


    python ../Scripts/Parse_UniProtDBs_datFile.py -i uniprot_sprot.dat -o uniprot_sprot.PARSED.dat.tsv
    python ../Scripts/Parse_UniProtDBs_datFile.py -i uniprot_trembl.dat -o  uniprot_trembl.PARSED.dat.tsv
    # Parsing the TrEMBL.dat file can take 3-4 hours.

    python ../Scripts/Check_KEGG_in_UniProt.py -UP uniprot_sprot.PARSED.dat.tsv -KO KEGG_List.tsv -o KEGG_Not_in_SwissProt.tsv