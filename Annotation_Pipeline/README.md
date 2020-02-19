# Annotation Pipeline for TrEMBL, UniProt, and KEGG

This is a workflow for annotating amino acid sequences in fasta format using
the three databases. Includes instructions for downloading and converting the
databases.

## Step 0: Required tools - Blast+ (Blastp) and Kofamscan.

Kofamscan and hmm databases can be downloaded from the Download section of [this
page](https://www.genome.jp/tools/kofamkoala/).
The publication is [here](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btz859/5631907)


Kofamscan can be easily installed using a [conda environment](https://docs.conda.io/en/latest/miniconda.html):
    ```bash
    conda create -n kofamscan hmmer parallel
    conda activate kofamscan
    conda install -c conda-forge ruby
    ```

Kofamscan can be run with default settings like this:
    ```bash
    ruby -o {outfile_name} {input_fasta}
    ```

Details of Blast+ can be found [here](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download).


Blast+ can be easily installed using a [conda environment](https://docs.conda.io/en/latest/miniconda.html):
    ```bash
    conda create -n blastplus
    conda activate blastplus
    conda install bioconda::blast=2.7.1 conda-forge::gnutls conda-forge::nettle
    ```

Blastp can be run with similar to this (the -outfmt 6 format order is required for downstream processing):
    ```bash
    blastp -task 'blastp' -evalue 0.01 -max_target_seqs 10 -num_threads 2 -db {pathto_db} -query {input_fasta} -out {outfile_name} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'
    ```

## Step 01: Download and parse databases.

There are various ways to download the databases. A straighforward approach is simply with wget:
    ```bash
    wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
    wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
    wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
    wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.dat.gz
    wget 'https://www.genome.jp/kegg-bin/download_htext?htext=ko00001&format=htext&filedir=' -O ko00001.keg
    ```

I wrote some python code to parse these database files for use downstream.
    ```python
    * python ../Scripts/Parse_UniProtDBs_datFile.py -i uniprot_sprot.dat -o uniprot_sprot.PARSED.dat.tsv
    * python ../Scripts/Parse_UniProtDBs_datFile.py -i uniprot_trembl.dat -o  uniprot_trembl.PARSED.dat.tsv
    * python ../Scripts/Check_KEGG_in_UniProt.py -UP uniprot_sprot.PARSED.dat.tsv -KO KEGG_List.tsv -o KEGG_Not_in_SwissProt.tsv
    ```
    * Parsing the TrEMBL.dat file can take 3-4 hours.