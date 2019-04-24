# domain classifier - Classify a sequence to its taxonomic domain using PFAM domains
## Install
pip install domain-classifier
## How to use
The input is usually a DNA fasta sequence which will be translated to proteins and will be searched for PFAM domains. These domains will be used to classify each sequence to one of four domains: Bacteria, Eukaryota, Archaea and Viruses. The main idea is for the user to be easily spot contaminations where the reference genomes are not known. Since protein domains are highly distinctive between the above taxonomic domain (with the exception of Viruses which kidnap their hosts' genes), the distinction is quite easy. 

First you'll have to download the data from PFAM ftp server:
```
module load diamond
build_domains_DB.py -d <database_dir>
```
If diamond is not installed or you wish not to use diamond for domain prediction (which is fine) you can add `--nodiamond`.

Now take a fasta file and run the classifier, you'll have to define two files for intermediate steps: one for the translated protein sequences and one for the hmmsearch results. Output will be written to STDOUT:
```
# You'll need hmmer and EMBOSS loaded
module load hmmer emboss
predict_domain.py -d <database_dir> -i input.fna -t input.faa -b input.hmm > output_table.txt
```
The script does some funny things like concatenating all the protein sequences and selecting only 100 proteins for each sequence. These are done to save running time of course, you can change the number of proteins using `--limit [int]` flag, I found it unnecessary with the sequences I tested. 
Other options are `--diamond` which will use diamond instead of hmmsearch to find domains, not recommended unless the sequences you have are well known. `--threads` to use another number of threads (default is 20), `--getorf` and `--hmmsearch` to define other paths to EMBOSS getorf and hmmer hmmsearch. 
`--pseudocounts` allows you to introduce more pseudocounts to the Naive-Bayes classifier initial counts (number of genomes the domain was found in) to introduce some uncertainty in the results, the default is 1 (just to avoid log of zero).

## Output
The output table will contain the maximum a-posterior (MAP) which is the most probable domain and the log probability of each of the four domains, unnormalized. You can use some filtering for the minimal number of domains or the difference between the maximal domain and the second best.

## Filtering kraken2 database
The first step would be to run each library.fna file of the downloaded kraken2 database. The next step would be to use the script `filter_kraken_db.py` with the input fna file and the table to get the records that match their taxonomic domain. Run:
```
filter_kraken_db.py --dbfile <taxonomy.sqlite> --taxonomy <kraken_db/taxonomy> <input.fna> <output_table.txt>
```
The script will generate an sqlite3 database in the file given by `--dbfile` if file doesn't exist, if omitted it will be written in memory. The `--taxonomy` flag is the taxonomy directory under the kraken2 database dir, if `--dbfile` exists it can be omitted. Some criteria must be met to filter a sequence:
 - At least `--mindomains` are found for the record (default 5)
 - The record's taxid doesn't match the predicted domain
 - The difference between the MAP domain and the second best is > `--midiff` (1 by default)

If there is no prediction or the prediction didn't meet the criteria the record will be written to the output (which goes to STDOUT)

&copy; 2019 The Jackson Laboratory
