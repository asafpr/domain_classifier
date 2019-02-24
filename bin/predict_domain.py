#!/usr/bin/env python
"""
Use PFAM domains to predict a DNA or protein sequence taxonomic domain (Eukaryota, Bacteria, Viruses or Archaea)
Should run build_db.py prior to running this to download files and build diamond db
Implement a Naive-Bayes classifier and return the posterior probability for each domain
"""

from __future__ import division
import domain_classifier
import sys
import os
import os.path
import argparse
from collections import defaultdict
from Bio import SeqIO
import gzip

def process_command_line(argv):
    """
    Return settings object.
    `argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
    """
    if argv is None:
        argv = sys.argv[1:]
    parser = argparse.ArgumentParser(
        description='Script description here',
        formatter_class=argparse.HelpFormatter)
    parser.add_argument("-d", "--dir", default=".",
                        help="Database directory")
    parser.add_argument("-p", "--protein", default=False, action='store_true',
                        help="Input data is protein")
    parser.add_argument("-a", "--all", action='store_true', default=False,
                        help='Predict domain for all the sequences together. By default return a prediction for each fasta record')
    parser.add_argument("-i", "--input",
                        help="input file in fasta format, might have gz extension")
    parser.add_argument("-b", "--blout", 
                        help="diamond output file")
    parser.add_argument("-t", "--translate", default=None,
                        help="Write ORFs to this file")
    parser.add_argument("--limit", type=int, default=100,
                        help="Limit the number of proteins from each sequence to this number. Default is 100, will only apply to hmmsearch, not diamond")
    parser.add_argument("--diamond", default=False, action='store_true',
                        help='Use diamond to find domains, default is hmmsearch')
    parser.add_argument("--threads", type=int, default=20, help="Number of threads to use with hmmsearch")
    parser.add_argument("--pseudocounts", type=int, default=1,
                        help='Add pseudocounts to the number of domains to implement Laplace smoothing. One by default (i.e. Lidstone Smoothing)') 
    parser.add_argument(
        '--hmmsearch', default='hmmsearch',
        help='hmmsearch executable, any other options like --mpi can be set here, default hmmsearch')
    parser.add_argument(
        '--getorf', default='getorf',
        help='getorf EMBOSS executable, default getorf')
    settings = parser.parse_args(argv)
    return settings


def main(argv=None):
    settings = process_command_line(argv)
    # dictionary: sequence -> list(domains)
    if settings.diamond:
        domains = domain_classifier.find_domains(settings.input, settings.dir, settings.blout, settings.all, settings.protein, threads=settings.threads)
    else:
        domains = domain_classifier.find_domains_hmm(
            settings.input, settings.dir, settings.translate, settings.blout, settings.all, settings.protein, 
            hmmsearch=settings.hmmsearch, getorf=settings.getorf, threads=settings.threads, shuffle=settings.limit)
    # likels: dictionary domain -> taxdomain -> likelihood lorder: list(taxdomain)
    (likels, lorder) = domain_classifier.read_likelihoods(settings.dir, settings.pseudocounts)
    # dictionary sequence -> taxdomain -> posterior prob
    posteriors = domain_classifier.compute_post(domains, likels)
    print "\t".join(["Record", "Description", "Length", "Number of Domains", "MAP"] + lorder)
    # Get the sequences lengths
    slens = defaultdict(str)
    desc = defaultdict(str)
    if (settings.input.endswith('gz')):
        fain = gzip.open(settings.input)
    else:
        fain = open(settings.input)
    for sr in SeqIO.parse(fain, 'fasta'):
        slens[sr.id] = len(sr.seq)
        desc[sr.id] = sr.description
    fain.close()
    for k in posteriors.keys():
        maxp = max(posteriors[k], key=posteriors[k].get)
        print "\t".join([str(y) for y in [k, desc[k], slens[k], len(domains[k]), maxp] + [posteriors[k][x] for x in lorder]])
    return 0

if __name__ == '__main__':
    status = main()
    sys.exit(status)


