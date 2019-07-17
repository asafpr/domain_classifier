#!/usr/bin/env python
"""
Assume you have diamond and wget
"""

from __future__ import division
import sys
import os
import os.path
import argparse
import subprocess
import domain_classifier
from collections import defaultdict

def process_command_line(argv):
    """
    Return settings object.
    `argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
    """
    if argv is None:
        argv = sys.argv[1:]
    parser = argparse.ArgumentParser(
        description='Build PFAM database for classifier',
        formatter_class=argparse.HelpFormatter)
    parser.add_argument(
        '-d', '--dir', default='.', 
        help="database directory")
    parser.add_argument(
        '--nodiamond', default=False, action='store_true',
        help="Don't build diamond DB, predict_domain.py with --diamond will fail")
    settings = parser.parse_args(argv)
    return settings


def main(argv=None):
    settings = process_command_line(argv)
    outdir = os.path.abspath(settings.dir)
    hmm_in = "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
    tax_base = "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/database_files/pfamA_tax_depth.txt.gz"
    alltax = "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/database_files/ncbi_taxonomy.txt.gz" 
    subprocess.check_call("wget -q -P {} {} {} {} ".format(outdir, tax_base, hmm_in, alltax), shell=True)
    # Run diamond build db
    if (not settings.nodiamond):
        fasta_in = "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.fasta.gz"
        subprocess.check_call("wget -q -P {} {}".format(outdir, fasta_in), shell=True)
        subprocess.check_call("diamond makedb --in {}/Pfam-A.fasta.gz -d {}/Pfam-A".format(outdir, outdir), shell=True)
    subprocess.check_call("gunzip {}/Pfam-A.hmm.gz".format(outdir), shell=True)
    # Read the distribution file and select the accessions to write
    with open("{}/domainchisq.txt".format(outdir), "wb") as writef:
        acc = domain_classifier.select_domains("{}/pfamA_tax_depth.txt.gz".format(outdir), writef=writef)
    domain_classifier.filter_hmm("{}/Pfam-A.hmm".format(outdir), "{}/Pfam-A-filter.hmm".format(outdir), acc)
#    subprocess.check_call("hmmpress -f {}/Pfam-A.hmm".format(outdir), shell=True)
    
    return 0

if __name__ == '__main__':
    status = main()
    sys.exit(status)


