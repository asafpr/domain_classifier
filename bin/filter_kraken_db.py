#!/usr/bin/env python
"""
Given an output file of predict_domain.py and an input kraken2 database fasta file, filter the kraken file
to include entries that match their classification 
"""

from __future__ import division, print_function
import sys
import os
import os.path
import argparse
import logging
import tempfile
import csv
import sqlite3
from collections import defaultdict
from Bio import SeqIO
import urllib
import subprocess
import shutil
import glob
import apsw

py3=sys.version_info >= (3, 0)
def inext(v):  # next value from iterator
    return next(v) if py3 else v.next()

def process_command_line(argv):
    """
    Return settings object.
    `argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
    """
    if argv is None:
        argv = sys.argv[1:]
    parser = argparse.ArgumentParser(
        description='Filter a kraken DB file according to its domain predictions',
        formatter_class=argparse.HelpFormatter)
    parser.add_argument(
        '--dbfile', default = ":memory:",
        help='Database file. If not given the database will be written to memory. If file exists use the data in the file')
    parser.add_argument(
        '--mindomains', type=int, default=5,
        help='Minimal number of domains to consider. Less than that will pass')
    parser.add_argument(
        '--mindiff', type=float, default=1.0,
        help='Minimal log probability difference between the MAP and the second to best')
    parser.add_argument(
        '--taxonomy',
        help="path taxonomy directory of kraken2 DB, should contain names.dmp, nodes.dmp and *accession2taxid files")
    parser.add_argument(
        '--filter_virus', default=False, action='store_true',
        help='By default keep all sequences from viral genomes, use this to filter them out')
    parser.add_argument(
        '--filter_enviro', default=False, action='store_true',
        help='Set to remove environmental samples')
    parser.add_argument(
        'fasta',
        help='Input fasta file. Headers must be accession numbers or in the format: >kraken:taxid|214684|NC_006670.1')
    parser.add_argument(
        'table',
        help='The predicted domains table, output of predict_domain.py')
    settings = parser.parse_args(argv)
    return settings


def build_db(dbfile, taxdir):
    logging.info(taxdir)
    conn = sqlite3.connect(dbfile)
    c = conn.cursor()
    c.execute('''CREATE TABLE nodes(
                 taxid INTEGER PRIMARY KEY,
                 parent INTEGER,
                 rank text,
                     FOREIGN KEY (parent) REFERENCES nodes(taxid))''')
    c.execute('''CREATE TABLE names(
                 taxid INTEGER,
                 name text NOT NULL,
                 class text NOT NULL,
                     FOREIGN KEY (taxid) REFERENCES nodes(taxid))''')
    c.execute('''CREATE TABLE acc2taxid(
                 acc text NOT NULL,
                 taxid INTEGER NOT NULL,
                     FOREIGN KEY (taxid) REFERENCES nodes(taxid))''')
    c.execute('''CREATE INDEX acc_index ON acc2taxid(acc ASC)''')
    names = "{}/names.dmp".format(taxdir)
    nodes = "{}/nodes.dmp".format(taxdir)
    logging.info("Inserting values to taxonomy tables")
    nodestp = []
    with open(nodes, 'rb') as nodesf:
        for line in nodesf:
            spl = line.strip().split("\t|\t")
            nodestp.append((spl[0], spl[1], spl[2]))
    c.executemany("INSERT INTO nodes VALUES (?,?,?)", nodestp)
    namestp = []
    with open(names, 'rb') as namesf:
        for line in namesf:
            spl = line.strip().split("\t|\t")
            namestp.append((spl[0], spl[1], spl[3]))
    c.executemany("INSERT INTO names VALUES (?,?,?)", namestp)
    accs = glob.glob("{}/*accession2taxid".format(taxdir))
    allaccs = []
    for fname in accs:
        with open(fname, 'rb') as fin:
            _ = fin.next()
            for line in fin:
                spl = line.strip().split()
                allaccs.append((spl[0], spl[2]))
    c.executemany("INSERT INTO acc2taxid VALUES (?,?)",allaccs)
    conn.commit()
    c.execute("select count(*) from nodes")
    nnodes = c.fetchone()[0]
    if (nnodes != len(nodestp)):
        logging.warn("Number of nodes in database is different than in input file ({} vs {})".format(nnodes, len(nodestp)))
    c.execute("select count(*) from names")
    nnames = c.fetchone()[0]
    if (nnames != len(namestp)):
        logging.warn("Numer of names in database is different than in input file ({} vs {})".format(nnames, len(namestp)))
    conn.close()


def match_txid(domain, taxid, curr, keepvir=False, keepenv=False):
    """
    Return True if the txid is under the domain
    """
    while (taxid != domain):
        try:
            (p,) = curr.execute("select parent from nodes where taxid=?", (taxid, )).next()
        except:
            logging.warn("Can't find parent for taxonomy {}".format(taxid))
            break
        if p == taxid:
            break
        if (taxid == 10239 and not keepvir):
            return True
        if (taxid == 48479 and not keepenv):
            return False
        taxid = p

    return taxid == domain 


def main(argv=None):
    settings = process_command_line(argv)
    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)
    if (settings.dbfile == ':memory:' or not os.path.exists(settings.dbfile)):
        build_db(settings.dbfile, settings.taxonomy)
    #dbconn = sqlite3.connect(settings.dbfile)
    diskconn = apsw.Connection(settings.dbfile)
    dbconn = apsw.Connection(":memory:")
    with dbconn.backup("main", diskconn, "main") as backup:
        backup.step()

    curr = dbconn.cursor()
    # Read the results file
    accdom = {}
    domains = {'Eukaryota': 2759, 'Viruses': 10239, 'Bacteria': 2, 'Archaea': 2157}
    with open(settings.table, 'rb') as tbin:
        for row in csv.DictReader(tbin, delimiter="\t"):
            if int(row['Number of Domains']) < settings.mindomains: continue
            best_score = float(row[row['MAP']])
            if best_score - max([float(row[x]) for x in (set(domains.keys())-set((row['MAP'],)))]) < settings.mindiff: continue
            accdom[row['Record']] = domains[row['MAP']]
    # Read the fasta and print the records passed filter to STDOUT
    for record in SeqIO.parse(settings.fasta, 'fasta'):
        if record.id.startswith("kraken"):
            (_, txid, acc) = record.id.split("|",2)
        else:
            acc = record.id
            nover = acc.split(".")[0]
            try:
                (txid,) = curr.execute("select taxid from acc2taxid where acc=?", (nover,)).next()
            except:
                taid = None
#            txid = curr.fetchone()
#            if txid: txid = txid[0]
        if txid and acc in accdom:
            if not match_txid(accdom[acc], txid, curr, settings.filter_virus, not settings.filter_enviro): continue
        if not txid:
            logging.warn("Can't find taxonomy ID for sequence: {}".format(record.id))
        sys.stdout.write(record.format('fasta'))
          
    return 0

if __name__ == '__main__':
    status = main()
    sys.exit(status)


