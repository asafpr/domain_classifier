from __future__ import division
import subprocess
from collections import defaultdict
import gzip
import tempfile
import sys
import math
import random
import os
import concurrent.futures as cf

def worker(hmmcmd, hmmout, faafile, dbdir, fnum):
#    sys.stderr.write(hmmcmd + " --domtblout {}/hmmout_{}.out {}/Pfam-A.hmm {}/prot_{}.faa\n".format(hmmout, fnum, dbdir, faaout, fnum))
    os.system(hmmcmd + " --domtblout {}/hmmout_{}.out {}/Pfam-A.hmm {}/prot_{}.faa".format(hmmout, fnum, dbdir, faafile, fnum))
#    os.system(r'echo -e "' + fasta + r'" > ' +  "{}/test_{}.txt".format( hmmout, fnum))
def find_domains(infile, dbdir, blout, all=False, protein=False, scov=50, minsim=20, threads=20):
    """
    Find the domains on the input file. For each sequence (or all if all==True) return a list of domains present on it
    Arguments:
    - `infile`: input fasta file or other file that diamond can read
    - `dbdir`: database dir, contains the diamond database
    - `blout`: Diamond output file
    - `all`: Return one list for the entire file
    - `protein`: The input file is protein, use blastp instead of blastx
    - `scov`: Minimal subject coverage
    - `minsim`: Minimal sequence similarity
    - `threads`: number of threads to use
    """
    # Run diamond
    dcmd = "diamond {} -p {} -d {}/Pfam-A -f 6 -q {} --subject-cover {} --id {} |sort -k12gr >  {}".format(("blastp" if protein else "blastx"), threads, dbdir, infile, scov, minsim, blout)
    subprocess.call(dcmd, shell=True)
    # Read the translation of hits to domains
    htod = {}
    fain = gzip.open("{}/Pfam-A.fasta.gz".format(dbdir), 'r')
    for line in fain:
        if line.startswith(">"):
            spl = line.strip().split()
            htod[spl[0][1:]] = spl[2][:7]    
    # Read diamond results. Return all the hits, even if they overlap
    # Save all the matches sequence -> list(domains)
    alldomains = defaultdict(list)
    # for each sequence keep the positions divided by 10 that are covered by a previous domain
    # Since results are sorted, the better matches will be selected
    covered_pos = defaultdict(set)
    with open(blout, 'r') as rb:
        for line in rb:
            spl  = line.strip().split("\t")
            # Test if region is already covered
            mrange = set(range(int(min(int(spl[6]), int(spl[7]))/10), int(max(int(spl[6]), int(spl[7]))/10)+1))
            ncov = len(mrange & covered_pos[spl[0]])
            if ncov < len(mrange)/2:
                sname = spl[0]
                if protein:
                    sname = spl[0].rsplit("_", 1)[0]
                if all:
                    sname = 'all'
                alldomains[sname].append(htod[spl[1]])
                for pos in mrange:
                    covered_pos[spl[0]].add(pos)
    return alldomains

def find_domains_hmm(infile, dbdir, faafile, hmmout, all=False, protein=False, incscore=20, threads=20, shuffle=100, hmmsearch='hmmsearch', getorf='getorf'):
    """
    Use hmmsearch to look for PFAM domains. Accuarte but longer runtime
    Arguments:
    - `infile`: Input fasta file
    - `dbdir`: database dir, assume Pfam-A.hmm in it
    - `faafile`: Write the protein fasta file here
    - `hmmout`: Write hmmsearch output here
    - `all`: Return one list for the entire file
    - `protein`: The input file is protein, use blastp instead of blastx
    - `incscore`: pass to --incdomT parameter
    - `shuffle`: Select random 1000 (or as defined) proteins for each DNA sequence
    - `hmmsearch`: hmmsearch executable
    - `getorf`: getorf executable
    """
    # Translate the DNA to proteins:
    if not protein:
        nfn = infile
        if infile.endswith(".gz"):
            nffile = tempfile.NamedTemporaryFile()
            subprocess.call("zcat {} > {}".format(infile, nffile.name), shell=True)
            nfn = nffile.name
        # Print the proteins into this file, then concatenate them
        tmppt = tempfile.NamedTemporaryFile(mode='w+')
        orfcmd = "{} -table 1 -find 1 -minsize 300 -sequence {} -outseq {}".format(getorf, nfn, tmppt.name)
        sys.stderr.write("Running: {}\n".format(orfcmd))
        subprocess.check_call(orfcmd, shell=True)
        if infile.endswith(".gz"):
            nffile.close()
        # Concat the proteins
#        os.makedirs(faafile, exist_ok=True)
#        os.makedirs(hmmout, exist_ok=True)
        fnum = 0
        with open(faafile, 'w') as ptout:
            sname = None
            sdesc = ''
            seqbuff = []
            allseqs = []
            for line in tmppt:
                if line.startswith(">"):
                    s = line.split()[0].rsplit("_",1)[0][1:]
                    if s != sname:
                        if sname:
                            # Write the previous proteins
                            buffer = "XXX".join(random.sample(seqbuff, min(len(seqbuff), shuffle)))
#                            with open("{}/prot_{}.faa".format(faafile, str(fnum)), 'w') as ptout:
                            ptout.write(">{} {}\n{}\n".format(sname, sdesc, buffer))
#                            allseqs.append(r">{} {}\n{}\n".format(sname, sdesc, buffer))
                            fnum += 1
                        sname = s
                        sdesc = ''
                        try:
                            sdesc = line.strip().split(" ",4)[4]
                        except IndexError:
                            pass
                        seqbuff = []
                else:
                    seqbuff.append(line.strip())
            if sname:
                buffer = "XXX".join(random.sample(seqbuff, min(len(seqbuff), shuffle)))  
 #               with open("{}/prot_{}.faa".format(faafile, str(fnum)), 'w') as ptout:            
                ptout.write(">{} {}\n{}\n".format(sname, sdesc, buffer))
      #      allseqs.append(r">{} {}\n{}\n".format(sname, sdesc, buffer))

        tmppt.close()                
    else:
        faafile = infile
    # Run hmmsearch
    hmmcmd = "{} --cpu {}  -o /dev/null --incT {} -T {} --domtblout {} {}/Pfam-A.hmm {}".format(hmmsearch, threads, incscore, incscore, hmmout, dbdir, faafile)
    sys.stderr.write("Running:{} {}\n".format(fnum, hmmcmd))
    subprocess.run(hmmcmd, shell=True)
# Ruuning the faa files with a pool of threads
#    with cf.ProcessPoolExecutor(threads) as pp:
#        for f in range(fnum+1):
#            sys.stderr.write("{} {} {} {}\n".format(hmmcmd, hmmout, dbdir, f))
#            pp.submit(worker, hmmcmd, hmmout, faafile, dbdir, f)
   
#    subprocess.check_call(hmmcmd, shell=True, stderr=subprocess.STDOUT)
    # Parse the output same way as with diamond
    alldomains = defaultdict(list)
    covered_pos = defaultdict(set)
 
 #   for f in range(fnum+1):   
    with open(hmmout, 'r') as hin:
        for line in hin:
            if line.startswith("#"):
                continue
            spl = line.strip().split()
            mrange = set(range(int(int(spl[19])/10), int(int(spl[20])/10)+1))
            ncov = len(mrange & covered_pos[spl[0]])
            if ncov < len(mrange) / 2:
                sname = spl[0]
                if all:
                    sname = 'all'
                alldomains[sname].append(spl[4].split(".")[0])
                for pos in mrange:
                    covered_pos[spl[0]].add(pos)
    return alldomains
            

def read_likelihoods(dbdir, pseudocounts=1):
    """
    Read the domain distribution in taxonomy and return the likelihood for each taxonomic domain
    Use pseudocounts to avoid zeros
    Arguments:
    - `dbdir`: Database dir, look for the file pfamA_tax_depth.txt.gz
    - `pseudocounts`: Laplace smmothing factor to use
    """
    # Keep the counts in this dict
    rawc = defaultdict(lambda: defaultdict(int))
    taxin = gzip.open("{}/pfamA_tax_depth.txt.gz".format(dbdir))
    lorder = []
    taxsum = defaultdict(int)
    tain = gzip.open("{}/ncbi_taxonomy.txt.gz".format(dbdir))
    for l in tain:
        spl = l.decode().strip().split("\t")
        try:
            king = spl[2].split(';')[0]
            taxsum[king] += 1
        except IndexError:
            pass
    for line in taxin:
        spl = line.decode().strip().split("\t")
        rawc[spl[0]][spl[1]] = int(spl[2])
   #     taxsum[spl[1]] += int(spl[2])
        if spl[1] not in lorder:
            lorder.append(spl[1])
    # Compute likelihood
    likel = defaultdict(lambda: defaultdict(float))
    for dom in rawc.keys():
#        domsum = 0 #sum(rawc[dom].values()) + pseudocounts*len(lorder)
        for tx in lorder:
            likel[dom][tx] = float(rawc[dom][tx] + pseudocounts)/taxsum[tx]
#        domsum = sum(likel[dom].values())
#        for tx in lorder:
#            likel[dom][tx] = likel[dom][tx]/domsum
    return (likel, lorder)

def compute_post(domains, likel):
    """
    compute posterior probabilities using the domains in each sequence and likelihoods
    Arguments:
    - `domains`: A dictionary with domain names in each sequence
    - `likel`: Likelihoods of domains taxonomies
    """
    post = dict()
    for sname in domains.keys():
        mult = defaultdict(lambda: 0)
        for domain in domains[sname]:
            for tx in likel[domain].keys():
                mult[tx] += math.log(likel[domain][tx])
        if (mult):
            post[sname] = mult
    return post


def select_domains(infile, statthr=6.25, writef=None):
    """
    Given a distribution file of domains select the domains that distinguish between the four tax domains
    Using chi-square test between the number of domains observed and the expected using all the domains in the
    database.
    Return a list of domains that passed the test.
    Arguments:
    - `infile`: an input file with counts of domains for each taxonomy
    - `statthr`: The minimal chi-square statist to include
    """
    domcount = defaultdict(lambda: defaultdict(int))
    alldoms = defaultdict(int)
    fin = gzip.open(infile)
    for l in fin:
        (acc, tax, num, d1, d2) = l.strip().split("\t")
        domcount[acc][tax] = int(num)
        alldoms[tax] += int(num)
    
    # Comput the test for each domain
    asum = sum(alldoms.values())
    for k in alldoms:
        alldoms[k] /= asum
    ret_domains = []
    for domain in domcount.keys():
        csq = 0
        dsum = sum(domcount[domain].values())
        for d in alldoms.keys():
            csq += pow(alldoms[d]*dsum - domcount[domain][d],2)/(alldoms[d]*dsum)
        if csq >= statthr:
            ret_domains.append(domain)
        if writef:
            writef.write("{}\t{}\n".format(domain, csq))
    return(ret_domains)


def filter_hmm(hmmfile, outfile, acc):
    """
    Filter hmm file to contain only the acc in the list
    Arguments:
    - `hmmfile`: input hmm file
    - `outfile`: output hmm file
    - `acc`: list of accessions to keep
    """
    with open(hmmfile, 'r') as inf, open(outfile, 'w') as outf:
        hmbuf = []
        keep = False
        for l in inf:
            hmbuf.append(l)
            if (l.startswith('ACC')):
                # ACC   PF10417.9
                (lab, pf) = l.strip().split()
                (a, v) = pf.split(".")
                if a in acc:
                    keep = True
            if (l.strip()=="//"):
                if keep:
                    for l2 in hmbuf:
                        outf.write(l2)
                keep = False
                hmbuf = []

