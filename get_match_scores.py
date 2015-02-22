#!/usr/bin/env python3
'''for each query and each mRNA read, compute their a similarity score. This program reports a summary.'''

from __future__ import print_function, division

from Bio import pairwise2, SeqIO
import csv
import time
import glob
import itertools
import os.path
import json

progress_filename = 'progress.json'

def print_scores(queries):
    '''for each query, print the number of matches found for each score. It's like a histogram.'''
    # column widths
    scores_widths = {}
    for query in queries:
        for score, freq in query['score_frequencies'].items():
            scores_widths[score] = max(scores_widths.get(score, 0), len(str(score)), len(str(freq)))

    # header row
    max_name_length = max(len(query['name']) for query in queries)
    print('{0:{1}} {2:<8}'.format('gene', 2+max_name_length, '(len)'), end='') #arcane format syntax for python2.6-compatibility
    for score, width in sorted(scores_widths.items()):
        print('{0:>{1}} '.format(score, width), end='')
    print('')

    # query rows
    for query in queries:
        # print name and length
        print('{0:{1}} {2:<8}'.format(query['name'],
                                      2+max_name_length,
                                      '({0})'.format(len(query['seq']))),
              end='')
        # print match-score frequencies
        for score, width in sorted(scores_widths.items()):
            print('{0:>{1}} '.format(query['score_frequencies'].get(score, '.'),
                                     width),
                  end='')
        print('')
    print('')

def store_progress(queries, num_reads_checked, start_time):
    print('\n{0} reads checked in {1:.0f} seconds ({2})'.format(num_reads_checked,
                                                                  time.time()-start_time,
                                                                  time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())))
    print_scores(queries)
    with open(progress_filename, 'w') as f:
        json.dump(queries, f, indent=0, sort_keys=True)


def local_similarity(s1, s2):
    '''Compute a match score. When a match score sets a new record, print it.'''
    s1, s2 = sorted((s1, s2), key=len)
    score = 0
    read_length = len(s1)
    for start_loc in range(0, len(s2), read_length): 
        results = pairwise2.align.localxs(s1, s2[start_loc:start_loc+2*read_length], -10, -10)
        new_score = int(results[0][2])
        if new_score > score:
            score = new_score
            if score > local_similarity.best_score:
                local_similarity.best_score = score
                print('new best alignment: ')
                print(pairwise2.format_alignment(*results[0]))
    return score
local_similarity.best_score = 0

# load the reads
fasta_files = glob.glob('*.fas')
reads = itertools.chain.from_iterable(SeqIO.parse(filename, 'fasta') for filename in fasta_files)

# load the queries.  Maybe they're already partially done.
if not os.path.isfile(progress_filename):
    with open('query_sequences.csv') as f:
        queries = list(csv.DictReader(f))
    for query in queries:
        query['seq'] = ''.join(c for c in query['seq'].upper() if c in 'ATCG')
        query['score_frequencies'] = {}
else:
    with open(progress_filename) as f:
        queries = json.load(f)
    # json turned all our int keys to strings. fix that.
    for query in queries:
        query['score_frequencies'] = dict((int(k),v) for k,v in query['score_frequencies'].items())
    # skip finished reads
    num_finished = sum(queries[0]['score_frequencies'].values())
    print('skipping {0} reads'.format(num_finished))
    for index, read in zip(range(num_finished), reads):
        pass
    # set best_score
    for query in queries:
        if max(query['score_frequencies'].keys(), default=0) > local_similarity.best_score:
                local_similarity.best_score = max(query['score_frequencies'].keys())
queries.sort(key=lambda query: len(query['seq']))
                
# compute!
start_time = time.time()
for index, read in enumerate(reads):
    # output progress
    if index % int(1e5/sum(len(r['seq']) for r in queries)) == 0:
        store_progress(queries, index, start_time)
        
    read_str = str(read.seq)
    for query in queries:
        h = local_similarity(read_str, query['seq'])
        query['score_frequencies'][h] = 1 + query['score_frequencies'].get(h, 0)

store_progress(queries, index, start_time)
