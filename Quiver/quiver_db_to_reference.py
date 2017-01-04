import argparse
import sqlite3
import os
import requests
import sys
import pandas as pd
import time


class SqlTable(object):
    """ SqlTable object
    takes a sql table name and populates the column names
    entries are added with the order of column names
    use add_column to add joined query data
    """
    def __init__(self, cursor, table_name):
        self.name = table_name
        cursor.execute('PRAGMA table_info({});'.format(table_name))
        self.table_info = cursor.fetchall()
        self.columns = [item[1] for item in self.table_info]
        self.entries = []
        self.index = {}

    def add_column(self, column):
        self.columns.append(column)

    def populate(self, entries, index=None):
        self.entries = [SqlEntry(self, entry) for entry in entries]
        if index is not None:
            self.index_table(index)

    def index_table(self, index):
        for i, entry_i in enumerate(self.entries):
            entry = entry_i.dict
            if entry[index] not in self.index:
                self.index[entry[index]] = [i]
            else:
                self.index[entry[index]].append(i)


class SqlEntry(object):
    """ One row from a sql table
    """
    def __init__(self, sql_table, tuple):
        if len(sql_table.columns) != len(tuple):
            raise ValueError('Column lengths in {} do not match length of tuple.  {}!={}'
                             ''.format(sql_table.name, len(sql_table.columns), len(tuple)))
        self.dict = {col: val for col, val in zip(sql_table.columns, tuple)}


class Fusion(object):
    def __init__(self, entry_array, transcript_dict=None, chrom_dict=None):
        self.middle = None
        self.right = None
        self.fusion_id = entry_array[0]['FusionID']
        self.valid = False
        max_inserted = 0
        for e in entry_array:
            if e['InsertedBases'] is not None:
                max_inserted = len(e['InsertedBases']) if len(e['InsertedBases']) > max_inserted else max_inserted
            if e['FusionOrder'] == 0:
                self.left = e
            elif e['FusionOrder'] == 2:
                if self.right is not None:
                    self.middle = self.right
                self.right = e
            elif e['FusionOrder'] == 1:
                self.right = e
            else:
                print 'Order out of range: {} - {}'.format(e['FusionOrder'], e['FusionID'])
        self.left_direction = self.left['Strand']
        self.right_direction = self.right['Strand']
        left_bp = self.left['GenomePosEnd'] if self.left_direction == 1 else self.left['GenomePosStart']
        right_bp = self.right['GenomePosEnd'] if self.right_direction == 1 else self.right['GenomePosStart']
        if self.left['GenomePosEnd'] is not None and self.left['GenomePosStart'] is not None:
            self.seq_bp = abs(int(self.left['GenomePosEnd']) - int(self.left['GenomePosStart']))
        else:
            self.seq_bp = -1
        self.left_gene = self.left['GeneName']
        self.right_gene = self.right['GeneName']
        if self.left['Chromosome'] is None or self.right['Chromosome'] is None:
            if chrom_dict is not None:
                self.left_bp = '{}:{}'.format(chrom_dict[self.left_gene].replace('chr', ''), left_bp)
                self.right_bp = '{}:{}'.format(chrom_dict[self.right_gene].replace('chr', ''), right_bp)
            else:
                return
        else:
            self.left_bp = '{}:{}'.format(self.left['Chromosome'].replace('chr', ''), left_bp)
            self.right_bp = '{}:{}'.format(self.right['Chromosome'].replace('chr', ''), right_bp)
        self.left_exon = self.left['ExonFusion3prime'] if self.left_direction == 1 else self.left['ExonFusion5prime']
        self.right_exon = self.right['ExonFusion5prime'] if self.right_direction == 1 else self.right['ExonFusion3prime']
        transcripts = [e['TranscriptID'] for e in entry_array]
        if transcript_dict is None:
            refseq_ids = [get_ref_seq_id(trs) for trs in transcripts]
        else:
            refseq_ids = [transcript_dict[trs] for trs in transcripts]
        genes = entry_array[0]['GenesFused'].split(':')
        exons = [e['ExonFusion5prime'] for e in entry_array]
        self.fusion_suggested = ':'.join(['{}-{}-Exon{}'.format(rs, g, ex) for rs, g, ex in zip(refseq_ids, genes, exons)])
        self.fusion_type = 'Intragenic' if genes.count(genes[0]) == len(genes) else 'Intergenic'
        left_bp_type = self.left['BreakpointType5prime'] if self.left_direction == 1 else self.left['BreakpointType3prime']
        right_bp_type = self.right['BreakpointType3prime'] if self.right_direction == 1 else self.right['BreakpointType5prime']
        if left_bp_type is None or right_bp_type is None:
            self.fusion_status = ''
        else:
            self.fusion_status = '{}:{}'.format(left_bp_type.replace(' Unknown', ''), right_bp_type.replace(' Unknown', ''))  # ex: Exon:Exon
        self.sequence = entry_array[0]['Sequence']
        self.n_inserted = max_inserted
        self.valid = True

    def detango_string(self):
        """ return string for detango output """
        if not self.valid:
            return ''
        rlim = 0 if self.seq_bp - 50 < 0 else self.seq_bp - 50
        llim = self.seq_bp + self.n_inserted + 50
        llim = len(self.sequence) if llim > len(self.sequence) else llim
        # 'Fusion_Isoform' 'Fusion_Suggested' 'Fusion_Coverage' 'Isoform_fusionStatus' 'Isoform_fusionType'
        detango = [self.fusion_suggested, self.fusion_suggested, '', self.fusion_status, self.fusion_type]
        # 'Left_BreakPoint' 'Right_BreakPoint' 'minBaseOff_ExonBoundary' 'BaseOff_ExonBoundary' 'Frame' 'Frame_Detail
        detango.extend([self.left_bp, self.right_bp, '', '', '', ''])
        # 'Left_Contig' 'Left_Amino' 'N_Insertion' 'Right_Contig' 'Right_Amino' 'Supporting_Read'
        detango.extend([self.sequence[rlim: self.seq_bp], '', self.n_inserted,
                       self.sequence[self.seq_bp + self.n_inserted: llim], '', ''])
        detango = [str(detab).replace('\t', '-') for detab in detango]
        return '\t'.join(detango)


def get_ref_seq_id(ens_transcript_id):
    """ Get RefSeq ID from ensembl rest services """
    server = "http://rest.ensembl.org"
    ext = "/xrefs/id/{}?".format(ens_transcript_id)
    try:
        r = requests.get(server+ext, headers={"Content-Type": "application/json"}, timeout=5)
    except requests.exceptions.ReadTimeout as timeout:
        print 'Timed out\n{}'.format(timeout.message)
        return ens_transcript_id
    requests_left = r.headers['X-RateLimit-Remaining']
    #print(requests_left)
    if int(requests_left) == 0:
        sleep_time = int(r.headers['X-RateLimit-Reset'])
        print 'Sleeping for {} seconds'.format(sleep_time)
        time.sleep(sleep_time)
    if not r.ok:
        return ens_transcript_id
    data = r.json()
    dfr = pd.DataFrame(data)
    dval = dfr[dfr['db_display_name'] == 'RefSeq mRNA']['primary_id'].values
    if len(dval) > 0:
        return dval[0]
    else:
        return ens_transcript_id


def get_sql_table(db_file):
    """ Extract info from Quiver database file"""
    conn = sqlite3.connect(db_file)
    c = conn.cursor()
    table_name = 'Fusion_Arm_Details'
    joined_table = 'Sequence'
    joined_table_2 = 'Known_Fusions_content'
    sql_table = SqlTable(c, table_name)
    #c.execute('SELECT * FROM Fusion_Arm_Details INNER JOIN Sequence on Fusion_Arm_Details.FusionID==Sequence.FusionID INNER JOIN (SELECT docid, c0FusionID, c1ApprovedGeneOrder FROM Known_Fusions_content) ON Fusion_Arm_Details.FusionID==c0FusionID')
    c.execute('SELECT * FROM {0} INNER JOIN {1} on {0}.FusionID=={1}.FusionID INNER JOIN '
              '(SELECT docid, c0FusionID, c1ApprovedGeneOrder FROM {2}) ON {0}.FusionID==c0FusionID;'
              ''.format(table_name, joined_table, joined_table_2))
    fusion_details = c.fetchall()
    sql_table.add_column('FusionID')
    sql_table.add_column(joined_table)
    sql_table.add_column('docid')
    sql_table.add_column('FusionID')
    sql_table.add_column('GenesFused')
    sql_table.populate(fusion_details, 'FusionID')
    conn.close()
    return sql_table


def get_chroms(table):
    """ Some chrom entries are blank in Archer DB so create a default chrom dict"""
    chroms = {}
    for entry_i in table.entries:
        entry = entry_i.dict
        if entry['GeneName'] is not None and entry['Chromosome'] is not None:
            chroms[entry['GeneName']] = entry['Chromosome']
    return chroms


def split_sequence_lines(seq, lim=80):
    """ Certain number of characters per line in sequence section """
    split_seq = ''
    while len(seq) > 0:
        split_seq = '{}{}\n'.format(split_seq, seq[:lim])
        seq = seq[lim:]
    return split_seq


def detango_annotation(entries, transcript_dict, chroms=None):
    fusion_obj = Fusion(entries, transcript_dict, chroms)
    return fusion_obj.detango_string()


def make_reference(table, out):
    """
    Write Fusion information to reference file
    If strand is -1 the gene sequence is the reverse complement of corresponding hg19 location,
       start at higher pos and go backwards
    """
    chrom_dict = get_chroms(table)
    breakpoint_dict = {}
    count = 0
    with open(out, 'w') as ref:
        for entry_i in table.entries:
            count += 1
            if count % 500 == 0:
                print 'Processed {}'.format(count)
            entry = entry_i.dict
            if entry['FusionOrder'] != 0:  # Take first entry only to find correct breakpoint
                continue
            chrom = entry['Chromosome']
            if chrom is None and entry['GeneName'] is not None:
                chrom = chrom_dict.get(entry['GeneName'], 'chrUnknown')
            chrom = chrom.replace('chr', '')
            fusion_id = entry['FusionID']
            transcript_start = entry['TranscriptPosStart']
            transcript_end = entry['TranscriptPosEnd']
            pos = entry['GenomePosStart']
            breakpoint_pos = abs(int(transcript_end) - int(transcript_start)) if \
                transcript_end is not None and transcript_start is not None else ""
            fusion_abbr = '{}_{}'.format(entry['docid'], entry['GenesFused'])
            sequence = split_sequence_lines(entry['Sequence'])
            if len(sequence) > 0:
                #  > chromosome dna:fusion id:reference name:chr:pos:num
                #  >1_CCDC6:RET dna:"CCDC6{ENST00000263102}:r.1_535_RET{ENST00000355710}:r.2369_5659":QuiverFusions:10:10:61665880/535:1
                ref.write('>{4} dna:"{1}":QuiverFusions:{0}:{0}:{2}/{5}:1\n{3}'
                          ''.format(chrom, fusion_id, pos, sequence, fusion_abbr, breakpoint_pos))
                breakpoint_dict[fusion_abbr] = breakpoint_pos
    return breakpoint_dict


def write_breakpoints(bp, bp_file):
    with open(bp_file, 'w') as bpf:
        for fusion in bp:
            bpf.write('{}\t{}\n'.format(fusion, bp[fusion]))


def write_detango_annotation(dt, dt_file):
    detango_flds = ['Fusion_Isoform',
                    'Fusion_Suggested',
                    'Fusion_Coverage',
                    'Isoform_fusionStatus',
                    'Isoform_fusionType',
                    'Left_BreakPoint',
                    'Right_BreakPoint',
                    'minBaseOff_ExonBoundary',
                    'BaseOff_ExonBoundary',
                    'Frame',
                    'Frame_Detail',
                    'Left_Contig',
                    'Left_Amino',
                    'N_Insertion',
                    'Right_Contig',
                    'Right_Amino',
                    'Supporting_Read',
                    ]
    with open(dt_file, 'w') as dtf:
        dtf.write('Fusion_Abbr\t{}\n'.format('\t'.join(detango_flds)))
        for fusion in dt:
            dtf.write('{}\n'.format(dt[fusion]))


def make_detango_dict(table, transcript_dict):
    """ Make dictionary with detango fields """
    size = len(table.index.keys())
    count = 0
    detango_dict = {}
    chrom = get_chroms(table)
    for fusion_id in table.index:
        count += 1
        if count % 500 == 0:
            print '{} of {} detango annotations processed'.format(count, size)
        fus_id_entries = [table.entries[x].dict for x in table.index[fusion_id]]
        if fus_id_entries[0]['Sequence'] == '':
            continue
        fusion_abbr = '{}_{}'.format(fus_id_entries[0]['docid'], fus_id_entries[0]['GenesFused'])
        detango_dict[fusion_id] = '{}\t{}'.format(fusion_abbr, detango_annotation(fus_id_entries, transcript_dict, chrom))
        if detango_dict[fusion_id] == '':
            del detango_dict[fusion_id]
    return detango_dict


def compile_new_transcript_file(table, out):
    trs_dict = {}
    size = len(table.entries)
    count = 0
    for entry_i in table.entries:
        count += 1
        if count % 500 == 0:
            print '{} of {} transcripts found'.format(count, size)
        entry = entry_i.dict
        t_id = entry['TranscriptID']
        if t_id is None:
            continue
        if t_id not in trs_dict:
            trs_dict[t_id] = get_ref_seq_id(t_id) if not t_id.startswith('NM') else t_id
    return trs_dict


def write_transcript_file(tdict, out):
    with open(out, 'wb') as tfile:
        for t_id in tdict:
            tfile.write('{}\t{}\n'.format(t_id, tdict[t_id]))


def read_transcript_file(infile):
    tdict = {}
    with open(infile, 'rb') as ifile:
        for line in ifile:
            t_id, t_ref_seq = line.strip().split('\t')
            tdict[t_id] = t_ref_seq
    return tdict


def main(db_file, out_file, bp_file=None, dt_file=None, trs_out=None, trs_in=None):
    """ Main entry point """
    sql_table = get_sql_table(db_file)
    breakpoints = make_reference(sql_table, out_file)
    trs_dict = None
    print 'Reference written to {}'.format(out_file)
    if bp_file is not None:
        print 'Writing breakpoint file...'
        write_breakpoints(breakpoints, bp_file)
        print 'Breakpoints for {} written to {}'.format(out_file, bp_file)
    else:
        print 'Breakpoint file NOT written'
    if trs_out is not None:
        print 'Compiling new transcript file.  This will take a while...'
        trs_dict = compile_new_transcript_file(sql_table, trs_out)
        write_transcript_file(trs_dict, trs_out)
        print 'Transcript file written to {}'.format(trs_out)
    if trs_dict is None and trs_in is not None:
        trs_dict = read_transcript_file(trs_in)
    if trs_dict is None:
        print 'Specify --transcript-out-file OR --transcript-in-file to create or use a transcript dictionary to write' \
              ' detango annotation file'

    if dt_file is not None and trs_dict is not None:
        print 'Compiling DeTango annotation file...'
        detango = make_detango_dict(sql_table, trs_dict)
        print 'Writing DeTango annotation file...'
        write_detango_annotation(detango, dt_file)
        print 'Detango annotation for {} written to {}'.format(out_file, dt_file)
    else:
        print 'Detango annotation file NOT written'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Quiver Database to Fusion Reference File",
                                     description="Specify a path to a sqlite db file to create a fusion reference file"
                                                 "\nhttp://archerdx.com/software/quiver")
    parser.add_argument("--db-file", required=True, help="Path to a .db file")
    parser.add_argument("--out-file", required=True, help="Path to output reference file")
    parser.add_argument("--bp-out-file", default=None, help="Path to breakpoints file")
    parser.add_argument("--dt-out-file", default=None, help="Path to detango annotation file")
    parser.add_argument("--transcript-out-file", default=None, help="Path to transcript annotation file")
    parser.add_argument("--transcript-in-file", default=None, help="Path to transcript annotation file to use")
    args = parser.parse_args()
    if not os.path.exists(args.db_file):
        err_txt = 'Database file {} does not exist'.format(args.db_file)
        raise ValueError(err_txt)
    main(args.db_file, args.out_file, args.bp_out_file, args.dt_out_file, args.transcript_out_file,  args.transcript_in_file)