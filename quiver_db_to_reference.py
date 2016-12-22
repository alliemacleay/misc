import argparse
import sqlite3
import os


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

    def add_column(self, column):
        self.columns.append(column)

    def populate(self, entries):
        self.entries = [SqlEntry(self, entry) for entry in entries]


class SqlEntry(object):
    """ One row from a sql table
    """
    def __init__(self, sql_table, tuple):
        if len(sql_table.columns) != len(tuple):
            raise ValueError('Column lengths in {} do not match length of tuple.  {}!={}'
                             ''.format(sql_table.name, len(sql_table.columns), len(tuple)))
        self.dict = {col: val for col, val in zip(sql_table.columns, tuple)}


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
    sql_table.populate(fusion_details)
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


def make_reference(table, out):
    """
    Write Fusion information to reference file
    If strand is -1 the gene sequence is the reverse complement of corresponding hg19 location,
       start at higher pos and go backwards
    """
    chrom_dict = get_chroms(table)
    breakpoint_dict = {}
    with open(out, 'w') as ref:
        for entry_i in table.entries:
            entry = entry_i.dict
            if entry['FusionOrder'] != 0:  # Take first entry only to find correct breakpoint
                continue
            chrom = entry['Chromosome']
            if chrom is None and entry['GeneName'] is not None:
                chrom = chrom_dict.get(entry['GeneName'], 'chrUnknown')
            chrom = chrom.replace('chr', '')
            fusion_id = entry['FusionID']
            pos = entry['GenomePosStart']
            end = entry['GenomePosEnd']
            breakpoint_pos = abs(int(pos) - int(end)) if end is not None and pos is not None else ""
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


def main(db_file, out_file, bp_file=None):
    """ Main entry point """
    sql_table = get_sql_table(db_file)
    breakpoints = make_reference(sql_table, out_file)
    print 'Reference written to {}'.format(out_file)
    if bp_file is not None:
        write_breakpoints(breakpoints, bp_file)
        print 'Breakpoints for {} written to {}'.format(out_file, bp_file)
    else:
        print 'Breakpoint file NOT written'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Quiver Database to Fusion Reference File",
                                     description="Specify a path to a sqlite db file to create a fusion reference file"
                                                 "\nhttp://archerdx.com/software/quiver")
    parser.add_argument("--db-file", required=True, help="Path to a .db file")
    parser.add_argument("--out-file", required=True, help="Path to output reference file")
    parser.add_argument("--bp-out-file", default=None, help="Path to breakpoints file")
    args = parser.parse_args()
    if not os.path.exists(args.db_file):
        err_txt = 'Database file {} does not exist'.format(args.db_file)
        raise ValueError(err_txt)
    main(args.db_file, args.out_file, args.bp_out_file)