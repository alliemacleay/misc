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
    sql_table = SqlTable(c, table_name)
    c.execute('SELECT * FROM {0} INNER JOIN {1} on {0}.FusionID=={1}.FusionID;'.format(table_name, joined_table))
    fusion_details = c.fetchall()
    sql_table.add_column('FusionID')
    sql_table.add_column(joined_table)
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
    """ Write Fusion information to reference file """
    chrom_dict = get_chroms(table)
    with open(out, 'w') as ref:
        for entry_i in table.entries:
            entry = entry_i.dict
            chrom = entry['Chromosome']
            if chrom is None and entry['GeneName'] is not None:
                chrom = chrom_dict.get(entry['GeneName'], 'chrUnknown')
            chrom = chrom.replace('chr', '')
            fusion_id = entry['FusionID']
            pos = entry['GenomePosStart']
            sequence = split_sequence_lines(entry['Sequence'])
            if len(sequence) > 0:
                # > chromosome dna:fusion id:reference name:chr:pos:num
                ref.write('>{0} dna:{1}:QuiverFusions:{0}:{2}:1\n{3}'.format(chrom, fusion_id, pos, sequence))


def main(db_file, out_file):
    """ Main entry point """
    sql_table = get_sql_table(db_file)
    make_reference(sql_table, out_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Quiver Database to Fusion Reference File",
                                     description="Specify a path to a sqlite db file to create a fusion reference file"
                                                 "\nhttp://archerdx.com/software/quiver")
    parser.add_argument("--db-file", required=True, help="Path to a .db file")
    parser.add_argument("--out-file", required=True, help="Path to output reference file")
    args = parser.parse_args()
    if not os.path.exists(args.db_file):
        err_txt = 'Database file {} does not exist'.format(args.db_file)
        raise ValueError(err_txt)
    main(args.db_file, args.out_file)