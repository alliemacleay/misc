import argparse
import pandas as pd

def get_transcript_df(bp_transcripts):
    gene_ids = []
    transcript_ids = []
    fusion_ids = []
    assay_types = []
    for index, row in bp_transcripts.iterrows():
        gene, trn, fus, assay = '', '', '', ''
        for data in row.annotation.split(';')[:-1]:
            col, val = data.lstrip().split(' ')
            val = val.replace('"', '')
            if col == 'gene_id':
                gene = val
            elif col == 'transcript_id':
                trn = val
            elif col == 'fusion_id':
                fus = val
            elif col == 'assay_type':
                assay = val
            else:
                print 'no column for {}'.format(col)
        gene_ids.append(gene)
        transcript_ids.append(trn)
        fusion_ids.append(fus)
        assay_types.append(assay)
    bp_transcripts['Gene_Id'] = gene_ids
    bp_transcripts['Transcript_Id'] = transcript_ids
    bp_transcripts['Fusion_Id'] = fusion_ids
    bp_transcripts['Assay_Type'] = assay_types
    return bp_transcripts

def get_breakpoint_dict(bp_trs):
    breaks = {}
    for id in bp_trs['Fusion_Id'].unique():
        print id
        fus_partners = bp_trs[bp_trs['Fusion_Id']==id]
        for i, ix in enumerate(fus_partners.index):
            if i == 0:
                breaks[id] = []
            if i == fus_partners.ndim - 1:
                # ignore last location
                continue
            breaks[id].append(fus_partners.ix[ix]['end'])
    return breaks

def print_bp_output(bp_dict, map, out):
    with open(out, 'wb') as bp_o:
        for key in bp_dict:
            abbr = map.get(key, '')
            if abbr == '':
                print 'No value for {}'.format(key)
            bp_o.write('{}\t{}\t{}\n'.format(key, ','.join([str(x) for x in bp_dict[key]]), abbr))

def get_ref_map(ref_file):
    map = {}
    with open(ref_file, 'rb') as iref:
        for line in iref:
            if line.startswith('>'):
                data = line.lstrip('>').split('"')
                map[data[1]] = data[0].split(' ')[0]
            else:
                continue
    return map


def main(bp_file, out, ref_file):
    bp_in = pd.read_csv(bp_file, sep='\t', names=['id', 'version', 'type', 'start', 'end', 'order', 'strand', 'transcript_end', 'annotation'])
    bp_trs = get_transcript_df(bp_in[bp_in['type'] == 'FusionSeq'])
    ref_map = get_ref_map(ref_file)
    print_bp_output(get_breakpoint_dict(bp_trs), ref_map, out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Quiver Database Breakpoints",
                                     description="Specify a path to a text file containing quiver breakpoints"
                                                 "\nhttp://archerdx.com/software/quiver")
    parser.add_argument("--out", required=True, help="Path to a output file")
    parser.add_argument("--breakpoint-file", required=True, help="Path to breakpoints file")
    parser.add_argument("--fusion-ref", required=True, help="Path to reference file")
    args = parser.parse_args()
    main(args.breakpoint_file, args.out, args.fusion_ref)