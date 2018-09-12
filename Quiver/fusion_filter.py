import pysam
import argparse
import pandas as pd
import numpy as np

CIGAR_CODES = {0: 'M',
               1: 'I',
               2: 'D',
               3: '',
               4: 'S',
               5: 'H'}


def match_spans_bp(bp, pos, cigar):
    """
    Does the sequence match on both sides of the breakpoint?
    :param bp: breakpoint
    :param cigar: read.cigartuple
    :return: boolean
    """
    for cig in cigar:
        if pos > bp:
            break
        if CIGAR_CODES[cig[0]] == 'M':
            if pos < bp < (pos + cig[1]):
                return True
        if CIGAR_CODES[cig[0]] not in ['H', 'S']:
            pos += cig[1]
    return False


def match_spans(bps, pos, cigar):
    """ Make sure all breakpoints in list bps are covered"""
    for bp in bps:
        if bp == '':
            continue
        if match_spans_bp(int(bp), pos, cigar):
            return True
    return False


def get_breakpoint_dict(bpf, fusd):
    """ Create breakpoint dictionary from file """
    bp_dict = {}
    with open(bpf, 'rb') as bp:
        for line in bp:
            tdat = line.strip().split('\t')
            fusion, breakpoint = tdat[0], tdat[1]
            breakpoint = breakpoint.strip()
            fus_abbr = fusd.get(fusion, '')
            if fus_abbr == '':
                print '{} not found'.format(fusion)
                breakpoint = ''
            if not breakpoint == '':
                bp_dict[fus_abbr] = breakpoint
            else:
                print 'No valid breakpoint for {}'.format(fusion)
    return bp_dict

def get_breakpoint_dict_meta(meta_df):
    """ Create breakpoint dictionary from file """
    bp_info = meta_df[meta_df['type'] == 'FusionSeq']
    bps = []
    for id in set(bp_info['id']):
        fsid = bp_info[bp_info['id'] == id]
        lfs = len(fsid)
        bps.append(','.join([str(bp) for bp in list(fsid['end'])[:lfs - 1]]))
    return dict(zip(set(bp_info['id']), bps))

def get_fusion_dict_meta(meta_df):
    fus_slice = meta_df[meta_df['start'] == 1][meta_df['type'] == 'FusionSeq']
    fusion_ids = [ann.split('fusion_id')[1].split(';')[0].lstrip() for ann in list(fus_slice['annotation'])]
    record_ids = list(fus_slice['id'])
    return dict(zip(record_ids, fusion_ids))


def write_bam(sam, read_dict, out):
    """ Write filtered bam to output file """
    cov_dict = {}
    sam_out = pysam.AlignmentFile(out, 'wb', template=sam)
    for fusion in read_dict:
        cov_dict[fusion] = len(read_dict[fusion])
        for read in read_dict[fusion]:
            sam_out.write(read)
    sam_out.close()
    return cov_dict


def covered(bps, pos, end, distance):
    """ Check each breakpoint in bps list for coverage at a certain distance """
    for bp in bps:
        if covered_bps(int(bp), pos, end, distance):
            return True
    return False


def covered_bps(bp, pos, end, distance):
    """ Check for coverage at N (distance) bases before and after breakpoint """
    if pos < bp - distance and end > bp + distance:
        return True
    return False


def detango_output(cov, outfile, annotation, fusion_ids):
    records = []
    ann = pd.read_csv(annotation, sep='\t')
    df = pd.DataFrame.from_dict(data=cov, orient='index')
    for ix in df.index:
        fusion_abbr = fusion_ids[ix].replace('"', '')
        rec = ann[ann['Fusion_Abbr'] == fusion_abbr]
        if len(rec) == 0:
            print '{} not found.'.format(fusion_abbr)
            continue
        records.append(rec.index[0])
    ndf = ann.iloc[records].to_dict()
    mkeys = ndf['Fusion_Suggested'].keys()
    cov_dict = {k: v for k, v in zip(mkeys, df[0].tolist())}
    ndf['Fusion_Coverage'] = cov_dict
    filtered = pd.DataFrame.from_dict(ndf)
    detango_fields = ['Fusion_Isoform',
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
                      'Supporting_Read']

    filtered.replace(np.nan, "", regex=True).to_csv(outfile, sep='\t', index=False, columns=detango_fields)


def analyze_coverage(cov, outfile):
    df = pd.DataFrame.from_dict(data=cov, orient='index')
    df.to_csv(outfile)

def get_fusion_dict(fus_id):
    fus_dict = {}
    with open(fus_id, 'rb') as fh:
        for line in fh:
            id, abbr = line.strip().split('\t')
            fus_dict[id] = abbr
    return fus_dict


def process_sam(sam, breakpoint_dict, min_coverage, distance):
    count = 0
    reads_spanning = {}
    distances = [3, 5, 10, 25, 35]
    cov = {k: 0 for k in distances}
    required_coverage = {key: 0 for key in breakpoint_dict}
    coverage_analysis = {key: cov.copy() for key in breakpoint_dict}
    tot = len(breakpoint_dict)
    txt_found = ''
    for fusion in breakpoint_dict:
        count += 1
        if len(reads_spanning) > 0:
            txt_found = ' {} valid fusions found'.format(len(reads_spanning))
        print 'Parsing reads for {}.  ({}/{}){}'.format(fusion, count, tot, txt_found)
        reads_spanning[fusion] = []
        breakpoint = breakpoint_dict[fusion]
        for read in sam.fetch(str(fusion), 1, 1000000):
            if read.cigartuples is not None and match_spans(breakpoint.split(','), read.pos, read.cigartuples):
                #  Aggregate coverage data at different distances for visualization
                for d in distances:
                    if covered(breakpoint.split(','), read.pos, read.aend, d):
                        coverage_analysis[fusion][d] += 1
                #  Check distance parameter coverage
                if covered(breakpoint.split(','), read.pos, read.aend, distance):
                    required_coverage[fusion] += 1
                reads_spanning[fusion].append(read)
        if len(reads_spanning[fusion]) == 0:
            del reads_spanning[fusion]
        if coverage_analysis[fusion][distances[0]] == 0:
            del coverage_analysis[fusion]
        if required_coverage[fusion] < min_coverage:
            del required_coverage[fusion]
    return reads_spanning, coverage_analysis, required_coverage


def main(bam, meta, out, det_ann, min_coverage, distance, fusion_out, analysis_out):
    """ Main entry point """
    meta_df = pd.read_csv(meta, sep='\t', names=['id', 'version', 'type', 'start', 'end', 'order', 'strand', 'transcript_end', 'annotation'])
    breakpoint_dict = get_breakpoint_dict_meta(meta_df)
    sam = pysam.AlignmentFile(bam, 'rb')
    reads_spanning, coverage_analysis, required_coverage = process_sam(sam, breakpoint_dict, min_coverage, distance)
    if out is not None:
        coverage = write_bam(sam, reads_spanning, out)
        print 'Filtered fusion reads printed to {}'.format(out)
        print coverage
    else:
        print 'No filtered bam output written.  Specify --output-file if needed.'
    sam.close()
    if analysis_out is not None:
        analyze_coverage(coverage_analysis, analysis_out)
        print 'Fusion statistics written to {}'.format(analysis_out)
    else:
        print 'No fusion statistics written.  Specify --analysis-out if needed.'
    if fusion_out is not None and det_ann is not None:
        fusion_ids = get_fusion_dict_meta(meta_df)
        detango_output(required_coverage, fusion_out, det_ann, fusion_ids)
        print 'Fusions passing filter written to {}'.format(fusion_out)
    else:
        print 'No fusions passing filter written.  Specify --fusion-out if needed.'


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Quiver Database to Fusion Reference File",
                                     description="Specify a path to a sqlite db file to create a fusion reference file"
                                                 "\nhttp://archerdx.com/software/quiver")
    parser.add_argument("--bam-file", required=True, help="Path to a bam file")
    parser.add_argument("--data-ref", required=True, help="Path to fasta metadata file (ex. Quiver_v3.17.gtf)")
    parser.add_argument("--output-file", default=None, help="Path filtered output bam")
    parser.add_argument("--detango-annotation-file", default=None, help="Path detango annotation file")
    parser.add_argument("--min-coverage", default=1, help="Minimum coverage required at N bases before and "
                                                          "after breakpoint.  (default=1)")
    parser.add_argument("--distance", default=5, help="Distance from breakpoint to analyze coverage. N bases. "
                                                      "(default=3)")
    parser.add_argument("--fusion-out", default=None, help="Path to save detango style output for "
                                                           "fusions passing filter")
    parser.add_argument("--analysis-out", default=None, help="Path to save analysis by distance")
    args = parser.parse_args()
    main(args.bam_file, args.data_ref, args.output_file, args.detango_annotation_file, int(args.min_coverage),
         args.distance, args.fusion_out, args.analysis_out)