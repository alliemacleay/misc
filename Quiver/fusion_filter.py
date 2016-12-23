import pysam
import argparse

CIGAR_CODES = {0: 'M',
               1: 'I',
               2: 'D',
               3: '',
               4: 'S',
               5: 'H'}


def match_spans(bp, pos, cigar):
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
        pos += cig[1]
    return False


def get_breakpoint_dict(bpf):
    """ Create breakpoint dictionary from file """
    bp_dict = {}
    with open(bpf, 'rb') as bp:
        for line in bp:
            fusion, breakpoint = line.split('\t')
            breakpoint = breakpoint.strip()
            if not breakpoint == '':
                bp_dict[fusion] = int(breakpoint)
            else:
                print 'No valid breakpoint for {}'.format(fusion)
    return bp_dict


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


def main(bam, bpf, out):
    """ Main entry point """
    reads_spanning = {}
    breakpoint_dict = get_breakpoint_dict(bpf)
    tot = len(breakpoint_dict)
    count = 0
    txt_found = ''
    sam = pysam.AlignmentFile(bam, 'rb')
    for fusion in breakpoint_dict:
        count += 1
        if len(reads_spanning) > 0:
            txt_found = ' {} valid fusions found'.format(len(reads_spanning))
        print 'Parsing reads for {}.  ({}/{}){}'.format(fusion, count, tot, txt_found)
        reads_spanning[fusion] = []
        breakpoint = breakpoint_dict[fusion]
        for read in sam.fetch(fusion, 1, 1000000):
            if read.cigartuples is not None and match_spans(breakpoint, read.pos, read.cigartuples):
                reads_spanning[fusion].append(read)
        if len(reads_spanning[fusion]) == 0:
            del reads_spanning[fusion]
    if out is not None:
        coverage = write_bam(sam, reads_spanning, out)
        print 'Filtered fusion reads printed to {}'.format(out)
        print coverage
    else:
        print 'No output written.  Specify --output-file if needed.'
    sam.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Quiver Database to Fusion Reference File",
                                     description="Specify a path to a sqlite db file to create a fusion reference file"
                                                 "\nhttp://archerdx.com/software/quiver")
    parser.add_argument("--bam-file", required=True, help="Path to a bam file")
    parser.add_argument("--breakpoint-file", required=True, help="Path to breakpoints file")
    parser.add_argument("--output-file", default=None, help="Path filtered output bam")
    args = parser.parse_args()
    main(args.bam_file, args.breakpoint_file, args.output_file)