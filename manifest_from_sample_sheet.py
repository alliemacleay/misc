import sys
import os
import pandas as pd

HEAD = """##fileformat=SOMATIC,1.0
IN_PIPELINE     OUT_PIPELINE    BATCH_ID        CID_ID  NAE_ID  SPECIMEN_ID     SPECIMEN_TYPE   P5_BARCODE      P7_BARCODE      RUN_FOLDER      R1      R2      R3      R4      BAM     VCF\n"""

def replace_chars(cell):
    return str(cell).replace('_', '')

def parse_sample_sheet(sheet, delim):
    cols = ['Sample_ID', 'Sample_Name', 'I7_Index_ID', 'index', 'I5_Index_ID', 'index2']
    skip = 1
    ext = sheet.split('.')[-1]
    new_sh = sheet.replace(ext, "2.{}".format(ext))
    with open(new_sh, 'w') as nsh:
        with open(sheet, 'r') as sh:
            for line in sh:
                nsh.write(line)
                if line.startswith("Project Name"):
                    batch_name = line.split(",")[1]
                if line.startswith("Description"):
                    folder = line.split(",")[1]
                if not line.startswith("[Data]"):
                    skip += 1
                else:
                    break
        df = pd.read_csv(sheet, sep=delim, skiprows=skip, dtype=str)
    df.applymap(replace_chars)
    df[cols].to_csv(new_sh, mode='a', sep=',', index=False)        
    return batch_name, folder, df

def make_manifest(ss_df, outf, run_folder, batch_id):
    inp = "bcl2fastq2"
    oup = "solid_fusion"
    with open(outf, 'w') as fh:
        fh.write(HEAD)
        for i, row in ss_df.iterrows():
            spec_id = str(row['Sample_ID'])
            cid_id = str(row['Sample_Name'])
            nae_id = str(row['Sample_Name'])
            p5 = str(row['I5_Index_ID'])
            p7 = str(row['I7_Index_ID'])
            # [u'Sample_ID', u'Sample_Name', u'Sample_Plate', u'Sample_Well', u'I7_Index_ID', u'index', u'I5_Index_ID', u'index2', u'Sample_Project', u'Description', u'GenomeFolder'],
            #5,Prof-5NY,Plate,Well,N702,CGTACTAG,A21,TGGAGAGG,170829_M03634_0273_000000000-B8LLD,Description,
            line = "{}\n".format("\t".join([inp, oup, batch_id, cid_id, nae_id, spec_id, "tumor", p5, p7, run_folder, "", "", "", "", "", ""]))
            fh.write(line)

if __name__ == "__main__":
    sample_sheet = sys.argv[1]
    ext = sample_sheet.split(".")[-1]
    delim = ","
    seq_folder = "/data/ngscid-clinical/sequencer"
    batch_name, folder_name, ss = parse_sample_sheet(sample_sheet, delim)
    run_folder = os.path.join(seq_folder, folder_name)
    make_manifest(ss, sample_sheet.replace(ext, "manifest"), run_folder, batch_name)
