import gzip

import pybedtools as pbt
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ngs_test_utils.base.fastq import FastqEntry, FastqFile



gtf = "/home/jure/Downloads/Homo_sapiens.GRCh38.92.gtf"
fasta = "/home/jure/Downloads/Homo_sapiens.GRCh38.dna.primary_assembly.fasta.gz"
bam = "/home/jure/Downloads/D701-D501_R1.bam"

gtf_out = "ann.gtf"
fasta_out = "genome.fa"
bam_out = "sub.bam"
bai_out = "sub.bam.bai"
mate1 = "mate1.fq"
mate2 = "mate2.fq"

# region
chrom = "1"
pos_min = 0
pos_max = 30000

# 1. Subsample GTF
segments = []
for segment in pbt.BedTool(gtf):
    if segment.chrom != chrom:
        continue
    if segment.end < pos_min:
        continue
    if segment.start > pos_max:
        continue
    if segment.attrs['gene_id'] not in ["ENSG00000223972", "ENSG00000227232"]:
        continue

    segments.append(segment)
pbt.BedTool(s for s in segments).saveas(gtf_out)

# 2. Subsample FASTA
with gzip.open(fasta, "rt") as handle1, open(fasta_out, "wt") as handle2:
    for entry in SeqIO.parse(handle1, "fasta"):
        if entry.name == chrom:
            entry = entry[pos_min:pos_max]
            SeqIO.write(entry, handle2, "fasta")
            break


# 3. Subsample BAM and convert to reads!

# Subsample bam!
reads ={}
with pysam.AlignmentFile(bam, mode="rb") as in_file:
    with pysam.AlignmentFile(bam_out, "wb", template=in_file) as out_file:
        for read in in_file.fetch(chrom, pos_min, pos_max):
            out_file.write(read)

            if read.qname not in reads:
                reads[read.qname] = {'R1': [], "R2": []}

            if read.is_read1:
                reads[read.qname]["R1"].append(read)
            elif read.is_read2:
                reads[read.qname]["R2"].append(read)


def write_entry(read):
    seq = read.query_sequence
    qualilties = ''.join([chr(33 + qual) for qual in read.query_qualities])
    if read.is_reverse:
        seq = str(Seq(seq).reverse_complement())
        qualilties = qualilties[::-1]

    return FastqEntry(
        seq_id="@" + read.qname,
        seq=seq,
        plus="+",
        quality=qualilties,
    )


# Write reads
ok = 0
fail = 0
fastq1 = FastqFile(mate1, "wt")
fastq2 = FastqFile(mate2, "wt")
for qname, pack in reads.items():
    r1s = pack["R1"]
    r2s = pack["R2"]
    if len(r1s) != 1 or len(r2s) != 1:
        fail += 1
        continue

    ok += 1

    read1 = r1s[0]
    read2 = r2s[0]

    entry1 = write_entry(read1)
    entry2 = write_entry(read2)

    fastq1.write([entry1])
    fastq2.write([entry2])

fastq1.close()
fastq2.close()
print("OK: ", ok)
print("fail: ", fail)
