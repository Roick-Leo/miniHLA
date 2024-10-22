import pysam
import argparse

def split_cigar(cigar):
    numl = []
    num = 0
    cigarstrl = []
    for i in cigar:
        if i.isdigit():
            num = num*10+int(i)
        else:
            numl.append(num)
            num = 0
            cigarstrl.append(i)
    return numl,cigarstrl

def SmallVar_signal_iter(ref_seq, read):
    result = []
    ref_start = int(read.reference_start)
    read_pos = 0
    ref_pos = ref_start

    # 需要储存的read信息
    ReadID = read.query_name
    cigar = read.cigarstring
    seq = read.query_sequence
    numls, cigarstrls = split_cigar(cigar)

    # numl = numl
    for num,cigarstr in zip(numls, cigarstrls):
        # 位于目标区段时
        if cigarstr == "X":
            for m in range(num):
                RefPos = ref_pos+1
                Ref = ref_seq[ref_pos]
                Alt = seq[read_pos]
                result.append("{}_SNP_{}_{}".format(RefPos,Ref,Alt))
                ref_pos += 1
                read_pos += 1
        elif cigarstr == "D":
            RefPos = ref_pos
            Ref = ref_seq[(ref_pos-1):(ref_pos+num)]
            Alt = ref_seq[ref_pos-1]
            result.append("DEL_{}_{}".format(Ref,Alt))
            pass
            ref_pos += num
        elif cigarstr == "I":
            RefPos = ref_pos
            Ref = ref_seq[ref_pos-1]
            Alt = ref_seq[ref_pos-1]+seq[read_pos:read_pos+num]
            result.append("INS_{}_{}".format(Ref,Alt))
            pass
            read_pos += num
        else:
            if cigarstr == "S":
                read_pos += num
                continue
            elif cigarstr == "N":
                ref_pos += num
                continue
            elif cigarstr in "M=":
                read_pos += num
                ref_pos += num
                continue
            else:
                continue
    return result

def main():
    parser = argparse.ArgumentParser(description='small variant signal processing')
    parser.add_argument('--bam', help='bam file')
    parser.add_argument('--fasta', help='fasta file')
    parser.add_argument('--output', help='output file')
    args = parser.parse_args()

    bamfile = args.bam
    ref_file = args.fasta
    out_put = args.output
    reference_name = ref_file.split('/')[-1].split('.')[0]

    with pysam.AlignmentFile(bamfile, "rb") as bam, pysam.FastaFile(ref_file) as ref, open(out_put, "w") as out:
        refseq = ref.fetch(reference_name)
        for read in bam.fetch():
            if read.flag == 0:
                read_id = read.query_name
                read_smallvar_signal = SmallVar_signal_iter(refseq, read)
                out.write(">{}\n{}\n".format(read_id, ";".join(read_smallvar_signal)))

if __name__ == "__main__":
    main()