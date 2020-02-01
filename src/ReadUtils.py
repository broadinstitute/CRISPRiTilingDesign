import pandas as pd
import numpy as np
from scipy import interpolate
import os
import os.path
from subprocess import check_call, check_output, PIPE, Popen, getoutput, CalledProcessError
import linecache
import traceback
import time


def runCommand(command, verbose=False):
    p = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
    if verbose: print("Running: " + command)
    (stdoutdata, stderrdata) = p.communicate()
    err = str(stderrdata, 'utf-8')
    return (stdoutdata, err)


def run_count_reads(target, output, bed_file, genome_sizes, use_fast_count, count5pEnds=False):
    if target.endswith(".bam"):
        count_bam(target, bed_file, output, genome_sizes=genome_sizes, use_fast_count=use_fast_count, count5pEnds=count5pEnds)
    elif target.endswith(".tagAlign.gz") or target.endswith(".tagAlign.bgz"):
        count_tagalign(target, bed_file, output, genome_sizes, count5pEnds=count5pEnds)
    elif isBigWigFile(target):
        if count5pEnds:
            raise ValueError("count5pEnds not supported for bigWig files")
        count_bigwig(target, bed_file, output)
    else:
        raise ValueError("File {} name was not in .bam, .tagAlign.gz, .bw".format(target))


def count_bam(bamfile, bed_file, output, genome_sizes, use_fast_count=True, verbose=False, count5pEnds=False):
    completed = True
        
    #Fast count:
    #bamtobed uses a lot of memory. Instead reorder bed file to match ordering of bam file. Assumed .bam file is sorted in the chromosome order defined by its header.
    #Then use bedtools coverage, then sort back to expected order
    #Requires an faidx file with chr in the same order as the bam file.
    if use_fast_count:
        temp_output = output + ".temp_sort_order"
        faidx_command = "awk 'FNR==NR {{x2[$1] = $0; next}} $1 in x2 {{print x2[$1]}}' {genome_sizes} <(samtools view -H {bamfile} | grep SQ | cut -f 2 | cut -c 4- )  > {temp_output}".format(**locals())
        command = "bedtools sort -faidx {temp_output} -i {bed_file}".format(**locals())
        if count5pEnds:
            command = command + " | " + collapseTo5pBaseCommand()
        command = command + " | bedtools coverage -g {temp_output} -counts -sorted -a stdin -b {bamfile} | awk '{{print $1 \"\\t\" $2 \"\\t\" $3 \"\\t\" $NF}}'  | bedtools sort -faidx {genome_sizes} -i stdin > {output}; rm {temp_output}".format(**locals())

        #executable='/bin/bash' needed to parse < redirect in faidx_command
        p = Popen(faidx_command, stdout=PIPE, stderr=PIPE, shell=True, executable='/bin/bash')
        if verbose: print("Running: " + faidx_command)
        (stdoutdata, stderrdata) = p.communicate()
        err = str(stderrdata, 'utf-8')

        p = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
        if verbose: print("Running: " + command)
        (stdoutdata, stderrdata) = p.communicate()
        err = str(stderrdata, 'utf-8')

        try:
            data = pd.read_table(output, header=None).ix[:,3].values
        except Exception as e:
            print("Fast count method failed to count: " + str(bamfile) + "\n")
            print(err)
            print("Trying bamtobed method ...\n")
            completed = False

    #Alternate counting method. Slower and requires more memory.
    # convert BAM to BED, filter to standard chromosomes, sort, then use the very fast bedtools coverage -sorted algorithm
    # Note: This requires that bed_file is also sorted and in same chromosome order as genome_sizes (first do bedtools sort -i bed_file -faidx genome_sizes)
    #         BEDTools will error out if files are not properly sorted
    # Also requires that {genome_sizes} has a corresponding {genome_sizes}.bed file
    if not use_fast_count or ("terminated"  in err) or ("Error" in err) or ("ERROR" in err) or not completed:
        command = "bedtools bamtobed -i {bamfile} | cut -f 1-3 | bedtools intersect -wa -a stdin -b {genome_sizes}.bed | bedtools sort -i stdin -faidx {genome_sizes} ".format(**locals()) 
        if count5pEnds:
            command = command + " | " + collapseTo5pBaseCommand()
        command = command + " | bedtools coverage -g {genome_sizes} -counts -sorted -a {bed_file} -b stdin | awk '{{print $1 \"\\t\" $2 \"\\t\" $3 \"\\t\" $NF}}' > {output}".format(**locals())
        
        p = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
        if verbose: print("Running: " + command)
        (stdoutdata, stderrdata) = p.communicate()

        try:
            data = pd.read_table(output, header=None).ix[:,3].values
        except Exception as e:
            print(e)
            print(stderrdata)
            completed = False

    # Check for successful finish -- BEDTools can run into memory problems
    #import pdb; pdb.set_trace()
    err = str(stderrdata, 'utf-8')
    if ("terminated" not in err) and ("Error" not in err) and ("ERROR" not in err) and any(data):
        print("BEDTools completed successfully. \n")
        completed = True
    else:
        print("BEDTools failed to count file: " + str(bamfile) + "\n")
        print(err)
        print("Trying using Java method ...\n")
        completed = False


def collapseTo5pBaseCommand():
    return "awk 'BEGIN {OFS = \"\t\" } { if ($6 == \"-\") $2 = $3 - 1; else $3 = $2 + 1; print $0 }'"


def count_tagalign(tagalign, bed_file, output, genome_sizes, count5pEnds=False):
    command1 = "tabix -B {tagalign} {bed_file} ".format(**locals())
    if count5pEnds:
        command1 = command1 + " | " + collapseTo5pBaseCommand()
    command1 = command1 + " | cut -f1-3"
    command2 = "bedtools coverage -counts -b stdin -a {bed_file} | awk '{{print $1 \"\\t\" $2 \"\\t\" $3 \"\\t\" $NF}}' ".format(**locals())
    p1 = Popen(command1, stdout=PIPE, shell=True)
    with open(output, "wb") as outfp:
        p2 = check_call(command2, stdin=p1.stdout, stdout=outfp, shell=True)

    if not p2 == 0:
        print(p2.stderr)


def count_bigwig(target, bed_file, output):
    from pyBigWig import open as open_bigwig    
    bw = open_bigwig(target)
    bed = read_bed(bed_file)
    with open(output, "wb") as outfp:
        for chr, start, end, *rest in bed.itertuples(index=False, name=None):
            # if isinstance(name, np.float):
            #     name = ""
            try:
                val = bw.stats(chr, int(start), int(max(end, start + 1)), "mean")[0] or 0
            except RuntimeError:
                print("Failed on", chr, start, end)
                raise
            val *= abs(end - start)  # convert to total coverage
            output = ("\t".join([chr, str(start), str(end), str(val)]) + "\n").encode('ascii')
            outfp.write(output)


def isBigWigFile(filename):
    return(filename.endswith(".bw") or filename.endswith(".bigWig") or filename.endswith(".bigwig"))



bed_extra_colnames = ["name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts"]
chromosomes = ['chr' + str(entry) for entry in list(range(1,23)) + ['M','X','Y']]   # should pass this in as an input file to specify chromosome order
def read_bed(filename, extra_colnames=bed_extra_colnames, chr=None, sort=False, skip_chr_sorting=False):
    skip = 1 if ("track" in open(filename, "r").readline()) else 0
    names = ["chr", "start", "end"] + extra_colnames
    result = pd.read_table(filename, names=names, header=None, skiprows=skip, comment='#')
    result = result.dropna(axis=1, how='all')  # drop empty columns
    assert result.columns[0] == "chr"

    result['chr'] = pd.Categorical(result['chr'], chromosomes, ordered=True)
    if chr is not None:
        result = result[result.chr == chr]
    if not skip_chr_sorting:
        result.sort_values("chr", inplace=True)
    if sort:
        result.sort_values(["chr", "start", "end"], inplace=True)
    return result


def read_bedgraph(filename):
    return read_bed(filename, extra_colnames=["score"], skip_chr_sorting=True)

def count_bam_mapped(bam_file):
    # Counts number of reads in a BAM file WITHOUT iterating.  Requires that the BAM is indexed
    chromosomes = ['chr' + str(x) for x in range(1,23)] + ['chrX'] + ['chrY']
    command = ("samtools idxstats " + bam_file)
    data = check_output(command, shell=True)
    lines = data.decode("ascii").split("\n")
    vals = list(int(l.split("\t")[2]) for l in lines[:-1] if l.split("\t")[0] in chromosomes)
    if not sum(vals) > 0:
        raise ValueError("Error counting BAM file: count <= 0")
    return sum(vals)

def count_tagalign_total(tagalign):
    #result = int(check_output("zcat " + tagalign + " | wc -l", shell=True))
    result = int(check_output("zcat {} | grep -E 'chr[1-9]|chr1[0-9]|chr2[0-2]|chrX|chrY' | wc -l".format(tagalign), shell=True))
    assert (result > 0)
    return result

def count_bigwig_total(bw_file):
    bw = open_bigwig(bw_file)
    result = sum(l * bw.stats(ch, 0, l, "mean")[0] for ch, l in bw.chroms().items())
    assert (abs(result) > 0)  ## BigWig could have negative values, e.g. the negative-strand GroCAP bigwigs
    return result

def count_total(infile):
    if infile.endswith(".tagAlign.gz") or infile.endswith(".tagAlign.bgz"):
        total_counts = count_tagalign_total(infile)
    elif infile.endswith(".bam"):
        total_counts = count_bam_mapped(infile)
    elif isBigWigFile(infile):
        total_counts = count_bigwig_total(infile)
    else:
        raise RuntimeError("Did not recognize file format of: " + infile)

    return total_counts
