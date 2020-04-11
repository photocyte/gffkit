#!/usr/bin/env python

##Got the multi-level argparse tip from https://chase-seibert.github.io/blog/2014/03/21/python-multilevel-argparse.html

import argparse
import sys
import pkg_resources
import Bio
import Bio.SeqIO
import re
import gzip

try:
    gffutils_version = pkg_resources.get_distribution("gffutils").version.split(".")[1]
except:
    print("This script requires gffutils version >= 0.9. See http://daler.github.io/gffutils")
    exit()

if float(gffutils_version) >= 8:
        import gffutils
else:
        print("This script requires gffutils version >= 0.9. See http://daler.github.io/gffutils")
        exit()


class GffKit(object):

    def __init__(self):
        parser = argparse.ArgumentParser(
            description='A small toolkit for common manipulations with gff3 files. By default, it prints the transformed gff3 file to the standard out, and prints informational and error messages to the standard error',
            usage='''gffkit.py <command> [<args>]

The commands are:
   grep                 Filter GFF3 features based on regular expressions.
   sort                 **Not implemented.** Use genometools (https://github.com/genometools/genometools) with 'gt gff3 -tidy -sort -retainids'
   rc                   Update the extent of GFF3 features to their reverse complement
   offset               Shift GFF3 features by a constant offset
   restart              For circular references, reset the GFF3 record break (see https://github.com/shenwei356/seqkit) 'seqkit restart' for the analogous operation on FASTA files.
   subgff		get features within a specific subregion
   add_ids              Add ID lines to features that lack it
   add_locus_tag        Take the ID attibute for a given gene feature and add it as an locus_tag attribute (used for NCBI table2asn_GFF)
   augustus_gtf_to_gff3 Convert augustus gtf format to GFF3
   add_name_to_fasta    Looks up the Name attribute from GFF3 & adds it to a fasta record
''')
        parser.add_argument('command', help='Subcommand to run')
        # parse_args defaults to [1:] for args, but you need to
        # exclude the rest of the args too, or validation will fail
        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            sys.stderr.write('Unrecognized command\n')
            parser.print_help()
            exit(1)
        # use dispatch pattern to invoke method with same name
        getattr(self, args.command)()

################################
################################
################################

    def grep(self):
        parser = argparse.ArgumentParser(
            description='Filter GFF features based on regular expressions')
        # prefixing the argument with -- means it's optional
        # now that we're inside a subcommand, ignore the first
        # TWO argvs, ie the command (git) and the subcommand (commit)
        group = parser.add_mutually_exclusive_group(required=True)
        group.add_argument('-p',metavar="string",required=False,type=str,default=None,help="The pattern to filter by")
        group.add_argument('-f',metavar="file",required=False,type=str,default=None,help="A file with patterns to filter by separated by newlines")
        parser.add_argument("file",help="The GFF3 file to filter")
        parser.add_argument("-t",metavar="type",type=str,default="gene",help="The GFF3 feature type to filter on (default=gene)")
        parser.add_argument("-v",action="store_true",default=False,help="Invert the sense of matching, to select non-matching lines.")

        args = parser.parse_args(sys.argv[2:])
        ##print 'Running git commit, amend=%s' % args.amend

        if args.p != None:
            patterns = set()
            patterns.add(args.p.strip())
        elif args.f != None and not args.f.endswith(".gz"):
            handle = open(args.f,"rU")
            patterns = set()
            for l in handle.readlines():
                patterns.add(l.strip())
        elif args.f != None and args.f.endswith(".gz"):
            handle = gzip.open(args.f,"rU")
            patterns = set()
            for l in handle.readlines():
                patterns.add(l.strip())

        db_path=args.file+".gffutils.db"
        sys.stderr.write("Reading GFF3 file: "+args.file+"\n")
        sys.stderr.write("Coverting to in-memory gffutils sqlite database.\n")
        sys.stderr.flush()

        db = gffutils.create_db(args.file,":memory:", force=True,merge_strategy="create_unique")
        sys.stderr.write("Done with conversion\n")

        ##This function doesn't mind >1 mRNA per gene
        if args.t == "gene":
            for g in db.features_of_type("gene"):
                if (args.v == True and g.id not in patterns) or (args.v == False and g.id in patterns):
                    print(g)
                    for c in db.children(g.id):
                        print(c)

        forbidden_genes = set()
        forbidden_mRNAs = set()
        genes_to_print = set()
        mRNAs_to_print = set()
        if args.t == "exon":
            sys.stderr.write("Warning. The exon mode should only be used for filtering with '-v', as currently the script does not correctly resize mRNAs/genes to reduced exons)\n")
        for e in db.features_of_type("exon"):
            if (args.v == True and e.id not in patterns):
                ##Print as per normal, but do it from an exon based perspective.
                e_parent = e.attributes["Parent"][0]
                mRNA = db[e_parent]
                m_parent = mRNA.attributes["Parent"][0]
                gene = db[m_parent]
                if gene.id not in genes_to_print:
                    genes_to_print.add(gene.id)
                if mRNA.id not in mRNAs_to_print:
                    mRNAs_to_print.add(mRNA.id)

                if (args.v == True and e.id in patterns):
                    e_parent = e.attributes["Parent"][0]
                    mRNA = db[e_parent]
                    m_parent = mRNA.attributes["Parent"][0]
                    gene = db[m_parent]
                    if gene.id not in forbidden_genes:
                        forbidden_genes.add(gene.id)
                    if mRNA.id not in forbidden_mRNAs:
                        forbidden_mRNAs.add(mRNA.id)

                elif (args.v == False and e.id in patterns):
                    sys.stderr.write("Only mode '-v' is currently supported\n")
                    break
	##Once all the exons have been examined, print those features that are not forbidden.
        allowed_genes = genes_to_print - forbidden_genes
        allowed_mRNAs = mRNAs_to_print - forbidden_mRNAs
        sys.stderr.write("Had "+str(len(genes_to_print))+" genes_to_print\n")
        sys.stderr.write("Had "+str(len(forbidden_genes))+" forbidden_genes\n")
        sys.stderr.write("Had "+str(len(mRNAs_to_print))+" mRNAs_to_print\n")
        sys.stderr.write("Had "+str(len(forbidden_mRNAs))+" forbidden_mRNAs\n")
        for g in db.features_of_type("gene"):
            if g.id in allowed_genes:
                print(g)
            else:
                continue

            for m in db.children(g.id,featuretype="mRNA"):
                if m.id in allowed_mRNAs:
                    print(m)
                    for c in db.children(m.id):
                        print(c)
                else:
                    sys.stderr.write("ERROR:This shouldn't happen.\n")
                    exit(11)
    sys.stderr.write("Done.\n")

################################
################################
################################

    def rc(self):
        parser = argparse.ArgumentParser(
            description='subcommand:rc Convert GFF features to their reverse complement')
        # NOT prefixing the argument with -- means it's not optional
        requiredNamed = parser.add_argument_group('required named arguments')
        requiredNamed.add_argument("-g",metavar="<GFF3_FILE>",required=True,help="The path to the GFF3 file to reverse complement")
        requiredNamed.add_argument("-f",metavar="<FASTA_FILE>",required=True,help="The path to the FASTA file to get record lengths from")
        parser.add_argument("--name_to_id",default=False,action='store_true',help="Workaround.  Set the name attribute to the ID attribute")
        parser.add_argument("--add_rc",default=False,action='store_true',help="Add '_rc' to the end of the chrom IDs, used for chrom id matching and feature modification")
        args = parser.parse_args(sys.argv[2:])

        ##print 'Running git fetch, repository=%s' % args.repository

        db_path=args.g+".gffutils.db"
        sys.stderr.write("Reading GFF3 file: "+args.g+"\n")
        sys.stderr.write("Coverting to in memory gffutils sqlite database\n")
        sys.stderr.flush()
        db = gffutils.create_db(args.g,":memory:", force=True,merge_strategy="create_unique")
        sys.stderr.write("Done converting. Now printing modified GFF3 to stdout...\n")
        sys.stderr.flush()

        errorOccurances = dict()
        errorOccurances[1] = 0

        reference_lengths = dict()
        handle = Bio.SeqIO.parse(args.f,"fasta")
        for record in handle:
            reference_lengths[record.id] = len(record.seq)

        sys.stdout.write("##gff-version 3\n")
        for f in db.all_features():

            new_attrs = []
            for a in f.attributes:
                 new_attrs.append(a+"="+f.attributes[a][0].replace("=",":")) ##Replacement because genome tools doesn't like extra '='
            new_attr_string = ";".join(new_attrs)

            if args.add_rc == True:
                f.chrom += '_rc'

            if f.chrom not in reference_lengths.keys():
                if errorOccurances[1] < 10:
                    sys.stderr.write("Warning: feature chrom "+f.chrom+" is not in fasta file. printing feature without modification\n")
                elif errorOccurances[1] == 10:
                    sys.stderr.write("Warning: (silencing rest of warnings)\n")
                if args.add_rc == True:
                    f.chrom = f.chrom[:-3] ##If '_rc' is added on to check, have to remove it since this line should be unmodified.
                gene_string = '\t'.join([f.chrom,f.source,f.featuretype,str(f.start),str(f.stop),f.score,f.strand,f.frame,new_attr_string])
                sys.stdout.write(gene_string+"\n")
                errorOccurances[1] += 1
                continue

            new_start = reference_lengths[f.chrom] - f.start + 1
            new_stop = reference_lengths[f.chrom] - f.stop + 1
            f.start = new_stop
            f.stop = new_start
            if args.name_to_id:
                f.attributes["Name"] = f.id

            if f.strand == "+":
                f.strand = "-"
            elif f.strand == "-":
                f.strand = "+"


            gene_string = '\t'.join([f.chrom,f.source,f.featuretype,str(f.start),str(f.stop),f.score,f.strand,f.frame,new_attr_string])
            sys.stdout.write(gene_string+"\n")

        sys.stderr.write("Conversion complete.\n")

#################
#################
#################

    def offset(self):
        parser = argparse.ArgumentParser(
            description='subcommand:offset Shift GFF features by a constant offset')
        requiredNamed = parser.add_argument_group('required named arguments')
        requiredNamed.add_argument("-g",metavar="example.gff3",help="The path to the GFF3 file to offset",required=True)
        requiredNamed.add_argument("-i",metavar="int",type=int,default=0,help="The number of nucleotides to offet the GFF by. Positive = right, negative = left",required=True)
        parser.add_argument("-o",default=False,action='store_true',help="Workaround.  Set the name attribute to the ID attribute")
        parser.add_argument("--name_to_id",default=False,action='store_true',help="Workaround.  Set the name attribute to the ID attribute")
        args = parser.parse_args(sys.argv[2:])

        db_path=args.g+".gffutils.db"
        sys.stderr.write("Reading GFF3 file: "+args.g+"\n")
        sys.stderr.write("Coverting to in memory gffutils sqlite database\n")
        sys.stderr.flush()
        db = gffutils.create_db(args.g,":memory:", force=True,merge_strategy="create_unique")
        sys.stderr.write("Done converting. Now printing modified GFF3 to stdout...\n")
        sys.stderr.flush()

        sys.stdout.write("##gff-version 3\n")
        for f in db.all_features():
            new_start = f.start + args.i
            new_stop = f.stop + args.i
            f.start = new_start
            f.stop = new_stop

            if args.name_to_id:
                f.attributes["Name"] = f.id
            new_attrs = []
            for a in f.attributes:
                new_attrs.append(a+"="+f.attributes[a][0])
            new_attr_string = ";".join(new_attrs)
            gene_string = '\t'.join([f.chrom,f.source,f.featuretype,str(f.start),str(f.stop),f.score,f.strand,f.frame,new_attr_string])
            sys.stdout.write(gene_string+"\n")


        sys.stderr.write("Conversion complete.\n")

#################
#################
#################

    def add_ids(self):
        parser = argparse.ArgumentParser(
            description='subcommand:add_ids add ID attributes to the features that don\'t have them')
        requiredNamed = parser.add_argument_group('required named arguments')
        requiredNamed.add_argument("-g",metavar="example.gff3",help="The path to the GFF3 file to offset",required=True)
        args = parser.parse_args(sys.argv[2:])

        db_path=args.g+".gffutils.db"
        sys.stderr.write("Reading GFF3 file: "+args.g+"\n")
        sys.stderr.write("Coverting to in memory gffutils sqlite database\n")
        sys.stderr.flush()
        db = gffutils.create_db(args.g,":memory:", force=True,merge_strategy="create_unique")
        sys.stderr.write("Done converting. Now printing modified GFF3 to stdout...\n")
        sys.stderr.flush()

        sys.stdout.write("##gff-version 3\n")
        z=0
        for f in db.all_features():
            if "ID" not in f.attributes.keys():
                parentID = list(db.parents(f.id))[0].id
                if f.strand == "+":
                    siblings = list(db.children(parentID, featuretype=f.featuretype, order_by='start'))
                elif f.strand == "-":
                    siblings = list(db.children(parentID, featuretype=f.featuretype, order_by='end', reverse=True))
                siblingsTotal = len(siblings)
                featureIndex = siblings.index(f)+1
                f.attributes["ID"] = parentID+"_"+f.featuretype+str(featureIndex)
            ##if args.name_to_id:
            ##    f.attributes["Name"] = f.id

            new_attrs = []
            for a in f.attributes:
                new_attrs.append(a+"="+f.attributes[a][0])
            new_attr_string = ";".join(new_attrs)
            gene_string = '\t'.join([f.chrom,f.source,f.featuretype,str(f.start),str(f.stop),f.score,f.strand,f.frame,new_attr_string])
            sys.stdout.write(gene_string+"\n")


        sys.stderr.write("Conversion complete.\n")

#################
#################
#################

    def add_locus_tag(self):
        parser = argparse.ArgumentParser(
            description='subcommand:add_locus_tag Take the ID attibute for a given gene feature and add it as an locus_tag attribute (used for NCBI table2asn_GFF)')
        requiredNamed = parser.add_argument_group('required named arguments')
        requiredNamed.add_argument("-g",metavar="example.gff3",help="The path to the GFF3 file to offset",required=True)
        args = parser.parse_args(sys.argv[2:])

        db_path=args.g+".gffutils.db"
        sys.stderr.write("Reading GFF3 file: "+args.g+"\n")
        sys.stderr.write("Coverting to in memory gffutils sqlite database\n")
        sys.stderr.flush()
        db = gffutils.create_db(args.g,":memory:", force=True,merge_strategy="create_unique")
        sys.stderr.write("Done converting. Now printing modified GFF3 to stdout...\n")
        sys.stderr.flush()

        sys.stdout.write("##gff-version 3\n")
        z=0
        for f in db.all_features():
            new_attrs = []
            for a in f.attributes:
                new_attrs.append(a+"="+f.attributes[a][0])
            if f.featuretype == "gene":
                f.attributes["locus_tag"] = f.id
                new_attrs.append("locus_tag="+f.id)
            else:
                pass
            new_attr_string = ";".join(new_attrs)
            gene_string = '\t'.join([f.chrom,f.source,f.featuretype,str(f.start),str(f.stop),f.score,f.strand,f.frame,new_attr_string])
            sys.stdout.write(gene_string+"\n")

        sys.stderr.write("Conversion complete.\n")

#################
#################
#################

    def subgff(self):
        parser = argparse.ArgumentParser(
            description='subcommand:subgff get features by region of the GFF3')
        requiredNamed = parser.add_argument_group('required named arguments')
        requiredNamed.add_argument("-g",metavar="example.gff3",help="The path to the GFF3 file to get regions from",required=True)
        requiredNamed.add_argument("-r",metavar="range",type=str,default=0,help="Range to filter by, colon delimited. E.g. 3400:3900",required=True)
        requiredNamed.add_argument("-x",default=False,action='store_true',help="Reset coordinate numbering so that range[0] = 1")
        args = parser.parse_args(sys.argv[2:])

        db_path=args.g+".gffutils.db"
        sys.stderr.write("Reading GFF3 file: "+args.g+"\n")
        sys.stderr.write("Coverting to in memory gffutils sqlite database\n")
        sys.stderr.flush()
        db = gffutils.create_db(args.g,":memory:", force=True,merge_strategy="create_unique")
        sys.stderr.write("Done converting. Now printing modified GFF3 to stdout...\n")
        sys.stderr.flush()

        range_re = re.compile("([0-9]+):([0-9]+)")
        range_result = range_re.search(args.r)
        if range_result == None:
            sys.stderr.write("subgff: Couldn't parse region string.")
            sys.stderr.flush()
            exit()
        subgff_start = int(range_result.group(1))
        subgff_end = int(range_result.group(2))

        assert subgff_start < subgff_end


        sys.stdout.write("##gff-version 3\n")
        for f in db.all_features():
            if f.start > subgff_start and f.end < subgff_end:
                new_attrs = []
                for a in f.attributes:
                    new_attrs.append(a+"="+f.attributes[a][0])
                    new_attr_string = ";".join(new_attrs)
                if args.x == True:
                    f.start = f.start - subgff_start + 1
                    f.end = f.end - subgff_start + 1
                gene_string = '\t'.join([f.chrom,f.source,f.featuretype,str(f.start),str(f.stop),f.score,f.strand,f.frame,new_attr_string])
                sys.stdout.write(gene_string+"\n")


        sys.stderr.write("Conversion complete.\n")

################################
################################
################################

    def augustus_gtf_to_gff3(self):
        parser = argparse.ArgumentParser(
            description='subcommand:augustus_gtf_to_gff3 Converts the augustus gtf (default) format, to GFF3')
        requiredNamed = parser.add_argument_group('required named arguments')
        requiredNamed.add_argument("-g",metavar="example.gff3",help="The path to the GFF3 file to offset",required=True)
        parser.add_argument("--make_long_id",default=False,action='store_true',help="Append the record name to the gene & transcript ID")
        args = parser.parse_args(sys.argv[2:])
##Example GFF3
##Ilumi1.1_Scaffold12837	EVM	gene	2913	3323	.	+	.	ID=ILUMI_05824;Name=EVM prediction Ilumi1.1_Scaffold12837.1

#Example GTF from augustus
#NODE_64342589_length_12508_cov_11.070071        AUGUSTUS        gene    1814    8945    0.01    -       .       g2
#NODE_64342589_length_12508_cov_11.070071        AUGUSTUS        transcript      1814    8945    0.01    -       .       g2.t1
#NODE_64342589_length_12508_cov_11.070071        AUGUSTUS        tts     1814    1814    .       -       .       transcript_id "g2.t1"; gene_id "g2";
#NODE_64342589_length_12508_cov_11.070071        AUGUSTUS        exon    1814    2912    .       -       .       transcript_id "g2.t1"; gene_id "g2";
#NODE_64342589_length_12508_cov_11.070071        AUGUSTUS        stop_codon      2050    2052    .       -       0       transcript_id "g2.t1"; gene_id "g2";
#NODE_64342589_length_12508_cov_11.070071        AUGUSTUS        intron  2913    3279    0.42    -       .       transcript_id "g2.t1"; gene_id "g2";
#NODE_64342589_length_12508_cov_11.070071        AUGUSTUS        intron  3434    3661    0.53    -       .       transcript_id "g2.t1"; gene_id "g2";
#NODE_64342589_length_12508_cov_11.070071        AUGUSTUS        intron  4333    4671    0.62    -       .       transcript_id "g2.t1"; gene_id "g2";
#NODE_64342589_length_12508_cov_11.070071        AUGUSTUS        intron  4996    5244    0.99    -       .       transcript_id "g2.t1"; gene_id "g2";
#NODE_64342589_length_12508_cov_11.070071        AUGUSTUS        intron  5650    5751    0.41    -       .       transcript_id "g2.t1"; gene_id "g2";
#NODE_64342589_length_12508_cov_11.070071        AUGUSTUS        intron  5946    6076    0.53    -       .       transcript_id "g2.t1"; gene_id "g2";
#NODE_64342589_length_12508_cov_11.070071        AUGUSTUS        CDS     2050    2912    0.42    -       2       transcript_id "g2.t1"; gene_id "g2";
#NODE_64342589_length_12508_cov_11.070071        AUGUSTUS        CDS     3280    3433    0.69    -       0       transcript_id "g2.t1"; gene_id "g2";
#NODE_64342589_length_12508_cov_11.070071        AUGUSTUS        exon    3280    3433    .       -       .       transcript_id "g2.t1"; gene_id "g2";
#NODE_64342589_length_12508_cov_11.070071        AUGUSTUS        CDS     3662    4332    0.63    -       2       transcript_id "g2.t1"; gene_id "g2";
#NODE_64342589_length_12508_cov_11.070071        AUGUSTUS        exon    3662    4332    .       -       .       transcript_id "g2.t1"; gene_id "g2";
#NODE_64342589_length_12508_cov_11.070071        AUGUSTUS        CDS     4672    4995    0.64    -       2       transcript_id "g2.t1"; gene_id "g2";
#NODE_64342589_length_12508_cov_11.070071        AUGUSTUS        exon    4672    4995    .       -       .       transcript_id "g2.t1"; gene_id "g2";
#NODE_64342589_length_12508_cov_11.070071        AUGUSTUS        CDS     5245    5649    0.42    -       2       transcript_id "g2.t1"; gene_id "g2";
#NODE_64342589_length_12508_cov_11.070071        AUGUSTUS        exon    5245    5649    .       -       .       transcript_id "g2.t1"; gene_id "g2";
#NODE_64342589_length_12508_cov_11.070071        AUGUSTUS        CDS     5752    5945    0.53    -       1       transcript_id "g2.t1"; gene_id "g2";
#NODE_64342589_length_12508_cov_11.070071        AUGUSTUS        exon    5752    5945    .       -       .       transcript_id "g2.t1"; gene_id "g2";
#NODE_64342589_length_12508_cov_11.070071        AUGUSTUS        CDS     6077    8754    0.47    -       0       transcript_id "g2.t1"; gene_id "g2";
#NODE_64342589_length_12508_cov_11.070071        AUGUSTUS        exon    6077    8945    .       -       .       transcript_id "g2.t1"; gene_id "g2";
#NODE_64342589_length_12508_cov_11.070071        AUGUSTUS        start_codon     8752    8754    .       -       0       transcript_id "g2.t1"; gene_id "g2";
#NODE_64342589_length_12508_cov_11.070071        AUGUSTUS        tss     8945    8945    .       -       .       transcript_id "g2.t1"; gene_id "g2";

        # Based off Jorvis's convert_augustus_to_gff3.py script:
        #   If GTF is detected, let's start by transforming the 9th column into GFF so the
        #   libraries can use it
        #   g1  ->  ID=g1
        #   g1.t1  ->  ID=g1.t1;Parent=g1
        #   transcript_id "g1.t1"; gene_id "g1";  ->  ID=g1.t1.cds;Parent=g1.t1
        col8_gene_regex = re.compile('(g\d+)') ##Group 1 == gene id
        col8_transcript_regex = re.compile('((g\d+).t\d+)') ##Group 1 == transcript id, group 2 == gene id
        col8_trancript_id_regex = re.compile('transcript_id "(g\d+.t\d+)"; gene_id "g\d+";') ##Group 1 == transcript id

        gene_dict = dict()

        sys.stdout.write("##gff-version 3\n")
        handle = open(args.g,"rU")
        for line in handle:
            if line[0] == "#":
                continue
            splitline = line.split("\t")
            if splitline[2] not in ["exon","CDS","transcript","gene"]:
                continue
            if splitline[2] == "gene":
                gene_match = col8_gene_regex.search(splitline[8])
                gene_id = gene_match.group(1)
                if args.make_long_id == True:
                    gene_id = splitline[0]+"-"+gene_id
                attr_str = "ID="+gene_id
                sys.stdout.write("\t".join(splitline[0:8]+[attr_str])+"\n")

            if splitline[2] == "transcript":
                splitline[2] = "mRNA" ##GFF3 is mRNA not transcript
                transcript_match = col8_transcript_regex.search(splitline[8])
                transcript_id = transcript_match.group(1)
                gene_id = transcript_match.group(2)
                if args.make_long_id == True:
                    gene_id = splitline[0]+"-"+gene_id
                    transcript_id = splitline[0]+"-"+transcript_id
                attr_str = ";".join(["ID="+transcript_id,"Parent="+gene_id])
                sys.stdout.write("\t".join(splitline[0:8]+[attr_str])+"\n")

            if splitline[2] == "exon" or splitline[2] == "CDS":
                transcript_id_match = col8_trancript_id_regex.search(splitline[8])
                transcript_id = transcript_id_match.group(1)
                if args.make_long_id == True:
                    transcript_id = splitline[0]+"-"+transcript_id
                attr_str = "Parent="+transcript_id
                sys.stdout.write("\t".join(splitline[0:8]+[attr_str])+"\n")

################################
################################
################################

    def add_name_to_fasta(self):
        parser = argparse.ArgumentParser(
            description='subcommand:add_name_to_fasta Looks up the name from GFF3 & adds it to a fasta record')
        requiredNamed = parser.add_argument_group('required named arguments')
        requiredNamed.add_argument("-g",metavar="example.gff3",help="The path to the GFF3 file to parse",required=True)
        requiredNamed.add_argument("-f",metavar="example.fasta",help="The path to the Fasta file to add the GFF3 names to",required=True)
        args = parser.parse_args(sys.argv[2:])
        sys.stderr.write("subcommand incomplete. exiting...+\n")
        db_path=args.g+".gffutils.db"
        sys.stderr.write("Reading GFF3 file: "+args.g+"\n")
        sys.stderr.write("Coverting to in memory gffutils sqlite database\n")
        sys.stderr.flush()
        db = gffutils.create_db(args.g,":memory:", force=True,merge_strategy="create_unique")
        sys.stderr.write("Done converting\n")
        sys.stderr.flush()

        handle = Bio.SeqIO.parse(args.f,"fasta")
        for record in handle:
            try:
                db[record.id]
            except:
                pass
            sys.stdout.write(str(db[record.id])+'\n')

        sys.stderr.write("Renaming complete.\n")



if __name__ == '__main__':
    GffKit()
