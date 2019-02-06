# gffkit
A small lightweight toolkit for common manipulations on GFF3 files

## Usage

```
usage: gffkit.py <command> [<args>]

The commands are:
   grep                 Filter GFF3 features based on regular expressions. 
   sort                 **Not implemented.** Use genometools (https://github.com/genometools/genometools) with 'gt gff3 -tidy -sort -retainids'
   rc                   Update the extent of GFF3 features to their reverse complement
   offset               Shift GFF3 features by a constant offset
   restart              For circular references, reset the GFF3 record break (see https://github.com/shenwei356/seqkit) 'seqkit restart' for the analogous operation on FASTA files.
   subgff		get features within a specific subregion
   add_ids              Add ID lines to features that lack it 
   augustus_gtf_to_gff3 Convert augustus gtf format to GFF3
   add_name_to_fasta    Looks up the Name attribute from GFF3 & adds it to a fasta record
```

## Dependencies
* https://github.com/daler/gffutils
