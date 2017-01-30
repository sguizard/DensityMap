# DensityMap

It can be usefull to visualize the presence and the density of genomics feature along chromosomes.

DensityMap is perl tool that automatically compute the percentage of feature covered bases on genomic window along chromosomes (or any sequences).

It produce a map of density for each sequences and for each feature.
The output SVG can be customized by setting margins, densityMaps widths, etc ...

If a fasta file of the plotted sequence is given as input, the script can also plot a densityMap the percentage of GC by window.

## Requierements
DensityMap rely only on two Perl modules:
* PerlIO::gzip
* GD::SVG

The modules can be installed throught linux package manager:
```
# On debian based systeme (ubuntu, Linux mint, ...)
sudo apt-get install libgd-svg-perl libperlio-gzip-perl

# On Red Hat based systeme (Fedora, ...)
sudo dnf install perl-GD-SVG perl-PerlIO-gzip
```

## Installation
First grab DensityMap source code on github:
[Release v1.0](https://github.com/sguizard/DensityMap/releases/tag/v1.0)

Unzip and untar the archive
```
tar -zxvf DensityMap-1.0.tar.gz
```

DensityMap can launched from uncompressed directory or can be installed by executing the installation script:
```
sudo ./install.sh
```
It will create a directory DensityMap in /usr/local/share and copy all files in it and create a symbolic link into /usr/local/bin.

## Quick start
DensityMap is delivered this test files:
* dmel.fa: the genomic sequence of drosophila melanogaster
* dmel.gff3: Annotation file of drosophilia melanogaster genes, CDS, transposables elements
* dmel.bed: region file defining which region of the chromosome to plot and the order of the chromosomes

```
head dmel.fa
>2L
CGACAATGCACGACAGAGGAAGCAGAACAGATATTTAGATTGCCTCTCATTTTCTCTCCCATATTATAGGGAGAAATATG
ATCGCGTATGCGAGAGTAGTGCCAACATATTGTGCTCTTTGATTTTTTGGCAACCCAAAATGGTGGCGGATGAACGAGAT
GATAATATATTCAAGTTGCCGCTAATCAGAAATAAATTCATTGCAACGTTAAATACAGCACAATATATGATCGCGTATGC
GAGAGTAGTGCCAACATATTGTGCTAATGAGTGCCTCTCGTTCTCTGTCTTATATTACCGCAAACCCAAAAAGACAATAC
ACGACAGAGAGAGAGAGCAGCGGAGATATTTAGATTGCCTATTAAATATGATCGCGTATGCGAGAGTAGTGCCAACATAT
TGTGCTCTCTATATAATGACTGCCTCTCATTCTGTCTTATTTTACCGCAAACCCAAATCGACAATGCACGACAGAGGAAG
CAGAACAGATATTTAGATTGCCTCTCATTTTCTCTCCCATATTATAGGGAGAAATATGATCGCGTATGCGAGAGTAGTGC
CAACATATTGTGCTCTTTGATTTTTTGGCAACCCAAAATGGTGGCGGATGAACGAGATGATAATATATTCAAGTTGCCGC
TAATCAGAAATAAATTCATTGCAACGTTAAATACAGCACAATATATGATCGCGTATGCGAGAGTAGTGCCAACATATTGT
...

head dmel.gff3
##gff-version 3
##sequence-region 2L 1 23513712
##species http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=7227
2L      RefSeq  gene    7529    9484    .       +       .       ID=gene2671;Dbxref=FLYBASE:FBgn0031208,GeneID:33155;Name=CG11023;gbkey=Gene;gene=CG11023;gene_biotype=protein_coding;gene_synonym=Dmel\CG11023;locus_tag=Dmel_CG11023
2L      RefSeq  gene    9839    21376   .       -       .       ID=gene2672;Dbxref=FLYBASE:FBgn0002121,GeneID:33156;Name=l(2)gl;description=lethal (2) giant larvae;gbkey=Gene;gene=l(2)gl;gene_biotype=protein_coding;gene_synonym=CG2671,D-LGL,dlgl,Dmel\CG2671,gl,l(2),l(2) giant larva,l(2)giant larvae,L(2)gl,L(2)GL,l-gl,lgl,Lgl,LGL,l[[2]]gl,MENE (2L)-B,p127,p127l(2)gl,p127[l(2)gl];locus_tag=Dmel_CG2671
2L      RefSeq  gene    21823   25155   .       -       .       ID=gene2673;Dbxref=FLYBASE:FBgn0031209,GeneID:33157;Name=Ir21a;description=Ionotropic receptor 21a;gbkey=Gene;gene=Ir21a;gene_biotype=protein_coding;gene_synonym=CG2657,CT8983,DmelIR21a,Dmel\CG2657,ir21a,IR21a;locus_tag=Dmel_CG2657                                                                                                                             
2L      RefSeq  gene    21952   24237   .       +       .       ID=gene2674;Dbxref=FLYBASE:FBgn0263584,GeneID:12797867;Name=CR43609;gbkey=Gene;gene=CR43609;gene_biotype=ncRNA;gene_synonym=Dmel\CR43609;locus_tag=Dmel_CR43609                                                                                                                                                                                                     
2L      RefSeq  gene    25402   65404   .       -       .       ID=gene2675;Dbxref=FLYBASE:FBgn0051973,GeneID:33158;Name=Cda5;description=Chitin deacetylase-like 5;gbkey=Gene;gene=Cda5;gene_biotype=protein_coding;gene_synonym=BcDNA:RH43162,CG2761,CG2776,CG31973,DmCDA5,Dmel\CG31973;locus_tag=Dmel_CG31973                                                                                                                    
2L      RefSeq  gene    54817   55767   .       +       .       ID=gene2676;Dbxref=FLYBASE:FBgn0267987,GeneID:26067564;Name=CR46254;gbkey=Gene;gene=CR46254;gene_biotype=ncRNA;gene_synonym=Dmel\CR46254;locus_tag=Dmel_CR46254                                                                                                                                                                                                     
2L      RefSeq  gene    65999   66242   .       +       .       ID=gene2677;Dbxref=FLYBASE:FBgn0266878,GeneID:19834983;Name=CR45339;gbkey=Gene;gene=CR45339;gene_biotype=ncRNA;gene_synonym=Dmel\CR45339;locus_tag=Dmel_CR45339  
...

cat dmel.bed
2L      0       23513712                                                                                                                                                                                          
2R      0       25286936                                                                                                                                                                                          
3L      0       28110227                                                                                                                                                                                          
3R      0       32079331                                                                                                                                                                                          
4       0       1348131                                                                                                                                                                                           
X       0       23542271                                                                                                                                                                                          
Y       0       3667352 
```

You can test the program by executing the following command:
```
DensityMap.pl -i dmel.gff3 -re dmel3.bed -o dmel.svg -ty "type=LINE&strand=fused&lab=LINE" -ty "key=Target&val=Gypsy&strand=fused&cs=10&ro=ceil&lab=Gypsy" -ba white -sc 20000 -sh 100 -title "LTR and LINE retrotransposon in Dmel genome" -for -v
gffs: dmel.gff3                                                                                                                                                                                                   
region_file: dmel3.bed                                                                                                                                                                                            
output_img_name: dmel.svg                                                                                                                                                                                         

Reading region file ...
Reading fasta dmel.fa ...
Create Picture ...
Load colors ...
Draw background ...
Drawing title ...
Drawing Scale ...
Reading GFF files ...
Ploting process ...
=> Processing 2L ...
==> Processing type LINE ...
==> Processing type Target=Gypsy ...
=> Processing 2R ...
==> Processing type LINE ...
==> Processing type Target=Gypsy ...
=> Processing 3L ...
==> Processing type LINE ...
==> Processing type Target=Gypsy ...
=> Processing 3R ...
==> Processing type LINE ...
==> Processing type Target=Gypsy ...
=> Processing 4 ...
==> Processing type LINE ...
==> Processing type Target=Gypsy ...
=> Processing X ...
==> Processing type LINE ...
==> Processing type Target=Gypsy ...
=> Processing Y ...
==> Processing type LINE ...
==> Processing type Target=Gypsy ...
Image saved ! 
```

Command explaination:
* inputs and outputs options:
 * -i dmel.gff3 -re dmel3.bed -o dmel.svg
* Features declaration
 * -ty "type=LINE&strand=fused" -ty "key=Target&val=Gypsy&strand=fused&cs=10&ro=ceil&lab=GYPSY"
 * It declare two features
  * LINE, selected on the third column of the gff with no distinction between intervals on minus or plus strand
  * Gypsy, selected on the ninth column with no distinction between intervals on minus or plus strand, the density colors will be defined by the color scale 10, the densities will rounded to the lower integer, the name Gypsy will be replaced by GYPSY
* Graphical options
 * -ba white -sc 20000 -sh 100 -title "LTR and LINE retrotransposon in Dmel genome" -for -v
 * The densities will be computed on genomic windows of 2000 kbp


This command will produce this image: 

![The obtained image](https://github.com/sguizard/DensityMap/blob/dev/dmel/dmel_1.svg)

It's possible to plot a densityMap of GC% by adding the option -gc {color_scale number} and to put the dasta sequence as input by using the -fasta option:

```
DensityMap.pl -i dmel.gff3 -re dmel.bed -fa dmel.fa -o dmel.svg -ty "type=LINE&strand=fused&lab=LINE" -ty "key=Target&val=Gypsy&strand=fused&cs=10&ro=ceil&lab=GYPSY" -ba white -sc 20000 -sh 100 -title "LTR and LINE retrotransposon in Dmel genome" -gc 12 -for
```

![The obtained image](https://github.com/sguizard/DensityMap/blob/dev/dmel/dmel_2.svg)

TO WRITE: Region Ploting
TO WRITE: Color scale customization


## Detailed options

### Inputs
* 1 or severall gff files (-input)
 * formatted following the [SequenceOntolgy](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md) rules
 * It define the position and caracteristics of features along chromosomes
 * **Mandatory option**
* 1 bed file (-region_file)
 * formatted folowing rules defined at [UCSC](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)
 * For each sequence/chromosomes, define the start and the end of the region to plot
 * It also define in which order the sequence will be plotted
 * **Mandatory option**
* 1 or severall fasta file (-fasta)
 * Contain the nucleotides sequences of chromosomes to plot
 * Optional
 
### Ouput
* An SVG image (-o)
* A table file (tabulate seperated value file)
 * Describe the density (0-100 %) for each feature, for each window

### General option
* -h | help
 * Show the help message
* -v | verbose
 * Script will describe what it do
* -for | force
 * Automatically validate the image size
* -ty | type_to_draw
 * This option is for defining the feature and their density option
 * Format: **"**Type**&**Rounding_Method**&**Strand&colorScale**&**Label**"**
   * Type:
     * Gff 3rd column -> type
     * Gff 9rd column -> key & val
   * Strand:
     * \- = strand -
     * \+ = strand +
     * both = strand - and strand +
     * fused = Combination of strand - and strand +
     * all = strand - and strand + and fused
   * Rounding method
     * floor rounding density to the lowest integer
     * ceil rounding density to the upper integer
   * Label
     * Overide feature name
 * Example 1: -ty "type=match&str=all" -ty "type=gene&str=both" -ty "type=CDS&str=fused&label=cds"
 * Example 2: -ty "key=ID&val=transposon&str=all&cs=7&ro=ceil" -ty "type=gene&str=both&cs=7"

### Density options
* -c | colour_scale
 * Color scale to use for the densities
 * Default = 7
* -sc | scale_factor
 * window length in bp
 * Default = 1000
* -a | auto_scale_factor
 * Max image height in pixel
* -ro | rounding_method (for all features)
 * floor or ceil
 * Default = floor
* -gc
 * add a density map of the GC%
 * Default = 12
 * fasta option must be set

### Graphical options
* -ti | title
 * Title to print on the picture
* -w | win_size
 * Height of window in pixel
 * Default = 1
* -sh | show_scale
 * Draw scale, the integer indicate the maximum tick to print on the scale
 * Default = 50
* -str_w | str_width
 * Strand width in pixel
 * Default = 50
* -str_s | str_space
 * Space between strands in pixel
 * Default = 30
* -sp | space_chr
 * Space between chromsomes
 * Default = 50
* -lm | lmargin
 * Left margin in pixel
 * Default = 50
* -rm | rmargin
 * Rigth margin in pixel
 * Default = 50
* -tm | tmargin
 * Top margin in pixel
 * Default = 50
* -bm | bmargin
 * Bottom margin in pixel
 * Default = 50
* -ba | background color
 * Fill Background
 * Default = no Background
* -la | label_strand_rotation
 * Rotation degree of strand label
 * Default = 0
* -ft_f | ft-family
 * Font to use for text
 * Default = Helvetica
* -ft_s | ft-size
 * Size of the font
 * Default = 16


## EXTERNAL LINKS
DensityMap can be used online with a graphical interface : [DensityMap GUI](http://chicken-repeats.inra.fr/launchDM_form.php)

DensityMap as been used to show the satellite (DNA repeats) in chicken genome: [Study using DensityMap](http://chicken-repeats.inra.fr)
