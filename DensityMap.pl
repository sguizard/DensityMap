#!/usr/bin/perl -w

use strict;
use diagnostics;
use warnings;
use Getopt::Long;
#use Term::ANSIColor;
#use Data::Dumper;
use PerlIO::gzip;
use GD::SVG;
use POSIX;
use Cwd 'abs_path';
use File::Basename;

#> Setting Parameters


##> Define outputs colors
#print STDOUT color 'blue';
#print STDERR color 'red';


##> Define options
my %config;
my @gffFiles;
GetOptions (\%config,
            'input=s' => \@gffFiles,
            'region_file=s',
			'fasta=s',
            'type_to_draw=s',
            'output_img_name=s',
            'rounding_method=s',
            'show_scale=i',
            'str_width=i',
            'str_space=i',
            'space_chr=i',
            'scale_factor=i',
            'auto_scale_factor=i',
            'background=s',
            'transparency=s', 
            'lmargin=i',
            'rmargin=i',
            'tmargin=i',
            'bmargin=i',
            'label_strand_rotation=i',
            'title=s',
            'force',
            'win_size=i',
            'colour_scale=i',
            'gc=i',
            'ft_family=s',
            'ft_size=i',
			'verbose',
            'debug',
			'help');


##> Print USAGE if --help
if ($config{help}) {printUsage(1);}


##> Check if gff file exist, if no mandatory parameter are missing
if (!scalar(@gffFiles))               {printError ("\n!!!! input is MANDATORY !!!!! \n\n\n", 0);                 printUsage(1);}
if (!exists $config{output_img_name}) {printError ("\n!!!! output_img_name is MANDATORY !!!!! \n\n\n", 0);       printUsage(1);}
if (!exists $config{region_file})     {printError ("\n!!!! region_file is MANDATORY !!!!! \n\n\n", 0);           printUsage(1);}
if (!exists $config{type_to_draw})    {printError ("\n!!!! type_to_draw options is  MANDATORY !!!!! \n\n\n", 0); printUsage(1);}
#if (! -e $config{input})                  {printError ("gff $config{input} not exist ! \n"); printUsage(1);}

##> Setting Global Variables
my %color;
my %gffTypes;
my %margin;
my %order;
$order{'-'}     = "-";
$order{'+'}     = "+";
$order{'both'}  = "-;+";
$order{'fused'} = "-+";
$order{'all'}   = "-;-+;+";
my %rand;
my %offset;
$offset{'x'}  = 0;
$offset{'y'}  = 0;
$offset{'x1'} = 0;
$offset{'y1'} = 0;
$offset{'x2'} = 0;
$offset{'y2'} = 0;
my %listChr;

$| = 1;

# Files
#my $input_file			  =  $config{input_file};
my $region_file			  =  $config{region_file};
my $fasta_file			  =  $config{fasta};

# Graphic options
my $scale_factor          = ($config{auto_scale_factor})     ? 1 : ($config{scale_factor}) ? $config{scale_factor} : 1000;
my $numMaxTicks           = ($config{show_scale})            ? $config{show_scale}            : 50;
my $strand_width          = ($config{str_width})             ? $config{str_width}             : 50;
my $strand_space          = ($config{str_space})             ? $config{str_space}             : 50;
my $space_chr             = ($config{space_chr})             ? $config{space_chr}             : 50;
   $margin{l}             = ($config{lmargin})               ? $config{lmargin}               : 50;
   $margin{r}             = ($config{rmargin})               ? $config{rmargin}               : 50;
   $margin{t}             = ($config{tmargin})               ? $config{tmargin}               : 50;
   $margin{b}             = ($config{bmargin})               ? $config{bmargin}               : 50;
my $label_strand_rotation = ($config{label_strand_rotation}) ? $config{label_strand_rotation} : 0;
my $colour_scale          = ($config{colour_scale})          ? $config{colour_scale}          : 7;
my $gc_cs                 = ($config{gc})                    ? $config{gc}                    : 7;
my $win_size              = ($config{win_size})              ? $config{win_size}              : 1;
my $rounding_method       = ($config{rounding_method})       ? $config{rounding_method}       : "floor";
my $fts                   = ($config{ft_size})               ? $config{ft_size}               : 16;
my $title                 = ($config{title})                 ? $config{title}                 : 0;
my $v                     = ($config{verbose})               ? 1                              : 0;
my $d                     = ($config{debug})                 ? 1                              : 0;
my $paternType = "(";
my $scaleAddWidth = 100; 
my $numTicks;
my $chr_length;
my $chr_length_reel;
my $picHeight;
my $picWidth;
my $count = 0;
my $countGff = -1;
my $font;



printv("gffs: ".join " ", @gffFiles);
printv("region_file: $region_file");
printv("output_img_name: $config{output_img_name}\n");

##> Setting parameters

# -1. Load region if defined
printv("Reading region file ...");

my %region;
my @chrOrder;
if ($region_file) {
	my $fh_region_file = openr($region_file);
	
	while (<$fh_region_file>) {
		chomp;
		my @bed = split /\t/;
		
		push @chrOrder, $bed[0];
		
		$region{$bed[0]}{length} = $bed[2] - $bed[1];
		$region{$bed[0]}{start}  = $bed[1];
		$region{$bed[0]}{end}    = $bed[2];
	}
	close $fh_region_file;
}

# 0. Read fasta file
if ($fasta_file){
	printv("Reading fasta file ...");
	my $seqHeader;
	
	my $fh_fasta_file = openr($fasta_file);
	
	while ($fh_fasta_file){
		chomp;
	    if (/>(.+)/) {$seqHeader = $1;}
		else         {$region{$seqHeader}{seq} .= $_;}
	}
	close $fh_fasta_file;
}

# 1. Get Image size
## 1.1 Height
printd("Searching Max Sequence Length ... ");

my $numOfGff          = scalar(@chrOrder);
my $maxSequenceLength = 0;

foreach my $k (keys %region){
    $maxSequenceLength = ($region{$k}{length} > $maxSequenceLength) ? $region{$k}{length} : $maxSequenceLength;
}


if ($title) {$margin{'t'} += 40};


if ($config{auto_scale_factor}) {
    while (1) {
        $picHeight = $margin{'t'}
                   + $margin{'b'}
                   + (floor(($maxSequenceLength/$scale_factor)*$win_size));
        last if ($picHeight < $config{auto_scale_factor});
        $scale_factor *= 10;
    }
    print "Selected Scale Factor = $scale_factor\n";
}
else{
    $picHeight = $margin{'t'}
               + $margin{'b'}
               + (floor(($maxSequenceLength/$scale_factor)*$win_size));
}


## 1.2 Width (+Check type option)
my $numOfStrand = 0;
my @type_array;
my %type_valid;

foreach my $type_group (split(/ /, $config{'type_to_draw'})){
	my $typeOption   = "";
	my $keyOption    = "";
	my $valOption    = "";
	my $strandOption = "";
	my $csOption     = "";
	my $roOption     = "";
	
	foreach my $option (split(/&/, $type_group)){
		$option =~ /([^=]+)=([^=]+)/;
		
		my $key = lc($1);
		my $val = $2;
		
		if    ($key eq "type")  {$typeOption   = $val;}
		elsif ($key eq "key")   {$keyOption    = $val;}
		elsif ($key eq "val")   {$valOption    = $val;}
		elsif ($key eq "strand"){$strandOption = $val;}
		elsif ($key eq "cs")    {$csOption     = $val;}
		elsif ($key eq "ro")    {$roOption     = $val;}
		else  {printError("Type option : INVALID KEY, check help please.\n", 1)}
	}
	
	if (!($typeOption ^ ($keyOption && $valOption))){
		printError("Type option : INVALID FORMAT, check help please.\n", 1);
	}
	
	my $type;
	if    ($typeOption){$type = $typeOption;}
	elsif ($keyOption) {$type = "$keyOption=$valOption";}
	
	push @type_array, $type;
	$type_valid{$type}++;
	
	$gffTypes{$type}{strand}   = $strandOption;
	$gffTypes{$type}{colour}   = ($csOption) ? $csOption : $colour_scale;
	$gffTypes{$type}{rounding} = ($roOption) ? $roOption : $rounding_method;
	
	if    ($strandOption eq "-" or
		   $strandOption eq "+" or
		   $strandOption eq "fused"){$numOfStrand++;}
	elsif ($strandOption eq "both") {$numOfStrand += 2;}
	elsif ($strandOption eq "all")  {$numOfStrand += 3;}
	
}

$margin{'l'} = $margin{'l'} + $scaleAddWidth if ($config{show_scale});
$picWidth = $margin{'l'}
          + $margin{'r'}
          + (($numOfGff * $numOfStrand      ) * $strand_width)
          + (($numOfGff * ($numOfStrand - 1)) * $strand_space)
          + (($numOfGff - 1) * $space_chr);

if ($config{gc}){
	$picWidth += ($numOfGff * $strand_width) 
			  +  ($numOfGff * $strand_space);
}

## 1.3 Ask if image size is OK
if (!$config{force}) {
    while (1) {
        print "Is the picture size will be ok ? ($picWidth px x $picHeight px) [y/n] : ";
        chomp (my $response = <STDIN>);
        if    ($response eq "n") {exit;}
        elsif ($response eq "y") {last;}
    }
}

# 2 Create picture
printv("Create Picture ...");
my $image = GD::SVG::Image->new($picWidth, $picHeight);


# 3 Loading colors from colours.txt
printv("Load colors ...");
open(COLOR, "<".dirname(abs_path($0))."/colours.txt") or die "Can not open colours.txt";
while (<COLOR>) {
    next if /^#/;
    /([\d\w]+);(\d+);(\d+);(\d+)/;
    $color{lc($1)} = $image->colorAllocate($2,$3,$4);
}
close COLOR;

# 4 Draw Background
if ($config{'background'}) {
    printv("Draw background ...");
    $image->filledRectangle(0, 0, $picWidth, $picHeight, $color{$config{'background'}});
}


# 5 Drawing Title and Scale
if ($config{title}){
	printv("Drawing title ...");
    if ($config{show_scale}) {
        $image->string( gdLargeFont,
                        $scaleAddWidth + ($picWidth - $scaleAddWidth) / 2 - ((gdLargeFont->width * length $config{title})/2), 
                        20 - ((gdLargeFont->height) / 2),
                        $config{title},
                        $color{'black'});
    }
    else {
        $image->string( gdLargeFont,
                        $picWidth / 2 - ((gdLargeFont->width * length $config{title})/2), 
                        20 - ((gdLargeFont->height) / 2),
                        $config{title},
                        $color{'black'});
    }
}

if ($config{show_scale}) {
    printv("Drawing Scale ...");
    my $scale_size = floor($maxSequenceLength/$scale_factor);
    printd("Scale_size : $scale_size");
    drawScale(\$image, $scale_size, $maxSequenceLength, $scaleAddWidth, $scale_factor, $numMaxTicks, $margin{'t'}, $win_size);
}


# 6 Get Type/Strand
#printv("Parsing types ...");
#foreach (split(/;/, $config{'type_to_draw'})){
#    $_=~ /([^=]+)=([^=]+)=?(\d*)/;
#
#    printd("type = $1\tstrand = $2");
#	
#	$gffTypes{$1}{colour} = ($3) ? $3 : $colour_scale;
#    $gffTypes{$1}{strand} = $2;
#
#}
#printd();



# 7.-1 Load gff files
printv("Reading GFF files ...");

# An annotation will be stored only if it's in the region to plot (describe region bed file)
# 6 cases are possible: 
#   #1 Annot start and end are before region's start
#   #2 Annot start is before region's start, Annot end is in the region
#   #3 Annot start is before region's start, Annot end is after region end
#   #4 Annot start is in the region, Annot end is in the region
#   #5 Annot start is in the region, Annot end  is after region end
#   #6 Annot start and end are after region's end

#######################################################
################### Possibles Cases ################### x <- annot start 
####################################################### + <- annot end 
#          Start                         End          # <-<-<- Region
#######################################################
#1   x   +   |                            ]           # 
#2     x     |              +             ]           # 
#3     x     |                            |     +     # 
#4           |       x             +      ]           # 
#5           |              x             ]     +     # 
#6           |                            ]   x   +   # 
#######################################################

my %gffData;
my $types_pattern = join " ", @type_array;
foreach my $gffFile (@gffFiles){
	my $fh_g = openr($gffFile);
	while (<$fh_g>){
		last if /^##FASTA/;
		next if /^#/;
		
		chomp;
		
		my @line = split /\t/;
		my $chr       = $line[0];
		my $type      = $line[2];
		my $start     = $line[3];
		my $end       = $line[4];
		my $strand    = $line[6];
		my $attributs = $line[8];
		
		next if (!exists($type_valid{$type}));
		foreach (grep {/=/} keys(%type_valid)){
			next if ($attributs !~ /$_/)
		}
		
		if (!exists($gffData{$chr}{$type})){
			$gffData{$chr}{$type} = [];
		}
		else {
			if ($start < $region{$chr}{start}){ # 1 or 2 or 3
				if ($end < $region{$chr}{start}){ # 1
					next;
				}
				elsif ($end <= $region{$chr}{end}){ # 2
					push @{$gffData{$chr}{$type}}, [$strand, $region{$chr}{start}, $end];
				}
				elsif ($end > $region{$chr}{end}){ # 3
					push @{$gffData{$chr}{$type}}, [$strand, $region{$chr}{start}, $region{$chr}{end}];
				}
			}
			elsif ($start < $region{$chr}{end}){ # 4 or 5
				if ($end <= $region{$chr}{end}){
					push @{$gffData{$chr}{$type}}, [$strand, $start, $end]; # 4
				}
				elsif ($end > $region{$chr}{end}){
					push @{$gffData{$chr}{$type}}, [$strand, $start, $region{$chr}{end}]; # 5
				}
			}
			else {next;} # 6
		}
	}
	close $fh_g;
}

# 7.0 Foreach file draw strand(s)
my %centromere;

my $csv = $config{output_img_name};
$csv =~ s/.svg/.csv/g;
open(CSV, ">$csv") or die "Can not open $csv ! ";
print CSV "sequence\tfeature\tstart\tend\tdensity\n";


## TODO: reimplmente centromere support
printv("Ploting process ...");
foreach my $chr (@chrOrder) {
	printv("==> Processing $chr ...");
	
	processData($chr);
	$countGff++;
}

close CSV;


# 8 Applying rotation to labels and Saving picture
open(IMG, ">$config{output_img_name}") or die "Can not open $config{output_img_name} ! ";
binmode IMG;
my @text = split("\t", $image->svg);
#if (!$config{transparency}) {
    my $switchRotate = 0;
    foreach (@text){
        if    ($config{label_strand_rotation} and /<g id="rotate_\d+">/)     {$switchRotate++;}
        elsif ($config{label_strand_rotation} and /<\/g>/ and $switchRotate) {$switchRotate--;}
        
        if ($switchRotate) {s/x="([\d\.]+)" y="([\d\.]+)"/x="$1" y="$2" transform="rotate($label_strand_rotation, $1, $2)"/;}
		
        if ($config{ft_family}) {s/font="Helvetica"/font-family="$config{ft_family}"/;}
        if ($fts != 16)         {s/font-size="16"/font-size="$fts"/;}
        #s/stroke-opacity: 1.0; stroke-width: 1/stroke: none/g;
        #s/stroke-opacity: 1.0;/stroke: none;/g;
        s/stroke-opacity: 1.0;?//g;
        s/stroke-width: 1;?//g;
        s/stroke: rgb\(\d+,\d+,\d+\);?//g;
        s/ +/ /g;
        s/style=" /style="/g;
        s/;? " width/" width/g;
        print IMG $_;
    }    
#}

#WAS USED FOR ADDING TRANSPARENCY ON MERGED TRACKS
#else {
#    my $switch = 0;
#    my $switchRotate = 0;
#    
#    foreach (@text){
#        #s/stroke\-width\: 1/stroke\-width\: 0/g;
#        s/stroke-opacity: 1.0; stroke-width: 1/stroke: none/g;
#        if (/<g id="\-\+_\d+">/)    {$switch++;}
#        elsif (/<\/g>/ and $switch) {$switch--;}
#        
#        if (/<g id="rotate_\d+">/)         {$switchRotate++;}
#        elsif (/<\/g>/ and $switchRotate)  {$switchRotate--;}
#        
#        if ($switch) {
#            s/fill-opacity: 1.0/fill-opacity: $config{transparency}/g;
#            s/stroke-opacity: 1.0; stroke-width: 1/stroke: none/g;
#        }
#        if ($switchRotate) {s/x="(\d+)" y="(\d+)"/x="$1" y="$2" transform="rotate(-45, $1, $2)"/;}
#        print IMG $_;
#    }
#}

close IMG;
print "Image saved ! \n";


###########################################################################
################################ Fonctions ################################
###########################################################################
sub processData{
	my $chr = shift;
	my @minus;
	my @plus;
	my @minusPlus;
	
	printd("processData: ");
	
	foreach my $type (@type_array) {
		printv("====> Processing type $type ...");
		@minus     = grep {$_->[0] eq "-"} @{$gffData{$chr}{$type}};
		@plus      = grep {$_->[0] eq "+"} @{$gffData{$chr}{$type}};
		@minusPlus =                       @{$gffData{$chr}{$type}};
		
		@minus     = sort {$a->[1] <=> $b->[1]} @minus;
        @plus      = sort {$a->[1] <=> $b->[1]} @plus;
        @minusPlus = sort {$a->[1] <=> $b->[1]} @minusPlus;
		
		@minus     = removeIntervalRedundancy(@minus);
        @plus      = removeIntervalRedundancy(@plus);
        @minusPlus = removeIntervalRedundancy(@minusPlus);
		
	    # 7.4 Set Offset
        if ($gffTypes{$type}{'strand'} eq "-")     {
            $offset{$type}{'-'}{'x'} = $margin{'l'}
			                      + ($count * $strand_width)
								  + ($count * $strand_space)
								  - ($countGff * $strand_space)
								  + ($countGff * $space_chr);
            $offset{$type}{'-'}{'y'} = $margin{'t'};
            $count++;
            
            printd("processData: x = $offset{$type}{'-'}{'x'}\ny = $offset{$type}{'-'}{'y'}");
        }
        elsif ($gffTypes{$type}{'strand'} eq "+")     {
            $offset{$type}{'+'}{'x'} = $margin{'l'}
			                      + ($count * $strand_width)
								  + ($count * $strand_space)
								  - ($countGff * $strand_space)
								  + ($countGff * $space_chr);
            $offset{$type}{'+'}{'y'} = $margin{'t'};
            $count++;
            
            printd("processData: + = $offset{$type}{'+'}{'x'}\n+ = $offset{$type}{'+'}{'y'}");
        }
        elsif ($gffTypes{$type}{'strand'} eq "both")  {
            $offset{$type}{'-'}{'x'} = $margin{'l'}
			                      + ($count * $strand_width)
								  + ($count * $strand_space)
								  - ($countGff * $strand_space)
								  + ($countGff * $space_chr);
            $offset{$type}{'-'}{'y'} = $margin{'t'};
            $count++;
            
            $offset{$type}{'+'}{'x'} = $margin{'l'}
			                      + ($count * $strand_width)
								  + ($count * $strand_space)
								  - ($countGff * $strand_space)
								  + ($countGff * $space_chr);
            $offset{$type}{'+'}{'y'} = $margin{'t'};
            $count++;
            
            printd("processData: + = $offset{$type}{'-'}{'x'}\n+ = $offset{$type}{'-'}{'y'}\n- = $offset{$type}{'+'}{'x'}\n- = $offset{$type}{'+'}{'y'}");
        }
        elsif ($gffTypes{$type}{'strand'} eq "fused") {
            $offset{$type}{'-+'}{'x'} = $margin{'l'}
								   + ($count * $strand_width)
								   + ($count * $strand_space)
								   - ($countGff * $strand_space)
								   + ($countGff * $space_chr);
            $offset{$type}{'-+'}{'y'} = $margin{'t'};
            $count++;
            
            printd("processData: f = $offset{$type}{'-+'}{'x'}\nf = $offset{$type}{'-+'}{'y'}");
        }
        elsif ($gffTypes{$type}{'strand'} eq "all")   {
            $offset{$type}{'-'}{'x'} = $margin{'l'}
								  + ($count * $strand_width)
								  + ($count * $strand_space)
								  - ($countGff * $strand_space)
								  + ($countGff * $space_chr);
            $offset{$type}{'-'}{'y'} = $margin{'t'};
            $count++;
            
            $offset{$type}{'-+'}{'x'} = $margin{'l'}
								   + ($count * $strand_width)
								   + ($count * $strand_space)
								   - ($countGff * $strand_space)
								   + ($countGff * $space_chr);
            $offset{$type}{'-+'}{'y'} = $margin{'t'};
            $count++;
            
            $offset{$type}{'+'}{'x'} = $margin{'l'}
								  + ($count * $strand_width)
								  + ($count * $strand_space)
								  - ($countGff * $strand_space)
								  + ($countGff * $space_chr);
            $offset{$type}{'+'}{'y'} = $margin{'t'};
            $count++;
            
            printd("processData: - = $offset{$type}{'-'}{'x'}\n- = $offset{$type}{'-'}{'y'}\nf = $offset{$type}{'-+'}{'x'}\nf = $offset{$type}{'-+'}{'y'}\n+ = $offset{$type}{'+'}{'x'}\n+ = $offset{$type}{'+'}{'y'}");
        }
		
		printd("processData: typeToDraw = $type");
    
        foreach my $strand (split(";", $order{$gffTypes{$type}{'strand'}})){
            printd("processData: strand = $strand");
            
            my $cs = $gffTypes{$type}{colour};
			my $ro = $gffTypes{$type}{rounding};
			my $ref_tab;
			if    ($strand eq "-") {$ref_tab =\@minus;}
			elsif ($strand eq "+") {$ref_tab =\@plus;}
			elsif ($strand eq "-+"){$ref_tab =\@minusPlus;}
			
            drawPixels(\$image, \%rand, $cs, $chr, $region{$chr}{length}, $scale_factor,
					   $type, $strand, $strand, \%centromere, $ref_tab, $win_size, $ro);
        }
    }
    
    if ($config{gc}){
        $offset{'gc'}{'x'} = $margin{'l'}
		                   + ($count * $strand_width)
						   + ($count * $strand_space)
						   - ($countGff * $strand_space)
						   + ($countGff * $space_chr);
        $offset{'gc'}{'y'} = $margin{'t'};
        $count++;
		
		printError("Fasta sequence is not in the gff file ! ", 1) if (!$region{$chr}{seq});
        drawPixelsGC(\$image, \%rand, $gc_cs, $chr, $region{$chr}{length}, $scale_factor, $win_size, \$region{$chr}{seq});
    }
}

###########################################################################
sub removeIntervalRedundancy{
    
    # The purpose of this function is to remove or merging the feature of the gff
    # Input  : sorted array of start - end array
    # Output : sorted array of start - end array without crossing intervals
    
    #my ($ref_intervalsFused, @gff) = @_;
    my (@gff) = @_;
    my @reduced_gff;
    my $bool_firstInterval = 1;
    
	printd("removeIntervalRedundancy: ");
	
    printd("removeIntervalRedundancy: Start removeIntervalRedundancy ...");
    printd("removeIntervalRedundancy: \tStarting gff size = ".@gff);
    
    
    foreach my $annot (@gff) {
        
        # The intervals of the gff have been sorted by 'start' before using this function
        # so only three case are possebles :
        #   #1 Previous and current intervals crossing
        #   #2 Current interval included in previous interval
        #   #3 Current interval not inclued in previous interval
        
        #######################################################
        ################### Possibles Cases ###################
        ####################################################### x + <- Current start/end interval
        #          Start                         End          # <-<-<- Previous interval
        #######################################################
        #1           |            x               ]     +     # Intervales chevauchants
        #2           |            x         +     ]           # Intervales inclu dans le précédent
        #3           |                            |   x    +  # Intervales non-chevauchants
        #######################################################
        
        print "\t---> Current interval = $annot->[1], $annot->[2]\n"            if $config{insaneDebugMode};
        
        #Load first interval
        if ($bool_firstInterval) {
            print "\nLoad first interval\n"                                     if $config{'insaneDebugMode'};
            
            push(@reduced_gff, [$annot->[1], $annot->[2]]);
            $bool_firstInterval--;
            
            print "\tannot = $annot\n"                                          if $config{insaneDebugMode};
            print "\t\$annot->[1], \$annot->[2] = $annot->[1], $annot->[2]\n"   if $config{insaneDebugMode}; 
            print "\t\@reduced_gff size = ".@reduced_gff."\n\n"                 if $config{insaneDebugMode};
            
            next;
        }
        
        if ($annot->[1] <= $reduced_gff[$#reduced_gff]->[1] and     #1
            $annot->[2] >  $reduced_gff[$#reduced_gff]->[1]) {
            # Replace previous end
            $reduced_gff[$#reduced_gff]->[1] = $annot->[2];
        }
        elsif ($annot->[1] <= $reduced_gff[$#reduced_gff]->[1] and  #2
               $annot->[2] <=  $reduced_gff[$#reduced_gff]->[1]) {
            # Go to next interval
            next;
        }
        else {                                                      #3
            # Insert interval in reduced gff
            push(@reduced_gff, [$annot->[1], $annot->[2]]);
        }
    } #End foreach (@gff)
    
    printd("removeIntervalRedundancy: \tFinal gff size    = ".@reduced_gff);
    printd("removeIntervalRedundancy: Finished\n");

    return @reduced_gff;
}

############################################################################
sub drawScale{
    
    # Drawing the scale on the left side of graph
    # Input :
    #   - $ref_img          -> reference of the image
    #   - $chr_size         -> chromosome/sequence size divided by scale_factor
    #   - $chr_size_reel    -> chromosome/sequence size
    #   - $widthScale       -> Width of the scale
    #   - $basesPerPixel    -> Number of bases per pixels (scale_factor)
    #   - $maxTicks         -> Number Max of ticks
    #   - $marginTop        -> size of top margin
    # Output : none
    
    my ($ref_img, $chr_size, $chr_size_reel, $widthScale, $basesPerPixel, $maxTicks, $marginTop, $win_size) = @_;
    my $basesPerTicks = 10;

    # Search the number of tick to use
    printd("drawScale: Start searching num ticks");
    while (1) {
        $numTicks = floor($chr_size/$basesPerTicks);
        printd("basesPerTicks = $basesPerTicks \t numTicks = $numTicks");
        last if ($numTicks <= $maxTicks);
        $basesPerTicks*=10; 
    }
    printd("drawScale: Found num ticks");
    
    # Define the 10^x bases to use as unit
    my $power = floor(log10($basesPerTicks*$basesPerPixel));
    
    # Add label unit
    my $string = "(x10e$power bases)";
    $$ref_img->string(gdSmallFont,
                   $widthScale - 15 - (gdSmallFont->width * length $string),
                   $marginTop - 5,
                   $string,
                   $color{'black'});
    
    #Add ratio pixel base label (2 lines)
    $string = "(1 win = ";
    $$ref_img->string(gdSmallFont,
                   $widthScale - 15 - (gdSmallFont->width * length "$basesPerPixel bases)"),
                   $marginTop - 5 + gdSmallFont->height,
                   $string,
                   $color{'black'});
    
    $string = "$basesPerPixel bases)";
    $$ref_img->string(gdSmallFont,
                   $widthScale - 15 - (gdSmallFont->width * length $string),
                   $marginTop - 5 + (gdSmallFont->height * 2),
                   $string,
                   $color{'black'});
    
    # Draw scale
    printd("drawScale: \$\$ref_img->filledRectangle($widthScale - 10, $marginTop, $widthScale - 8, $marginTop + $chr_size, \$color{'black'})\;");
    $$ref_img->filledRectangle($widthScale - 10,
                            $marginTop,
                            $widthScale - 8,
                            $marginTop + $chr_size * $win_size,
                            $color{'black'});
    
    # Draw scale last tick
    $string = sprintf("%.2f", $chr_size_reel/(10**$power));
    $$ref_img->string(gdLargeFont,
                   #$widthScale - 15 - (gdLargeFont->width * length $string),
                   $widthScale - 15,
                   #$marginTop + $chr_size - 8,
                   $marginTop + $chr_size * $win_size,
                   $string,
                   $color{'black'});
    printd("drawScale: \$\$ref_img->filledRectangle($widthScale - 10, $marginTop + $chr_size - 1, $widthScale, $marginTop + $chr_size, \$color{'black'})\;");
    $$ref_img->filledRectangle($widthScale - 10,
                            $marginTop + $chr_size * $win_size - 1,
                            $widthScale,
                            $marginTop + $chr_size * $win_size,
                            $color{'black'});
    
    #Draw scale other ticks
    for (my $currentTick = 0 ; $currentTick <= $numTicks ; $currentTick++) {
        $string = sprintf("%.2f", ((($basesPerTicks * $currentTick) * $basesPerPixel) /10**$power));
        $string = 0 if ($currentTick == 0);
        $$ref_img->string(gdLargeFont,
                       $widthScale - 15 -(gdLargeFont->width * length $string),
                       $marginTop + ($basesPerTicks * $currentTick) * $win_size - 8,
                       $string,
                       $color{'black'});
         
        printd("drawScale: \$\$ref_img->filledRectangle($widthScale - 10, $marginTop + ($basesPerTicks * $currentTick), $widthScale, $marginTop + ($basesPerTicks * $currentTick) + 1, \$color{'black'})\;");
        $$ref_img->filledRectangle($widthScale - 10,
                                $marginTop + ($basesPerTicks * $currentTick) * $win_size,
                                $widthScale,
                                $marginTop + ($basesPerTicks * $currentTick) * $win_size + 1,
                                $color{'black'});
    }
	printd("drawScale: ");
}

###########################################################################
sub drawPixels{
    
    # Draw each pixel of strand
    # Input :
    #   - $ref_img      ->  ref of the image
    #   - $ref_rand     ->  ref on the hash of random numbers
    #   - $ref_colour_scale  ->  colour_scale to use for colors
    #   - $chr_size     ->  Chromosome/Sequence size
    #   - $scaleFactor  ->  Scale_factor
    #   - $type         ->  current to draw (label)
    #   - $strand       ->  current strand in process (label)
    #   - $strandColor  ->  color to use with strand
    #   - $ref_gff      ->  ref of the current strand gff
    # Output: none
    
    printv("======> Start Drawing pixels ...");
	printd("drawPixels: ");
    
    my ($ref_img, $ref_rand, $cs, $seqName,
		$chr_size, $scaleFactor, $type, $strand,
		$strandColor, $ref_centromere, $ref_gff, $win_size, $ro) = @_;
    my @gff    = @{$ref_gff};
    my %centro = %{$ref_centromere};
    my $randNum;
    my %intervals;
    my @previousBases = (0);
    
    # Search a unique random number
    while (1) {
        $randNum = int rand(1000);
        redo if $$ref_rand{$randNum}++;
        last;
    }
    
    printd("drawPixels: chr_size    = $chr_size");
    printd("drawPixels: scaleFactor = $scaleFactor");
    printd("drawPixels: type        = $type");
    printd("drawPixels: strand      = $strand\n");
    
    # Open chromosome/sequence group
    $$ref_img->startGroup("${strand}_${randNum}");
    
    # For each pixel of the chromosome/sequence
    my $posPic = 0;
    for (my $pos = 0 ; $pos <= $chr_size ; $pos++) {
        
    	next if ($config{region_file} && (($pos*$scaleFactor) < $region{$seqName}{start} || ($pos*$scaleFactor+$scaleFactor) > $region{$seqName}{end}));
		$posPic++;
        
		# Get number of base covered by the previous on pixel crosssing interval
        my $basesCoverred = shift @previousBases;
        $basesCoverred = 0 if (!defined $basesCoverred); # if @previousBases is empty
        
        printd("drawPixels: pos = $pos");
        printd("drawPixels: basesCoverred Tab = -".join(" ", @previousBases)."-");
        printd("drawPixels: basesCoverred = $basesCoverred (previous shift)");
        printd("drawPixels: Start while");
        printd("drawPixels: \$gff[0]->[0] = $gff[0]->[0]");
		printd("drawPixels: \$pos * \$scaleFactor = ".($pos * $scaleFactor));
        #printError( "\$gff[0]->[0] undef ! \n") if(!defined($gff[0]->[0]));
        
        # while the end of gff is not reached and the next start is in the current pixel
        while (defined($gff[0]->[0]) and $gff[0]->[0] < ($pos * $scaleFactor)) {
            my $ref_interval = shift(@gff);
            
            last if (!defined $ref_interval->[0]); # = defined($gff[0]->[0]) avoid warnings
            
            # get clearly start and end 
            $intervals{'start_reel'}  = $ref_interval->[0];
            $intervals{'end_reel'}    = $ref_interval->[1];

            printd("drawPixels: \tstart_reel = $intervals{'start_reel'}");
            printd("drawPixels: \tend_reel   = $intervals{'end_reel'}");
            printd("drawPixels: \tcheck position");
            
            # feature end is on current pixel
            if ($intervals{'end_reel'} < ($pos * $scaleFactor) ){ 
                printd("drawPixels: \tCase 1 : \$basesCoverred += $intervals{'end_reel'} - $intervals{'start_reel'}");
                printd("drawPixels: \tCase 1 : \$basesCoverred += ".($intervals{'end_reel'} - $intervals{'start_reel'}));
                $basesCoverred += $intervals{'end_reel'} - $intervals{'start_reel'}; # increment the base counting
            }
            
            # feature end is on next pixels
            else{
                # determine on how much pixel the interval is crossing
                my $numPixelImpliy = floor(($intervals{'end_reel'} - $intervals{'start_reel'} - (($pos * $scaleFactor) - $intervals{'start_reel'}))/$scaleFactor);
                
                # Treating current position
                $basesCoverred += ($pos * $scaleFactor) - $intervals{'start_reel'}; 
                printd("drawPixels: \tCase 2 : \$basesCoverred += ($pos * $scaleFactor) - $intervals{'start_reel'}");
                
                # Treating the next full pixels 
                for (my $i = 0 ; $i < $numPixelImpliy ; $i++){
                    printd("drawPixels: \tCase 2 : \$previousBases[$i] += $scaleFactor");
                    $previousBases[$i] += $scaleFactor;
                }
                
                # Treat the last not full pixel
                $previousBases[$numPixelImpliy] += (($intervals{'end_reel'} - $intervals{'start_reel'} - (($pos * $scaleFactor) - $intervals{'start_reel'}))%$scaleFactor);
                printd("drawPixels: \tCase 2 : \$previousBases[$numPixelImpliy] += (($intervals{'end_reel'} - $intervals{'start_reel'} - (($pos * $scaleFactor) - $intervals{'start_reel'}))%$scaleFactor);");
                
            }# feature end is on next pixels
            
            printd("drawPixels: \tbasesCoverred = $basesCoverred");
            
        }# while ($intervalsTabM[0]->[0] < ($pos * $scaleFactor))
        
        printd("drawPixels: Stop while");
        
        # Compute percentage of base coverage
        my $percentage;
        if    ($ro eq "floor") {$percentage = floor(($basesCoverred/$scaleFactor) * 100);}
        elsif ($ro eq "ceil" ) {$percentage = ceil (($basesCoverred/$scaleFactor) * 100);}
        
        printd("drawPixels: percentage = $basesCoverred/$scaleFactor = ".($basesCoverred/$scaleFactor));
        printd("drawPixels: percentage = $percentage % \n");
        print "$pos = ${cs}_heatmap$percentage\n" if $config{insaneDebugMode};
        
        # kill if more than 100 % 
        if    ($percentage > 100) {printError("Higher than 100 % ($percentage % )", 1);}
		elsif ($percentage < 0)   {printError("Lesser than 0 % ($percentage % )", 1);}
        
        # Draw the current pixel
	if ($config{region_file}) {
		$$ref_img->filledRectangle($offset{$type}{$strand}{'x'},                 $offset{$type}{$strand}{'y'} + $posPic * $win_size,
        	                       $offset{$type}{$strand}{'x'} + $strand_width, $offset{$type}{$strand}{'y'} + $posPic * $win_size + $win_size,
        	                       $color{"${cs}_heatmap$percentage"});
	}
	else {
		$$ref_img->filledRectangle($offset{$type}{$strand}{'x'},                 $offset{$type}{$strand}{'y'} + $pos * $win_size,
        	                       $offset{$type}{$strand}{'x'} + $strand_width, $offset{$type}{$strand}{'y'} + $pos * $win_size + $win_size,
        	                       $color{"${cs}_heatmap$percentage"});
	}

	my $st = $pos * $scaleFactor;
	my $en = $pos * $scaleFactor + $scaleFactor;
	print CSV "$seqName\t$type\t$st\t$en\t$percentage\n";
	
    }# End # For each pixel of the chromosome/sequence
    
    # Draw Centromere
    if (defined($centro{start})) {
        #Left Triangle
        my $poly = new GD::Polygon;
        $poly->addPt($offset{$type}{$strand}{'x'} - 1, ($offset{$type}{$strand}{'y'} + $$ref_centromere{start}/$scaleFactor));
        $poly->addPt($offset{$type}{$strand}{'x'} - 1, ($offset{$type}{$strand}{'y'} + $$ref_centromere{end}  /$scaleFactor));
        $poly->addPt(($offset{$type}{$strand}{'x'} + $strand_width/2), ($offset{$type}{$strand}{'y'} + ($$ref_centromere{start}/$scaleFactor + $$ref_centromere{end}/$scaleFactor)/2));
        $$ref_img->filledPolygon($poly, $color{$config{'background'}});
        
        #Right Triangle
        $poly = new GD::Polygon;
        $poly->addPt($offset{$type}{$strand}{'x'} + $strand_width + 1, ($offset{$type}{$strand}{'y'} + $$ref_centromere{start}/$scaleFactor));
        $poly->addPt($offset{$type}{$strand}{'x'} + $strand_width + 1, ($offset{$type}{$strand}{'y'} + $$ref_centromere{end}  /$scaleFactor));
        $poly->addPt(($offset{$type}{$strand}{'x'} + $strand_width/2), ($offset{$type}{$strand}{'y'} + ($$ref_centromere{start}/$scaleFactor + $$ref_centromere{end}/$scaleFactor)/2));
        $$ref_img->filledPolygon($poly, $color{$config{'background'}});
    }
    
    # Open label group for next rotation 
    $$ref_img->startGroup("rotate_$randNum");
    
    # Draw label Sequence Name
    $$ref_img->string(gdLargeFont,
                   $offset{$type}{$strand}{'x'} + (($strand_width/2) - (length($seqName) * gdLargeFont->width)/2), 
                   $offset{$type}{$strand}{'y'} - (3 * gdLargeFont->height), 
                   $seqName,
                   $color{'black'});
    
    # Draw label type
    $$ref_img->string(gdLargeFont,
                   $offset{$type}{$strand}{'x'} + (($strand_width/2) - (length($type) * gdLargeFont->width)/2), 
                   $offset{$type}{$strand}{'y'} - (2 * gdLargeFont->height), 
                   $type,
                   $color{'black'});
    # Close label group for next rotation
    $$ref_img->endGroup;
    
    # Draw strand label
    $$ref_img->string(gdLargeFont,
                   $offset{$type}{$strand}{'x'} + (($strand_width/2) - (length($strand) * gdLargeFont->width)/2), 
                   $offset{$type}{$strand}{'y'} - gdLargeFont->height, 
                   $strand,
                   $color{'black'});
    
    # close chromosome/sequence group
    $$ref_img->endGroup;
    printv("======>  Finish Drawing pixels.");
}

###########################################################################
sub drawPixelsGC{
    
    # Draw each pixel of strand
    # Input :
    #   - $ref_img      ->  ref of the image
    #   - $ref_rand     ->  ref on the hash of random numbers
    #   - $cs           ->  colour_scale to use for colors
    #	- $seqName		->	Name of the annotated sequence
    #   - $chr_size     ->  Chromosome/Sequence size
    #   - $scaleFactor  ->  Scale_factor
    #	- $win_size		->	Size of the printed window in pixel
    #	- $ref_sequence	->	sequence of the chromosome
    # Ouput : none
    
    print "Start Drawing pixels GC ...\n" if $config{'verbose'};
    
    my ($ref_img, $ref_rand, $cs, $seqName, $chr_size, $scaleFactor, $win_size, $ref_sequence) = @_;
    my $randNum;
    
    # Search a unique random number
    while (1) {
        $randNum = int rand(1000);
        redo if $$ref_rand{$randNum}++;
        last;
    }
    
    print "\tchr_size    = $chr_size\n"       if $config{'debug'};
    print "\tscaleFactor = $scaleFactor\n"    if $config{'debug'};
    
    # Open chromosome/sequence group
    $$ref_img->startGroup("+-_${randNum}");
    
    # For each pixel of the chromosome/sequence
    my $posPic = 0;
    for (my $pos = 0 ; $pos <= $chr_size ; $pos++) {
	#Get next window sequence
	my $win_seq = substr $$ref_sequence, 0, $scaleFactor, '';
	
	next if ($config{region_file} && (($pos*$scaleFactor) < $region{$seqName}{start} || ($pos*$scaleFactor+$scaleFactor) > $region{$seqName}{end}));
	$posPic++;

	#Count GC bases
	my $gcBases = $win_seq =~ tr/GC//;

        # Compute percentage of base coverage
        my $percentage;
        if    ($rounding_method eq "floor") {$percentage = floor(($gcBases/$scaleFactor) * 100);}
        elsif ($rounding_method eq "ceil" ) {$percentage = ceil (($gcBases/$scaleFactor) * 100);}
        
        # kill if more than 100 % 
        if ($percentage > 100) {printError("Higher than 100 % ($percentage % )", 1);}
        
        # Draw the current pixel
	if ($config{region_file}) {
        	$$ref_img->filledRectangle($offset{gc}{'x'},                 $offset{gc}{'y'} + $posPic * $win_size,
                	                   $offset{gc}{'x'} + $strand_width, $offset{gc}{'y'} + $posPic * $win_size + $win_size,
                	                   $color{"${cs}_heatmap$percentage"});
	}
	else {
		$$ref_img->filledRectangle($offset{gc}{'x'},                 $offset{gc}{'y'} + $pos * $win_size,
                                           $offset{gc}{'x'} + $strand_width, $offset{gc}{'y'} + $pos * $win_size + $win_size,
                                           $color{"${cs}_heatmap$percentage"});
	}

        my $st = $pos * $scaleFactor;
        my $en = $pos * $scaleFactor + $scaleFactor;

        print CSV "$seqName\tGC%\t$st\t$en\t$percentage\n";

    }# End # For each pixel of the chromosome/sequence
    
    # Open label group for next rotation 
    $$ref_img->startGroup("rotate_$randNum");
    
    # Draw label Sequence Name
    $$ref_img->string(gdLargeFont,
                   $offset{gc}{'x'} + (($strand_width/2) - (length($seqName) * gdLargeFont->width)/2), 
                   $offset{gc}{'y'} - (3 * gdLargeFont->height), 
                   $seqName,
                   $color{'black'});
    
    # Draw label type
    $$ref_img->string(gdLargeFont,
                   $offset{gc}{'x'} + (($strand_width/2) - (length("GC%") * gdLargeFont->width)/2), 
                   $offset{gc}{'y'} - (2 * gdLargeFont->height), 
                   "GC%",
                   $color{'black'});
    # Close label group for next rotation
    $$ref_img->endGroup;
    
    # Draw strand label
    $$ref_img->string(gdLargeFont,
                   $offset{gc}{'x'} + (($strand_width/2) - (length("-+") * gdLargeFont->width)/2),
                   $offset{gc}{'y'} - gdLargeFont->height, 
                   "-+",
                   $color{'black'});
    
    # close chromosome/sequence group
    $$ref_img->endGroup;
    print "Finish Drawing pixels.\n" if $config{'verbose'};
}

###########################################################################
sub openr{
    my $file = shift;
	my $fh;
	
	if ($file =~ /\.gz$/){open $fh, "<:gzip", $file or printError("ERROR : Cannot open $file ! \n", 1);}
    else                 {open $fh, "<",      $file or printError("ERROR : Cannot open $file ! \n", 1);}
	
    return $fh;
}

###########################################################################
sub openw{
    my $file = shift;
    
    open my $fh, ">", $file or printError("ERROR : Cannot open $file ! \n", 1);
    return $fh;
}

###########################################################################
sub sortedKeys{
    my $hashRef = shift;
    
    return sort {$a cmp $b} (keys %$hashRef);
}

###########################################################################
sub printv{
    my $str = shift;
    
    print "$str\n" if $v;
}

###########################################################################
sub printd{
    my $str = shift;
    $str = " " if (!$str);
    
    print STDERR "DEBUG: $str\n" if $d;
}
###########################################################################
sub printError{
    my $string = shift;
    my $exit = shift;
    
    print STDERR $string;
    exit if $exit;
}

###########################################################################
sub printUsage{
    my $exit = shift;
    
    print STDOUT
"USAGE: DensityMap.pl -i chromosome1.gff3 -i chromosome2.gff3 -ty 'match=all' -o Chromosome
    
Options:
    -i     | input                 [string]    Gff file                                (Mandatory)
    -re    | region_file           [string]    A bed file describing the regions to plot 
                                               for each sequence to plot. Allow to plot 
                                               sepecific region and not whole sequence.
                                               Exemple:
                                               seq1\t100000\t200000
                                               seq2\t200000\t350000
    -o     | output_img_name       [string]    Name of the output image                (Mandatory)
    -ty    | type_to_draw          [string]    List of type (column 3 of gff)          (Mandatory)
                                               to draw, strand to use and color scale
                                               Type: match, gene, CDS, ...
                                               Strand: - -        -> strand -
                                                       - +        -> strand +
                                                       - both     -> strand - and strand +
                                                       - fused    -> Combination of strand - and strand +
                                                       - all      -> strand - and strand + and fused
                                               Format: \"match=all;gene=both;CDS=fused\" or
                                                       \"match=all=7;gene=both=7;CDS=fused=10\"

Generic options: 
    -for   | force                 [booleen]   Automaticaly answer yes to picture size validation
    -v     | verbose               [booleen]   MORE text dude !!!!
    -h     | help                  [booleen]   This help

Density options: 
    -c     | colour_scale          [integer]   Color scale to use    (Default = 7)
    -sc    | scale_factor          [integer]   = window length in bp (Default = 1000)
    -a     | auto_scale_factor     [integer]   Max picture height in pixel
    -ro    | rounding_method       [string]    floor or ceil         (Default = floor)
    -gc    | gc                    [integer]   if set, add a density map of the GC%

Graphical options: 
    -ti    | title                 [string]    Title to print on the picture
    -w     | win_size              [integer]   Height of window in pixel       (Default = 1)
    -sh    | show_scale            [integer]   Draw scale, the integer indicate the maximum 
                                               tick to print on the scale      (Default = 50)
    -str_w | str_width             [integer]   Strand width in pixel           (Default = 50)
    -str_s | str_space             [integer]   Space between strands in pixel  (Default = 50)
    -sp    | space_chr             [integer]   Space between chromsomes        (Default = 50)
    -lm    | lmargin               [integer]   Left margin in pixel            (Default = 50)
    -rm    | rmargin               [integer]   Rigth margin in pixel           (Default = 50)
    -tm    | tmargin               [integer]   Top margin in pixel             (Default = 50)
    -bm    | bmargin               [integer]   Bottom margin in pixel          (Default = 50)
    -ba    | background color      [string]    Fill Background                 (Default = no Background)
    -la    | label_strand_rotation [integer]   Rotation degree of strand label (Default = 0)
    -ft_f  | ft-family             [string]    Font to use for text            (Default = Helvetica)
    -ft_s  | ft-size               [integer]   Size of the font                (Default = 16)


    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
\n\n";

#    -d   | debug                     Display debug info
#    -i   | insaneDebugMode           So much prints !!!! 
#    -tr | transparency n            Set transparency, integer between 0 and 1 included

    exit if $exit;
}
