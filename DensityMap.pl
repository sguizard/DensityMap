#!/usr/bin/perl -w

use strict;
use diagnostics;
use warnings;
use Getopt::Long;
#use Term::ANSIColor;
#use Data::Dumper;
use GD::SVG;
use POSIX;

#> Setting Parameters


##> Define outputs colors
#print STDOUT color 'blue';
#print STDERR color 'red';


##> Define options
my %config;
GetOptions (\%config,
            'gff=s',
            'type_to_draw=s', 
            'output_img_name=s',
            'rounding_method=s',
            'show_scale=i',
            'strand_width=i',
            'space_between_str=i',
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
            'yes',
            'win_size=f',
            'verbose',
            'debug',
            'insaneDebugMode',
            'colour_scale=i');


##> Print USAGE if --help
if ($config{help}) {printUsage(1);}


##> Check if gff file exist, if no mandatory parameter are missing
if (!exists $config{gff} or
    !exists $config{output_img_name} or
    !exists $config{type_to_draw})
                                        {printError ("\n!!!!gff, output_img_name, type_to_draw options are MANDATORY !!!!! \n\n\n", 0); printUsage(1);}
#if (! -e $config{gff})                  {printError ("gff $config{gff} not exist ! \n"); printUsage(1);}

##> Setting Global Variables
my $paternType;
my $scaleAddWidth = 100; 
my $numMaxTicks;
my $numTicks;
my $chr_length;
my $chr_length_reel;
my $scale_factor;
my $strand_width;
my $space_between_str;
my $label_strand_rotation;
my $picHeight;
my $picWidth;
my $colour_scale;
my $count = 0;
my $win_size;
my $rounding_method;

my %color;
my %gffTypes;
my %margin;
my %order;
$order{'-'} = "-";
$order{'+'} = "+";
$order{'both'} = "-;+";
$order{'fused'} = "-+";
$order{'all'} = "-;-+;+";
my %rand;
my %offset;
$offset{'x'} = 0;
$offset{'y'} = 0;
$offset{'x1'} = 0;
$offset{'y1'} = 0;
$offset{'x2'} = 0;
$offset{'y2'} = 0;


$| = 1;


if (!exists $config{strand_width})          {$strand_width = 50;}
else                                        {$strand_width = $config{strand_width};}

if (!exists $config{show_scale})            {$numMaxTicks = 50;}
else                                        {$numMaxTicks = $config{show_scale};}

if (!exists $config{space_between_str})     {$space_between_str = 50;}
else                                        {$space_between_str = $config{space_between_str};}

if    ($config{auto_scale_factor})          {$scale_factor = 1;}
elsif (!exists $config{scale_factor})       {$scale_factor = 1000;}
else                                        {$scale_factor = $config{scale_factor};}

if (!exists $config{lmargin})               {$margin{l} = 50;}
else                                        {$margin{l} = $config{lmargin};}

if (!exists $config{rmargin})               {$margin{r} = 50;}
else                                        {$margin{r} = $config{rmargin};}

if (!exists $config{tmargin})               {$margin{t} = 50;}
else                                        {$margin{t} = $config{tmargin};}

if (!exists $config{bmargin})               {$margin{b} = 50;}
else                                        {$margin{b} = $config{bmargin};}

if (!exists $config{label_strand_rotation}) {$label_strand_rotation = 0;}
else                                        {$label_strand_rotation = $config{label_strand_rotation};}

if (!exists $config{colour_scale})         {$colour_scale = 7;}
else                                        {$colour_scale = $config{colour_scale};}

if (!exists $config{win_size})              {$win_size = 1;}
else                                        {$win_size = $config{win_size};}

if (!exists $config{rounding_method})       {$rounding_method = "floor";}
else                                        {$rounding_method = $config{rounding_method};}


print "gffs : $config{gff}\n" if $config{verbose};
print "output_img_name : $config{output_img_name}\n" if $config{verbose};
##> Setting parameters


# 1. Get Image size
## 1.1 Height (+Check GFF validity)
print "Searching Max Sequence Length ... \n" if $config{debug};

my $numOfGff = 0;
my $maxSequenceLength = 0;

foreach my $file (split(/;/, $config{'gff'})){
    print "\tLooking at $file \n" if $config{debug};
    
    open(GFF, "<$file") or die "Can not open $file ! \n";
    $numOfGff++;
    
    my $first_line = <GFF>;
    printError ("Not GFF3 format ! (1)\n", 1) if $first_line !~ /##gff-version 3/;
    
    $first_line = <GFF>;
    printError ("Not GFF3 format ! (2)\n", 1) if $first_line !~ /##sequence-region\s+[^\s]+\s+\d+\s+(\d+)/;
    
    $maxSequenceLength = ($1 > $maxSequenceLength) ? $1 : $maxSequenceLength;
    close GFF;
}

$margin{'t'} += 40 if ($config{title});
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
foreach (split(/;/, $config{'type_to_draw'})){
    if($_=~ /([^=]+)=([^=]+)/){
        push @type_array, $1;
        
        if   ($2 eq "-" or $2 eq "+" or $2 eq "fused")  {$numOfStrand++;}
        elsif($2 eq "both")                             {$numOfStrand += 2;}
        elsif($2 eq "all")                              {$numOfStrand += 3;}
        else                                            {printError("Type option : INVALID STRAND FORMAT, check help please.\n", 1);}
    }
    else {printError("Type option : INVALID FORMAT, check help please.\n", 1);}
}

$margin{'l'} = $margin{'l'} + $scaleAddWidth if ($config{show_scale});
$picWidth = $margin{'l'}
          + $margin{'r'}
          #+ (($numOfGff * $numOfStrand) * ($strand_width + $space_between_str));
          + (($numOfGff * $numOfStrand    ) * $strand_width)
          + (($numOfGff * $numOfStrand - 1) * $space_between_str);

## 1.3 Ask if image size is OK
if (! $config{yes}) {
    while (1) {
        print "Is the picture size will be ok ? ($picWidth px x $picHeight px) [y/n] : ";
        chomp (my $response = <STDIN>);
        if    ($response eq "n") {exit;}
        elsif ($response eq "y") {last;}
    }
}

# 2 Create picture
print "Create Picture ...\n" if $config{'verbose'};
my $image = GD::SVG::Image->new($picWidth, $picHeight);


# 3 Loading colors from colours.txt
print "Load colours ...\n" if $config{'verbose'};
open(COLOR, "<./colours.txt") or die "Can not open colours.txt";
while (<COLOR>) {
    next if /^#/;
    /([\d\w]+);(\d+);(\d+);(\d+)/;
    $color{lc($1)} = $image->colorAllocate($2,$3,$4);
}
close COLOR;

# 4 Draw Background
if ($config{'background'}) {
    print "Draw background ...\n" if $config{'verbose'};
    $image->filledRectangle(0, 0, $picWidth, $picHeight, $color{$config{'background'}});
}


# 5 Drawing Title and Scale
if ($config{title}){
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
    print "Drawing Scale ...\n" if $config{verbose};
    my $scale_size = floor($maxSequenceLength/$scale_factor);
    print "Scale_size : $scale_size\n" if $config{debug};
    drawScale(\$image, $scale_size, $maxSequenceLength, $scaleAddWidth, $scale_factor, $numMaxTicks, $margin{'t'}, $win_size);
}


# 6 Get Type/Strand
foreach (split(/;/, $config{'type_to_draw'})){
    $_=~ /([^=]+)=([^=]+)/;

    print "type = $1\tstrand = $2\n" if $config{debug};

    $gffTypes{$1}{strand} = $2;
    $gffTypes{$1}{'-'}  = [];
    $gffTypes{$1}{'+'}  = [];
    $gffTypes{$1}{'-+'} = [];

    $paternType .= '\t'.$1.'\t|';
}
#$paternType =~ s/\|$//g;
$paternType .= 'centromere';
print "Patern Regexp = $paternType\n" if $config{debug};


# 7 Foreach file draw strand(s)
foreach my $file (split(/;/, $config{'gff'})){
    
    # 7.1 Load GFF
    print "Check GFF and load chromosome size ...\n" if $config{'verbose'};
    open(GFF, "<$file") or die "Cannot open $file ! \n";
    
    my $first_line = <GFF>;
       $first_line = <GFF>;
    
    $first_line =~ /##sequence-region\s+([^\s]+)\s+\d+\s+(\d+)/;
    $chr_length_reel = $2;
    $chr_length = floor($chr_length_reel/$scale_factor);
    
    my $seqName = $1;
    
    print "Loading $seqName ...\n"  if $config{verbose};
    
    my $boolOK = 1;
    
    # clean intervals sets
    foreach (keys(%gffTypes)){
        $gffTypes{$_}{'-'}  = [];
        $gffTypes{$_}{'+'}  = [];
        $gffTypes{$_}{'-+'} = [];
    }
    
    my %centromere;
    
    while (<GFF>) {
        next if $_ !~ /$paternType/;
        chomp;
        my @line = split(/\t/);
        
        if ($line[2] eq "centromere") {
            $centromere{start} = $line[3];
            $centromere{end}   = $line[4];
            #print "START ========> $centromere{start}\n";
            #print "END\  ========> $centromere{end}\n";
            next;
        }
        
        my $ref_tabM =  $gffTypes{$line[2]}{'-'};
        my $ref_tabP =  $gffTypes{$line[2]}{'+'};
        my $ref_tabMP = $gffTypes{$line[2]}{'-+'};
        
        if    ($line[6] eq "-") {push(@{$ref_tabM},  [$line[3], $line[4]]);}
        elsif ($line[6] eq "+") {push(@{$ref_tabP},  [$line[3], $line[4]]);}
                                 push(@{$ref_tabMP}, [$line[3], $line[4]]);
        
        #print $_."\n" if (!defined($line[3]) or !defined($line[4]) or !defined($line[6]));
        #$gffTab[line][0] = $line[3]; # gff start
        #$gffTab[line][1] = $line[4]; # gff end
        #$gffTab[line][2] = $line[6]; # gff strand
    }
    close GFF;
    
    # 7.2 Sorting intervals
    foreach (keys(%gffTypes)){
        my $ref_tabM =  $gffTypes{$_}{'-'};
        my $ref_tabP =  $gffTypes{$_}{'+'};
        my $ref_tabMP = $gffTypes{$_}{'-+'};
        
        @{$ref_tabM}  = sort {$a->[0] <=> $b->[0]} @{$ref_tabM};
        @{$ref_tabP}  = sort {$a->[0] <=> $b->[0]} @{$ref_tabP};
        @{$ref_tabMP} = sort {$a->[0] <=> $b->[0]} @{$ref_tabMP};
    }
    
    # 7.3 Reducing intervals
    foreach (keys(%gffTypes)){
        $gffTypes{$_}{'-'}  = removeIntervalRedundancy(@{$gffTypes{$_}{'-'}});
        $gffTypes{$_}{'+'}  = removeIntervalRedundancy(@{$gffTypes{$_}{'+'}});
        $gffTypes{$_}{'-+'} = removeIntervalRedundancy(@{$gffTypes{$_}{'-+'}});
    }
    
    # 7.4 Set Offset
    #foreach (keys(%gffTypes)){
    foreach (@type_array){
        if ($gffTypes{$_}{'strand'} eq "-")     {
            $offset{$_}{'-'}{'x'} = $margin{'l'} + ($count * ($strand_width + $space_between_str));
            $offset{$_}{'-'}{'y'} = $margin{'t'};
            $count++;
            
            print "x = $offset{$_}{'-'}{'x'}\n" if $config{debug};
            print "y = $offset{$_}{'-'}{'y'}\n" if $config{debug};
        }
        if ($gffTypes{$_}{'strand'} eq "+")     {
            $offset{$_}{'+'}{'x'} = $margin{'l'} + ($count * ($strand_width + $space_between_str));
            $offset{$_}{'+'}{'y'} = $margin{'t'};
            $count++;
            
            print "+ = $offset{$_}{'+'}{'x'}\n" if $config{debug};
            print "+ = $offset{$_}{'+'}{'y'}\n" if $config{debug};
        }
        if ($gffTypes{$_}{'strand'} eq "both")  {
            $offset{$_}{'-'}{'x'} = $margin{'l'} + ($count * ($strand_width + $space_between_str));
            $offset{$_}{'-'}{'y'} = $margin{'t'};
            $count++;
            
            $offset{$_}{'+'}{'x'} = $margin{'l'} + ($count * ($strand_width + $space_between_str));
            $offset{$_}{'+'}{'y'} = $margin{'t'};
            $count++;
            
            print "+ = $offset{$_}{'-'}{'x'}\n" if $config{debug};
            print "+ = $offset{$_}{'-'}{'y'}\n" if $config{debug};
            print "- = $offset{$_}{'+'}{'x'}\n" if $config{debug};
            print "- = $offset{$_}{'+'}{'y'}\n" if $config{debug};
        }
        if ($gffTypes{$_}{'strand'} eq "fused") {
            $offset{$_}{'-+'}{'x'} = $margin{'l'} + ($count * ($strand_width + $space_between_str));
            $offset{$_}{'-+'}{'y'} = $margin{'t'};
            $count++;
            
            print "f = $offset{$_}{'f'}{'x'}\n" if $config{debug};
            print "f = $offset{$_}{'f'}{'y'}\n" if $config{debug};
        }
        if ($gffTypes{$_}{'strand'} eq "all")   {
            $offset{$_}{'-'}{'x'} = $margin{'l'} + ($count * ($strand_width + $space_between_str));
            $offset{$_}{'-'}{'y'} = $margin{'t'};
            $count++;
            
            $offset{$_}{'-+'}{'x'} = $margin{'l'} + ($count * ($strand_width + $space_between_str));
            $offset{$_}{'-+'}{'y'} = $margin{'t'};
            $count++;
            
            $offset{$_}{'+'}{'x'} = $margin{'l'} + ($count * ($strand_width + $space_between_str));
            $offset{$_}{'+'}{'y'} = $margin{'t'};
            $count++;
            
            print "- = $offset{$_}{'-'}{'x'}\n"  if $config{debug};
            print "- = $offset{$_}{'-'}{'y'}\n"  if $config{debug};
            print "f = $offset{$_}{'-+'}{'x'}\n" if $config{debug};
            print "f = $offset{$_}{'-+'}{'y'}\n" if $config{debug};
            print "+ = $offset{$_}{'+'}{'x'}\n"  if $config{debug};
            print "+ = $offset{$_}{'+'}{'y'}\n"  if $config{debug};
        }  
    }
    
    # 7.5 Foearch type draw strands
    #foreach my $typeToDraw (keys(%gffTypes)){
        
        #WAS USED FOR MERGING + AND - TRACK AND THEN ADDING TRANSPARENCY ON THE MERGED TRACKS
        #print "typeToDraw = $typeToDraw\n"  if $config{'debug'};
        #
        #foreach my $strand (split(";", $order{$gffTypes{$typeToDraw}{'strand'}})){
        #    print "strand = $strand\n"      if $config{'debug'};
        #    
        #    if ($strand ne "-+") {
        #        my $ref_tab = $gffTypes{$typeToDraw}{$strand};
        #        drawPixels(\$image, \%rand, $chr_length, $scale_factor, $typeToDraw, $strand, $strand, $ref_tab);
        #    }
        #    elsif ($strand eq "-+"){
        #        my $ref_tab = $gffTypes{$typeToDraw}{"-"};
        #        drawPixels(\$image, \%rand, $chr_length, $scale_factor, $typeToDraw, $strand, "-", $ref_tab);
        #        
        #        $ref_tab = $gffTypes{$typeToDraw}{"+"};
        #        drawPixels(\$image, \%rand, $chr_length, $scale_factor, $typeToDraw, $strand, "+", $ref_tab);
        #    }
        #}
        
        foreach my $typeToDraw (@type_array){
            print "typeToDraw = $typeToDraw\n"  if $config{'debug'};
            
            foreach my $strand (split(";", $order{$gffTypes{$typeToDraw}{'strand'}})){
                print "strand = $strand\n"      if $config{'debug'};
                
                my $ref_tab = $gffTypes{$typeToDraw}{$strand};
                drawPixels(\$image, \%rand, \$colour_scale, $seqName, $chr_length, $scale_factor, $typeToDraw, $strand, $strand, \%centromere, $ref_tab, $win_size);
                
            }
        }
    #}   
} # END 7


# 8 Applying rotation to labels and Saving picture
open(IMG, ">$config{output_img_name}.svg") or die "Can not open $config{output_img_name} ! ";
binmode IMG;
my @text = split("\t", $image->svg);
#if (!$config{transparency}) {
    my $switchRotate = 0;
    foreach (@text){
        if    ($config{label_strand_rotation} and /<g id="rotate_\d+">/)    {$switchRotate++;}
        elsif ($config{label_strand_rotation} and/<\/g>/ and $switchRotate) {$switchRotate--;}
        
        if ($switchRotate) {s/x="(\d+)" y="(\d+)"/x="$1" y="$2" transform="rotate($label_strand_rotation, $1, $2)"/;}
        
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
sub removeIntervalRedundancy{
    
    # The purpose of this function is to remove or merging the feature of the gff
    # Input  : sorted array of start - end array
    # Output : sorted array of start - end array without crossing intervals
    
    #my ($ref_intervalsFused, @gff) = @_;
    my (@gff) = @_;
    my @reduced_gff;
    my $bool_firstInterval = 1;
    
    print "Start removeIntervalRedundancy ... \n"   if $config{'verbose'};
    
    print "\tStarting gff size = ".@gff."\n"        if $config{'verbose'};
    
    
    foreach my $annot (@gff) {
        
        # The intervals of the gff have been sorted by 'start' before using this function
        # so only three case are possibles :
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
        
        print "\t---> Current interval = $annot->[0], $annot->[1]\n"            if $config{insaneDebugMode};
        
        #Load first interval
        if ($bool_firstInterval) {
            print "\nLoad first interval\n"                                     if $config{'insaneDebugMode'};
            
            push(@reduced_gff, [$annot->[0], $annot->[1]]);
            $bool_firstInterval--;
            
            print "\tannot = $annot\n"                                          if $config{insaneDebugMode};
            print "\t\$annot->[0], \$annot->[1] = $annot->[0], $annot->[1]\n"   if $config{insaneDebugMode}; 
            print "\t\@reduced_gff size = ".@reduced_gff."\n\n"                 if $config{insaneDebugMode};
            
            next;
        }
        
        if ($annot->[0] <= $reduced_gff[$#reduced_gff]->[1] and     #1
            $annot->[1] >  $reduced_gff[$#reduced_gff]->[1]) {
            # Replace previous end
            $reduced_gff[$#reduced_gff]->[1] = $annot->[1];
        }
        elsif ($annot->[0] <= $reduced_gff[$#reduced_gff]->[1] and  #2
               $annot->[1] <=  $reduced_gff[$#reduced_gff]->[1]) {
            # Go to next interval
            next;
        }
        else {                                                      #3
            # Insert interval in reduced gff
            push(@reduced_gff, [$annot->[0], $annot->[1]]);
        }
    } #End foreach (@gff)
    
    print "\tFinal gff size    = ".@reduced_gff."\n" if $config{'verbose'};
    print "Finished\n" if $config{'verbose'};

    return \@reduced_gff;
}

###########################################################################
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
    print "Start searching num ticks\n" if $config{debug};
    while (1) {
        $numTicks = floor($chr_size/$basesPerTicks);
        print "basesPerTicks = $basesPerTicks \t numTicks = $numTicks\n"  if $config{debug};
        last if ($numTicks <= $maxTicks);
        $basesPerTicks*=10; 
    }
    print "Found num ticks\n" if $config{debug};
    
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
    print "\$\$ref_img->filledRectangle($widthScale - 10, $marginTop, $widthScale - 8, $marginTop + $chr_size, \$color{'black'})\;\n" if $config{debug};
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
    print "\$\$ref_img->filledRectangle($widthScale - 10, $marginTop + $chr_size - 1, $widthScale, $marginTop + $chr_size, \$color{'black'})\;\n" if $config{debug};
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
         
        print "\$\$ref_img->filledRectangle($widthScale - 10, $marginTop + ($basesPerTicks * $currentTick), $widthScale, $marginTop + ($basesPerTicks * $currentTick) + 1, \$color{'black'})\;\n" if $config{debug};
        $$ref_img->filledRectangle($widthScale - 10,
                                $marginTop + ($basesPerTicks * $currentTick) * $win_size,
                                $widthScale,
                                $marginTop + ($basesPerTicks * $currentTick) * $win_size + 1,
                                $color{'black'});
    }   
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
    #   - $ref_gff      ->  ref of the cureent strand gff
    # Ouput : none
    
    print "Start Drawing pixels ...\n" if $config{'verbose'};
    
    my ($ref_img, $ref_rand, $ref_colour_scale, $seqName, $chr_size, $scaleFactor, $type, $strand, $strandColor, $ref_centromere, $ref_gff, $win_size) = @_;
    my @gff    = @{$ref_gff};
    my %centro = %{$ref_centromere};
    my %strand2color;
    my $randNum;
    my %intervals;
    my @previousBases = (0);
    
    # Search a unique random number
    while (1) {
        $randNum = int rand(1000);
        redo if $$ref_rand{$randNum}++;
        last;
    }
    
    # Set colors for each strands
    $strand2color{'+'}  = $$ref_colour_scale;
    $strand2color{'-'}  = $$ref_colour_scale;
    $strand2color{'-+'} = $$ref_colour_scale;
    
    print "\tchr_size    = $chr_size\n"       if $config{'debug'};
    print "\tscaleFactor = $scaleFactor\n"    if $config{'debug'};
    print "\ttype        = $type\n"           if $config{'debug'};
    print "\tstrand      = $strand\n"         if $config{'debug'};
    
    # Open chromosome/sequence group
    $$ref_img->startGroup("${strand}_${randNum}");
    
    # For each pixel of the chromosome/sequence
    for (my $pos = 0 ; $pos <= $chr_size ; $pos++) {
        
        # Get number of base covered by the previous on pixel crosssing interval
        my $basesCoverred = shift @previousBases;
        $basesCoverred = 0 if (!defined $basesCoverred); # if @previousBases is empty
        
        print "pos = $pos\n"                                            if $config{'debug'};
        print "basesCoverred Tab = -".join(" ", @previousBases)."-\n"   if $config{'debug'};        
        print "basesCoverred = $basesCoverred (previous shift)\n"       if $config{'debug'};
        print "Start while\n"                                           if $config{'debug'};
        print "\$gff[0]->[0] = $gff[0]->[0]\n"                          if $config{'debug'};
        #printError( "\$gff[0]->[0] undef ! \n") if(!defined($gff[0]->[0]));
        
        # while the end of gff is not reached and the next start is in the current pixel
        while (defined($gff[0]->[0]) and $gff[0]->[0] < ($pos * $scaleFactor)) {
            my $ref_interval= shift(@gff);
            
            last if (!defined $ref_interval->[0]); # = defined($gff[0]->[0]) avoid warnings
            
            # get clearly start and end 
            $intervals{'start_reel'}  = $ref_interval->[0];
            $intervals{'end_reel'}    = $ref_interval->[1];
            
            print "\tstart_reel = $intervals{'start_reel'}\n"           if $config{'debug'};
            print "\tend_reel   = $intervals{'end_reel'}\n"             if $config{'debug'};
            print "\tcheck position\n"                                  if $config{'debug'};
            
            # feature end is on current pixel
            if ($intervals{'end_reel'} < ($pos * $scaleFactor) ){ 
                print "\tCase 1 : \$basesCoverred += $intervals{'end_reel'} - $intervals{'start_reel'}\n" if $config{'debug'};
                print "\tCase 1 : \$basesCoverred += ".($intervals{'end_reel'} - $intervals{'start_reel'})."\n" if $config{'debug'};
                $basesCoverred += $intervals{'end_reel'} - $intervals{'start_reel'}; # increment the base counting
            }
            
            # feature end is on next pixels
            else{
                # determine on how much pixel the interval is crossing
                my $numPixelImpliy = floor(($intervals{'end_reel'} - $intervals{'start_reel'} - (($pos * $scaleFactor) - $intervals{'start_reel'}))/$scaleFactor);
                
                # Treating current position
                $basesCoverred += ($pos * $scaleFactor) - $intervals{'start_reel'}; 
                print "\tCase 2 : \$basesCoverred += ($pos * $scaleFactor) - $intervals{'start_reel'}\n" if $config{'debug'};
                
                # Treating the next full pixels 
                for (my $i = 0 ; $i < $numPixelImpliy ; $i++){
                    print "\tCase 2 : \$previousBases[$i] += $scaleFactor\n" if $config{'debug'};
                    $previousBases[$i] += $scaleFactor;
                }
                
                # Treat the last not full pixel
                $previousBases[$numPixelImpliy] += (($intervals{'end_reel'} - $intervals{'start_reel'} - (($pos * $scaleFactor) - $intervals{'start_reel'}))%$scaleFactor);
                print "\tCase 2 : \$previousBases[$numPixelImpliy] += (($intervals{'end_reel'} - $intervals{'start_reel'} - (($pos * $scaleFactor) - $intervals{'start_reel'}))%$scaleFactor);\n" if $config{'debug'};
                
            }# feature end is on next pixels
            
            print "\tbasesCoverred = $basesCoverred\n"      if $config{'debug'};
            
        }# while ($intervalsTabM[0]->[0] < ($pos * $scaleFactor))
        
        print "Stop while\n"                                if $config{'debug'};
        
        # Compute percentage of base coverage
        my $percentage;
        if    ($rounding_method eq "floor") {$percentage = floor(($basesCoverred/$scaleFactor) * 100);}
        elsif ($rounding_method eq "ceil" ) {$percentage = ceil (($basesCoverred/$scaleFactor) * 100);}
        
        print "percentage = $basesCoverred/$scaleFactor = ".($basesCoverred/$scaleFactor)."\n"  if $config{'debug'};
        print "percentage = $percentage % \n\n"             if $config{'debug'};
        print "$pos = $strand2color{$strand}_heatmap$percentage\n" if $config{insaneDebugMode};
        
        # kill if more than 100 % 
        if ($percentage > 100) {printError("Higher than 100 % ($percentage % )", 1);}
        
        # Draw the current pixel
        $$ref_img->filledRectangle($offset{$type}{$strand}{'x'},                    $offset{$type}{$strand}{'y'} + $pos * $win_size,
                                   $offset{$type}{$strand}{'x'} + $strand_width,    $offset{$type}{$strand}{'y'} + $pos * $win_size + $win_size,
                                   $color{"$strand2color{$strandColor}_heatmap$percentage"});
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
                   $offset{$type}{$strand}{'x'}, 
                   $offset{$type}{$strand}{'y'} - (3.5 * gdLargeFont->height), 
                   $seqName,
                   $color{'black'});
    
    # Draw label type
    $$ref_img->string(gdLargeFont,
                   $offset{$type}{$strand}{'x'}, 
                   $offset{$type}{$strand}{'y'} - (2 * gdLargeFont->height), 
                   $type,
                   $color{'black'});
    # Close label group for next rotation
    $$ref_img->endGroup;
    
    # Draw strand label
    $$ref_img->string(gdLargeFont,
                   $offset{$type}{$strand}{'x'}, 
                   $offset{$type}{$strand}{'y'} - gdLargeFont->height, 
                   $strand,
                   $color{'black'});
    
    # close chromosome/sequence group
    $$ref_img->endGroup;
    print "Finish Drawing pixels.\n" if $config{'verbose'};
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
"USAGE: DensityMap.pl -gff chromosome.gff3 -ty 'match=all' -o Chromosome
    
Options:
    -g   | gff string                  Gff file (version gff3 only !)          			(Mandatory)
                                       Format: file.gff or \"file1.gff;file2.gff\"
    -o   | output_img_name string      Name of the output image                			(Mandatory)
    -ty  | type_to_draw string         List of type (column 3 of gff) to draw and strand to use (Mandatory)
                                       Type: match, gene, CDS, ...
                                       Strand:  - -        -> strand -
                                                - +        -> strand +
                                                - both     -> strand - and strand +
                                                - fused    -> Combination of strand - and strand +
                                                - all      -> strand - and strand + and fused
                                       Format: \"match=all;gene=both;CDS=fused\"
    -v   | verbose                     MORE text dude !!!!
    -hel | help                        This help

Density options: 
    -c   | coulour_scheme int          color scheme to use 	       (Default = 7)
    -sc  | scale_factor int            = window length in bp           (Default = 1000)
    -a   | auto_scale_factor int       Max picture height in pixel
    -ro  | rounding_method string      floor or ceil		       (Default = floor)

Graphical options: 
    -ti  | title string                Title to print on the picture
    -w   | win_size int                Height of window in pixel       (Default = 1)
    -sh  | show_scale int              Draw Scale, n = num max ticks   (Default = 50)
    -st  | strand_width int            Strand width in pixel           (Default = 50)
    -lm  | lmargin int                 Left margin in pixel            (Default = 50)
    -rm  | rmargin int                 Rigth margin in pixel           (Default = 50)
    -tm  | tmargin int                 Top margin in pixel             (Default = 50)
    -bm  | bmargin int                 Bottom margin in pixel          (Default = 50)
    -sp  | space_between_str int       Space between strands in pixel  (Default = 50)
    -ba  | background color            Fill Background                 (Default = no Background)
    -la  | label_strand_rotation int   Rotation degree of strand label (Default = 0)



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
