#!/usr/bin/perl -w

use strict;
use diagnostics;
use warnings;
use Getopt::Long;
use Term::ANSIColor;
use Data::Dumper;
use GD::SVG;
use POSIX;

#> Setting Parameters


##> Define outputs colors
print STDOUT color 'blue';
print STDERR color 'red';


##> Define options ---> NO OPTION


##> Print USAGE if --help


##> Check if gff file exist, if no mandatory parameter are missing


##> Setting Global Variables
my %color;

my $margin = 5;
my $x = $margin;
my $y = $margin;
my $height_scale = gdSmallFont->height;
my $space_between_scales = 5;
my $width_pixel = 2;
my $height;
my $width;

$| = 1;
##> Setting parameters

open(COLOR, "<colours.txt") or die printError("Cannot open colours.txt ! \n", 1);
while (<COLOR>) {
    next if !(/heatmap/);
    /(\d+)_heatmap(\d+);(\d+);(\d+);(\d+)/;
    $color{$1}{$2}{r} = $3;
    $color{$1}{$2}{g} = $4;
    $color{$1}{$2}{b} = $5;
    #print "$1 - $2 - $3 - $4 - $5\n";
}
close COLOR;

$height = $margin * 2 + scalar(keys(%color)) * $height_scale + (scalar(keys(%color)) - 1 ) * $space_between_scales;
$width  = $margin * 2 + $width_pixel * 101 + gdSmallFont->width * 3; 
my $image = GD::SVG::Image->new($width, $height);


my $white = $image->colorAllocate(255,255,255);
my $grey  = $image->colorAllocate(84,84,84);
my $black = $image->colorAllocate(0,0,0);


foreach my $sc (sort {$a <=> $b} keys(%color)){
    $image->string( gdSmallFont, $x, $y, $sc, $black );
    $x += gdSmallFont->width * 3;
    
    for (my $i = 0 ; $i < 101 ; $i++){
        my $c = $image->colorAllocate($color{$sc}{$i}{r}, $color{$sc}{$i}{g}, $color{$sc}{$i}{b});
        $image->filledRectangle($x, $y, $x + 1, $y + $height_scale, $c);
        $x += $width_pixel;
    }
    $x = $margin;
    $y += $height_scale + $space_between_scales;
}

open(IMG, ">scales.svg");
binmode IMG;
print IMG $image->svg;
close IMG;

###########################################################################
################################ Fonctions ################################
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
"USAGE : scaleColorDrawer.pl
Draw avaiable scales.
\n\n";
    exit if $exit;
}
