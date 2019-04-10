#!/usr/bin/perl -w 
#use strict;
###### Get folder , parameter, and temp. information   ########
open(INPUT,"<$ARGV[0]") || die "cannot find $ARGV[0]\n";

my $mol_ener ="-538.6603856" ;      ## clet MP3=-538.6603856
my $he_ener ="-2.8910899" ;         ## He   MP3=-2.8910899
my $ne_ener ="-128.7072031" ;       ## Ne   MP3=-128.7072031
my $rare_gas = 0; 
if ($ARGV[1] eq "he"){
  $rare_gas =  $he_ener;
}
if ( $ARGV[1] eq "ne")
{
 $rare_gas =  $ne_ener;
} 
 
my $line = "";
print "#### $ARGV[1]  $rare_gas\n";
while ( $line !~ /^ Summary of the potential surface scan/){
	$line = <INPUT> ; 
}
	$line = <INPUT> ; 
	$line = <INPUT> ; 
	$line = <INPUT> ;    ## first one; 
my %mp3_ener ;   
while ( $line !~ /^ ----/ ){
     my @temp_array = split(/\s+/, $line);

     my $inter = ($temp_array[5] - $mol_ener - $rare_gas)*627.50959; ## convert to kcal/mol
     $mp3_energ{sprintf("%.1f", $temp_array[2])} = $inter;
##     $mp3_energ{$temp_array[2]} = $temp_array[5];
     $line = <INPUT> ; 
}

my $i=1.5;
for ( $i=1.5; $i<= 6.0; $i=sprintf("%.1f", $i+0.1)){
print "$i  $mp3_energ{$i} \n";
}

