use strict;
use warnings;
use Text::ParseWords;
#-----------------------------------------------------------------------------80
# Global variables
#-----------------------------------------------------------------------------80
my $in  = "../rawdata/2020_001_Reyes/SCP548/expression/tmp_01.csv";
my $tmp = "../rawdata/2020_001_Reyes/SCP548/expression/tmp_02.csv";
my $out = "../rawdata/2020_001_Reyes/SCP548/expression/tmp_03.csv";
my @nreads = ();
#=============================================================================80
# Main 1
#=============================================================================80
open( IN,       $in  ) or die "$!";
open( TMP, ">", $tmp ) or die "$!";
#-----------------------------------------------50
# Header
#-----------------------------------------------50
my $header = <IN>;
chomp( $header );
print TMP $header, "\n";
my @parsed_line = &parse_line( ',', undef, $header );
my $ncol = @parsed_line;
push( @nreads, "0" );
push( @nreads, "nReads" );
foreach my $i (2..($ncol-1)){
  push( @nreads, 0 );
}
#-----------------------------------------------50
# Below header
#-----------------------------------------------50
while( <IN> ){
  chomp;
  @parsed_line = &parse_line( ',', undef, $_ );
  $nreads[0] += $parsed_line[0];
  foreach my $i (2..($ncol-1)){
    $nreads[$i] += $parsed_line[$i];
  }
  print TMP $_, "\n";
}
close( TMP );
close( IN  );
#=============================================================================80
# Main 2
#=============================================================================80
open( IN,       $tmp ) or die "$!";
open( OUT, ">", $out ) or die "$!";
#-----------------------------------------------50
# Header
#-----------------------------------------------50
$header = <IN>;
chomp( $header );
print OUT $header, "\n";
print OUT join( "\,", @nreads), "\n";
#-----------------------------------------------50
# Below header
#-----------------------------------------------50
while( <IN> ){
  chomp;
  print OUT $_, "\n";
}
close( TMP );
close( IN  );
#-----------------------------------------------50
# Output log.
#-----------------------------------------------50
print "The following files were generated: ", "\n";
print "* " . $tmp, "\n";
print "* " . $out, "\n\n";
