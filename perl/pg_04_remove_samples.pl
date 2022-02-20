use strict;
use warnings;
use Text::ParseWords;
#-----------------------------------------------------------------------------80
# Global variables
#-----------------------------------------------------------------------------80
my $in  = "../rawdata/2020_001_Reyes/SCP548/expression/tmp_04.csv";
my $out = "../rawdata/2020_001_Reyes/SCP548/expression/scp_gex_matrix_red.csv";
my $MIN_NREADS = 1000;
#=============================================================================80
# Main 1
#=============================================================================80
open( IN, $in ) or die "$!";
#-----------------------------------------------50
# Header
#-----------------------------------------------50
my $header = <IN>;
$header = <IN>;
chomp( $header );
my @nreads = &parse_line( ',', undef, $header );
close( IN );
#=============================================================================80
# Main 2
#=============================================================================80
open( IN,       $in  ) or die "$!";
open( OUT, ">", $out ) or die "$!";

while( <IN> ){
  chomp;
  my @parsed_line = &parse_line( ',', undef, $_ );
  my $ncol = @parsed_line;
  my @parsed_line_out = ();
  foreach my $i (0..1){
    push( @parsed_line_out, $parsed_line[$i] );
  }
  foreach my $i (2..($ncol-1)){
    if( $nreads[$i] >= $MIN_NREADS ){
      push( @parsed_line_out, $parsed_line[$i] );
    }else{
      next;
    }
  }
  print OUT join( "\,", @parsed_line_out ), "\n";
}
close( OUT );
close( IN  );
#-----------------------------------------------50
# Output log.
#-----------------------------------------------50
print "The following file was generated: ", "\n";
print "* " . $out, "\n\n";
