use strict;
use warnings;
use Text::ParseWords;
#-----------------------------------------------------------------------------80
# Global variables
#-----------------------------------------------------------------------------80
my $in = "../rawdata/2020_001_Reyes/SCP548/expression/scp_gex_matrix_raw.csv";
my $out = "../rawdata/2020_001_Reyes/SCP548/expression/tmp_01.csv";
#-----------------------------------------------------------------------------80
# Compute the number of non-zero expressing samples.
#-----------------------------------------------------------------------------80
sub comp_nsamples{
  my ($parsed_line) = @_;
  my $ncol = @$parsed_line;
  my $nsamples = 0;
  foreach my $i (1..($ncol-1)){
    $nsamples += $$parsed_line[$i];
  }
  return $nsamples;
}
#=============================================================================80
# Main
#=============================================================================80
open( IN,       $in  ) or die "$!";
open( OUT, ">", $out ) or die "$!";
#-----------------------------------------------50
# Header
#-----------------------------------------------50
my $header = <IN>;
chomp( $header );
my @parsed_line = &parse_line( ',', undef, $header );
my @parsed_line_out = ();
push( @parsed_line_out, "nSamples" );
foreach my $element (@parsed_line){
  push( @parsed_line_out, $element );
}
print OUT join( "\,", @parsed_line_out ), "\n";
#-----------------------------------------------50
# Below header
#-----------------------------------------------50
while( <IN> ){
  chomp;
  @parsed_line = &parse_line( ',', undef, $_ );
  @parsed_line_out = ();

  push( @parsed_line_out, comp_nsamples( \@parsed_line ) );
  foreach my $element (@parsed_line){
    push( @parsed_line_out, $element );
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
