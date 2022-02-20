use strict;
use warnings;
use Text::ParseWords;
#-----------------------------------------------------------------------------80
# Global variables
#-----------------------------------------------------------------------------80
my $in  = "../rawdata/2020_001_Reyes/SCP548/expression/tmp_03.csv";
my $out = "../rawdata/2020_001_Reyes/SCP548/expression/tmp_04.csv";
my $MIN_NSAMPLES = 100;
#=============================================================================80
# Main 1
#=============================================================================80
open( IN,       $in  ) or die "$!";
open( OUT, ">", $out ) or die "$!";
#-----------------------------------------------50
# Header
#-----------------------------------------------50
my $header = <IN>;
chomp( $header );
print OUT $header, "\n";
#-----------------------------------------------50
# Below header
#-----------------------------------------------50
while( <IN> ){
  chomp;
  my @parsed_line = &parse_line( ',', undef, $_ );
  if( $parsed_line[0] >= $MIN_NSAMPLES ){
    print OUT $_, "\n";
  }else{
    next;
  }
}
close( OUT );
close( IN  );
#-----------------------------------------------50
# Output log.
#-----------------------------------------------50
print "The following file was generated: ", "\n";
print "* " . $out, "\n\n";
