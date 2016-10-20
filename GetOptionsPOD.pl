
GetOptions (
         'rpath|r=s' => \$rpath,
         'genome|g=s' => \$g,
         'outdir|o=s' => \$out,
         'threads|t:i' => \$t,
         'insert-size|is:i' => \$meanIn,
         'reverse_seq|rs:s' => \$reverse_seq,
         'version|v' => \$opt_versions,
         'help|h' => \$help,
         'man' => \$man,
) or pod2usage({-verbose => 1, -message => "Invalid arguments.\n",}) && exit;

 #   Check for requests for help or for man (full documentation):
pod2usage( -verbose => 1 ) && exit if defined $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) && exit if defined $man;

#   Check for required variables.
if ( !($rpath) or !($g) ) { 
   pod2usage({-exitstatus => 2, -verbose => 1, -message=>"\n$0: Missing arguments.\n",});
   exit;
}

