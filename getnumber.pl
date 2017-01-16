#!/usr/bin/perl


use warnings;
#use strict;
use Fcntl;

sub main{
    $mypath="/home/dailey.110/analysis_oindree/";
    $j=0;
$HANDLE_in="filehandle_in";
$FilePath_in="$mypath/output.txt";
sysopen(HANDLE_in,$FilePath_in,O_RDONLY) or die "death"; # read only

while (<HANDLE_in>) {
    $tmp=$_;
    if ($tmp =~ /\bdailey\b/ig) {
	++$j;
    }
    

}
print  $j;
close(HANDLE_in);
return $j
}
main();
