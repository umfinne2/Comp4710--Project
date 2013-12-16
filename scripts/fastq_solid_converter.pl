#!/usr/bin/perl

use strict;

if (@ARGV != 1)
{
    print "DESCRIPTION:\n\tSOLiD reads have sequences enumerated between 0 - 3 rather than A, C, G, T\n";
    print "\tfor some stupid unreadable reason, so this script parses the file and converts it back\n";
    print "USAGE:\n\t$0 simulation-data/sim.fq \n";
    exit;
}

my ($file)=@ARGV;
my $count = 0;

open(my $fh, '<', $file) or die "Can't read the file '$file' [$!]\n";
while (my $line = <$fh>)
{
    chomp $line;
    
    if ($count % 4 == 1)
    {
        $line =~ s/^T//g;
        $line =~ s/0/a/g;
        $line =~ s/1/c/g;
        $line =~ s/2/g/g;
        $line =~ s/3/t/g;
    }

    print "$line\n";
    $count++;
}

