#!/usr/bin/perl

use strict;

if (@ARGV != 2)
{
    print "DESCRIPTION:\n\tSOLiD error_profile_multiplier multiplies the read error values \n";
    print "\tin the supplied ART SOLiD profile file by the constant provided and prints out the result,\n";
    print "\twhich can be piped to a new file\n";
    print "USAGE:\n\t$0 SOLiD_profile/profile_default 4\n";
    exit;
}

my ($file, $mult)=@ARGV;
my $count = 0;

open(my $fh, '<', $file) or die "Can't read the file '$file' [$!]\n";
while (my $line = <$fh>)
{
    chomp $line;
    
    my @fields = split(/\s+/, $line);

    #print scalar(@fields)."\n";
    if (scalar(@fields) == 5 && $count > 1) #count prevent trying to do this on second line
    {
        #increase error columns
        my $err1 = $fields[3] * $mult;
        my $err2 = $fields[4] * $mult;

        print "$fields[0]\t$fields[1]\t$fields[2]\t$err1\t$err2\n";
    }

    else
    {
        print "$line\n";
    }

    $count++;
}

