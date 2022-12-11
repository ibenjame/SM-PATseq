#Combine processed report files

use strict;

my @files=@ARGV;
my @ids=@files;

my %h;

foreach my $file (@files) {
    open(FILE,$file);
    while(<FILE>){
        chomp;
        my ($gid, $c, $m, $std, $min, $q1, $med, $q3, $max, @arr) = split(/\t/,$_);
        next if ($c =~ /Count/);
        $h{$gid}{$file}{'c'} = $c;
        $h{$gid}{$file}{'m'} = $m;
        $h{$gid}{$file}{'std'} = $std;
	$h{$gid}{$file}{'min'} = $min;
	$h{$gid}{$file}{'q1'} = $q1;
$h{$gid}{$file}{'med'} = $med;
	$h{$gid}{$file}{'q3'} = $q3;
	$h{$gid}{$file}{'max'} = $max;
        $h{$gid}{'samp'}++;
        $h{$gid}{'sum'} += $c;
    }
    close(FILE);
}

print "#gene\tSamples\tTotalCounts";
foreach my $file2 (@ids) {
    $file2 =~ s/.+\///;
    $file2 =~ s/\..+//;
    print "\tCount_${file2}\tMean_${file2}\tSTE_${file2}\tQ1_${file2}\tMedian_$file2\tQ3_$file2";
}
print "\n";
foreach my $gid (sort keys %h) {
    print "$gid\t$h{$gid}{'samp'}\t$h{$gid}{'sum'}";
    foreach my $file (@files) {
        print "\t", $h{$gid}{$file}{'c'} || '.', "\t", $h{$gid}{$file}{'m'} || '.', "\t", $h{$gid}{$file}{'std'} || '.';
	print "\t", $h{$gid}{$file}{'q1'} || '.', "\t", $h{$gid}{$file}{'med'} || '.', "\t", $h{$gid}{$file}{'q3'} || '.';
    }
    print "\n";
}
