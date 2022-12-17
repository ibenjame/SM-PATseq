#Process a blat alignment file and a text file of poly-A lengths
#Provides an extended per-read tabular report file
use strict;

my $gencode = $ARGV[2];

my %h;
my ($tot, $c);

#Minimap2 BAM file
open(FILE,"samtools view $ARGV[0] |");
while(<FILE>) {
    chomp;
    my $line = $_;
    my @arr = split(/\t/,$line);
    next unless ($line =~ /tp:A:P/);

    $h{$arr[0]}{"full"}=$arr[2];
    if($gencode eq "GENCODE"){
#Gencode
        my @v = split(/\|/,$arr[2]);
        $h{$arr[0]}{"gene"} = $v[5];
    } else {
#Ensembl; Cerevisiae
        $arr[2] =~ s/_mRNA.+//;
        $h{$arr[0]}{"gene"} = $arr[2];
    }
    if($line =~ /AS:i:(\d+)/) {
        $h{$arr[0]}{"max"} = $1;
    }
}
close(FILE);

#A-length scan file
open(FILE,"gunzip -c $ARGV[1] |");
while(<FILE>){
    chomp;
    my @arr = split(/\t/,$_);
    $tot++;
    $h{$arr[0]}{"length"} = $arr[1];
    if ($arr[1] > 0) {$c++;} else {$arr[3] = 0; $arr[4] = 0; $arr[5] = 0; $arr[6] = 0;}
    $h{$arr[0]}{"gene"} = "NA" unless exists $h{$arr[0]}{"gene"};
    $h{$arr[0]}{"full"} = "NA" unless exists $h{$arr[0]}{"full"};
    print join("\t", @arr, $h{$arr[0]}{"gene"}, $h{$arr[0]}{"max"}+0, $h{$arr[0]}{"full"}), "\n"; 
}
close(FILE);

print STDERR "$tot total reads\n$c poly-A length validated reads\n";
