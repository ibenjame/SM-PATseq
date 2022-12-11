#Process a minimap2 alignment file and a text file of poly-A lengths
#outputs a tabular report file joining the information at the gene assignment level with some stats
use strict;

my %h;
my ($tot, $c);

my $gencode = $ARGV[3];

#A-length scan file
open(FILE,"zcat $ARGV[1] |");
while(<FILE>){
    chomp;
    my @arr = split(/\t/,$_);
    $tot++;
    next if ($arr[1] == 0);
    $h{$arr[0]}{"length"} = $arr[1];
    $c++;
}
close(FILE);

print STDERR "$tot total reads\n$c poly-A length validated reads\n";
$c = 0;

#Minimap2 BAM file
open(FILE,"samtools view $ARGV[0] |");
while(<FILE>) {
    chomp;
    my $line = $_;
    my @arr = split(/\t/,$line);
    next unless ($line =~ /tp:A:P/);

    if($gencode eq "GENCODE"){
#Gencode gene string
        my @v = split(/\|/,$arr[2]);
        $h{$arr[0]}{"gene"} = $v[5];
    } else {
#Ensembl gene string; Cerevisiae
        $arr[2] =~ s/_mRNA.+//;
        $h{$arr[0]}{"gene"} = $arr[2];
    }

    if($line =~ /AS:i:(\d+)/) {
        $h{$arr[0]}{"max"} = $1;
    }
}
close(FILE);

#Iterators:
#read-level : bad=short/no map; nl=not long enough;
#gene-level (redundant) : ns=too short ; nc=not enough (1)
my ($g, $nl, $nc, $bad, $ns);

my %res;
foreach my $seq (keys %h) {
    if ($h{$seq}{"max"} < 50) {$bad++; next;}
    $c++;
    
    my $gene = $h{$seq}{"gene"};
    
    if ($h{$seq}{'length'} < 10) {$nl++; next;}
    
    $res{$gene}{'count'}++;
    $res{$gene}{'lens'} += $h{$seq}{'length'};
    #lengths array for stats
    push(@{$res{$gene}{"array"}}, $h{$seq}{'length'});
}

print "#Gene\tCount\tMeanPolyA\tSTDERR\tMin\tQ1\tMedian\tQ3\tMax\tListOfValues\n";
foreach my $gene (sort keys %res) {
    if ($res{$gene}{'lens'} < 10) { $ns++; next;}
    if ($res{$gene}{'count'} < 1) { $nc++; next;}
    my $avg = $res{$gene}{'lens'} / $res{$gene}{'count'};
    print "$gene\t$res{$gene}{'count'}\t";
    print "$avg";
    print "\t", std_dev($avg, @{$res{$gene}{"array"}})/sqrt($res{$gene}{'count'});
    print "\t", join("\t", stats(@{$res{$gene}{"array"}}));
    print "\t", join("\t", @{$res{$gene}{"array"}}), "\n";
    $g++;
}
print STDERR "$c total passed mapped reads\n$bad bad reads under 50bp mappable\n$nl mapped reads missing length\n$g total genes detected\n";
#print STDERR "$nc genes missing alignments\n$ns genes with too short poly-A\n";

sub std_dev {
        my ($average, @values) = @_;

        my $count = scalar @values;
        my $std_dev_sum = 0;
        $std_dev_sum += ($_ - $average) ** 2 for @values;

        return $count ? sqrt($std_dev_sum / $count) : 0;
}

sub stats
{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    my $min = $vals[0];
    my $max = $vals[-1];
    my $median;
    my ($q1, $q3);
    if($len%2) #odd?
    {
        $median = $vals[int($len/2)];
    }
    else #even
    {
        $median = ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
    if (($len/2)%2) #half odd
    {
        $q1 = $vals[int(0.25*$len)];
        $q3 = $vals[int(0.75*$len)];
    }
    else #half even
    {
        $q1 = ($vals[int(0.25*$len)] + $vals[int(0.25*$len)-1])/2;
        $q3 = ($vals[int(0.75*$len)] + $vals[int(0.75*$len)-1])/2;
    }
    return ($min, $q1, $median, $q3, $max);
}
