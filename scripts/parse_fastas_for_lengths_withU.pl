#Scan fastas of ccs sequences produced by PB poly-A seq
#Report back poly-A or poly-T tails (inital direction based)
#Allows for terminal uridylation and related stats
#Single base interruptions to tail (bracketed by minimum 2+ A) are tolerated and rate reported

my $file=$ARGV[0];


#Inital pattern generation: Up to two mismatches allowed
my $adapter = $ARGV[1];
my $Radapter = $ARGV[2];

my @subpats;
my @subpats2;
for my $i (0..length($adapter)-1) {
   for my $j ($i+1..length($adapter)-1) {
      my $subpat = join('',
         substr($adapter, 0, $i),
         '.?',  # or '\\w'
         substr($adapter, $i+1, $j-$i-1),
         '.?',  # or '\\w'
         substr($adapter, $j+1),
      );
      my $subpat2 = join('',
         substr($Radapter, 0, $i),
         '.?',  # or '\\w'
         substr($Radapter, $i+1, $j-$i-1),
         '.?',  # or '\\w'
         substr($Radapter, $j+1),
      );
      push @subpats, $subpat;
      push @subpats2, $subpat2;
   }
}

my $pat = join('|', @subpats);
my $pat2 = join('|', @subpats2);
#print STDERR "Ambiguity Adapters:\n$pat\n$pat2\n";



#Open file and begin search
#tracking values
my ($p, $m, $b, $c) = (0, 0, 0, 0);

print STDERR "working on $file\n";
my %h;

my $id;
my $seq;

open(FILE,"gunzip -c $file |") || die "Could not open file specified (fasta required)\n";

while(<FILE>){
    chomp;
    my $line = $_;
    if ($line =~ /^>(.+)/) {
        $id = $1;
        $c++;
    } else {
        $h{$id} = "$h{$id}$line";
    }
}
print STDERR "$c sequences in file\n";
foreach my $s (keys %h) {
    my ($l, $d, $cs, $gs, $ts, $us) = (0, "-", 0, 0, 0, 0);
    #l = length, d = direction, cs = # C substitutions, gs, ts = G and T subs. us = uridyls
    my @res = $h{$s} =~ /($pat2|GCAGAG)(A*TT+(.TT+)*)/;
    if (@res > 0) {
        my $string2 = pop(@res);
        my $string = pop(@res);
        if ($string =~ /^(A+)/) {
		$us = length($1);
		$string =~ s/^A+//;
	}
        while($string =~ m/G/g) { $cs++;}
        while($string =~ m/C/g) { $gs++;}
        while($string =~ m/A/g) { $ts++;}
        my $a;
	    while ($h{$s} =~ /(T+)/g) {
		    my $f = length($1);
		    if ($f > $a) {$a = $f;}
	    }
	    if ($a > $l) {$l += $a;}
	    $d = "-";
    }
    if ($h{$s} =~ /((AA+.)*AA+T*)($pat|CTCTGC)/) {
        my $string = $1;
	my ($cr, $gr, $tr, $ur) = (0,0,0,0);
	if ($string =~ /(T+)$/) {
		$ur = length($1);
		$string =~ s/T+$//;
	}
        my $l2 = length($string);
        while($string =~ m/G/g) { $gr++;}
        while($string =~ m/C/g) { $cr++;}
        while($string =~ m/T/g) { $tr++;}
        my $a;
	    while ($h{$s} =~ /(A+)/g) {
                my $f = length($1);
                if ($f > $a) {$a = $f;}
        }
        if ($a > $l2) {$l2 += $a;}
        if ($l2 > $l) {
            $l = $l2;
            $d = "+";
	    ($cs, $gs, $ts, $us) = ($cr, $gr, $tr, $ur);
        }
    }
    if ($l < 12) { $b++; print "$s\t0\tFail\n";}
    elsif ($d eq "+") {$p++; print "$s\t$l\t+\t$cs\t$gs\t$ts\t$us\n";}
    else { $m++; print "$s\t$l\t-\t$cs\t$gs\t$ts\t$us\n";}
    
}
close(FILE);

print STDERR "$p forward CCS with poly-A\n$m reverse CCS with poly-T\n$b unresolved seq\n";
