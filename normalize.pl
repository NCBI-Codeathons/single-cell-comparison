#!/usr/bin/env perl

use strict;
use warnings;


my $count = "/data/gene_count_length_files";
my $tpm = "/data/norm/tpm";
my $rpkm = "/data/norm/rpkm";

opendir my $ph, $count or die "can't read $count\n";
my @studies = readdir $ph;
closedir $ph;

foreach my $study(@studies) {
    next if $study =~ /^\./;
    opendir my $sd, "$count/$study" or die "can't read $count/$study\n";
    my @files = readdir $sd;
    closedir $sd;
    foreach my $file(@files) {
	next if $file =~ /^\./;
	process($study, $file);
    }
}

sub process() {
    my ($study, $fn) = @_;
    my $read_count = 0;
    my $rpk_sum = 0;
    my @lines;
    my @headers;
    open my $fh, "$count/$study/$fn" or die "can't open $count/$study/$fn: $!\n";
    while(<$fh>) {
	chop;
	if (/^\d+/) {
	    my @v = split/\t/;
	    $read_count += $v[2];
	    my $rpk = $v[2] == 0 ? 0 : $v[2]/($v[1]/1000);
	    $rpk_sum += $rpk;
	    push @v, $rpk;
	    push @lines, \@v;	    
	}
	else {
	    s/hits/hits\/norm/ if /^Gene/;
	    push @headers, $_;
	}
    }
    close $fh;
    
    my $read_scale = $read_count / 1000000;
    my $per_mil = $rpk_sum / 1000000;

    open my $rpkm_h, ">", "$rpkm/$study/$fn" or die "$rpkm/$study/$fn $!";
    open my $tpm_h, ">", "$tpm/$study/$fn" or die "$tpm/$study/$fn $!";
    
    print $rpkm_h join "\n", (@headers);
    print $rpkm_h "\n";    

    print $tpm_h join "\n", (@headers);    
    print $tpm_h "\n";    
    
    foreach my $v (@lines) {
	my @vals = @$v;
	my $rpm  = $vals[2] eq "0" ? 0 : $vals[2] / $read_scale;
	my $rpkm = $rpm eq "0"    ? 0 : sprintf("%0.3f", $rpm / ($vals[1]/1000));
	my $tpm  = $vals[3] eq "0" ? 0 : sprintf("%0.3f", $vals[3] / $per_mil);
	print $rpkm_h join "\t", (@vals[0..2], $rpkm, "\n");
	print $tpm_h join "\t", (@vals[0..2], $tpm, "\n");	
    }
    close $rpkm_h;
    close $tpm_h;
}

