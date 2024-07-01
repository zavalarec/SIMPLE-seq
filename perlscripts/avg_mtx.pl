#!/usr/bin/perl
use strict;
use warnings;

my %data;
my %peaks;
my %peaks_count;
my %peaks_total;

my $input = $ARGV[0];
my $cutoff = $ARGV[1];
my $output = $ARGV[0]."_avg_".$cutoff;


system("mkdir $output");
system("cp $input/barcodes.tsv $output/");

open IN, "$input/genes.tsv" or die $!;
my $idx = 0;
while(<IN>){
        $idx++;
        chomp;
        $peaks{$idx} = $_;
}
close IN;
open IN, "$input/matrix.mtx" or die $!;
my $read_title = 0;
my @title;
while(<IN>){
        chomp;
        next if $_ =~ m/^\%/;
        if($read_title == 0){
                $read_title = 1;
                @title = split/\s+/, $_;
                next;
        }
        my @tmp = split/\s+/, $_;
        my $peak_idx = $tmp[0];
        my $cell_idx = $tmp[1];
        my $value = $tmp[2];
        $data{$peak_idx}{$cell_idx} = $value;
        $peaks_count{$peak_idx} = 0 if not exists $peaks_count{$peak_idx};
        $peaks_count{$peak_idx} ++;
        $peaks_total{$peak_idx} = 0 if not exists $peaks_total{$peak_idx};
        $peaks_total{$peak_idx}+= $value;
}
close IN;
$idx = 0;
my %new_peaks;
open OUT, ">$output/genes.tsv" or die $!;
foreach my $peak_idx (sort{$a<=>$b} keys %peaks_count){
        next if $peaks_count{$peak_idx} < $title[1]*$cutoff;
        $idx++;
        $new_peaks{$idx} = $peak_idx;
        print OUT $peaks{$peak_idx}."\n";;
}
close OUT;

my $n_total_peaks = $idx;
my $n_total_features = $n_total_peaks * $title[1];

open OUT, ">$output/matrix.mtx" or die $!;
print OUT "\%\%MatrixMarket matrix coordinate real general\n";
print OUT "\%\n";
print OUT "$n_total_peaks $title[1] $n_total_features\n";
foreach my $i (1 .. $n_total_peaks){
        my $old_peak_idx = $new_peaks{$idx};
        my $avg = $peaks_total{$old_peak_idx} / $peaks_count{$old_peak_idx};
        foreach my $j (1 .. $title[1]){
                my $value = $avg;
                $value = $data{$i}{$j} if exists $data{$i}{$j};
                print OUT "$i $j $value\n";
        }
}
close OUT
