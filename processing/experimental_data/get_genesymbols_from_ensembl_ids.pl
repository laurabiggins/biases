#!/usr/bin/perl

# This was amended from the script create_gene_info_file_from_gtf.pl.
# This script requires 2 input files - one is a file with ensembl ids and associated scores.
# The other is a gtf file. The gtf file is used to find the ensembl id and access the official gene symbol.
# The official gene symbol is added to each line in the file.

use warnings;
use strict;
$|++;
use Getopt::Long;
use Cwd;
use feature 'say';

my $input_file;
my $gtf_file;
my $help;
my $outfile_suffix = '_gene_names.txt';

my $config_result = GetOptions(
"file=s" => \$input_file,
"gtf=s" => \$gtf_file
);
die "Could not parse options" unless ($config_result);

unless((defined($gtf_file)) && ($gtf_file =~ /.gtf$|.gtf.gz$/)){
	die "\ngtf file must be supplied using --file file_path\n";
}	

unless((defined($input_file)) && ($input_file =~ /.txt$|.txt.gz$/)){
	die "\ninput file must be supplied using --gtf file_path\n";
}	


#Check if the read file is compressed and open accordingly
if($gtf_file =~ /\.gz$/){
	print "\nusing gtf file: $gtf_file\n";
	print "----------------------------------\n";
	open (IN_GTF,"zcat \'$gtf_file\' |") or die "Can't read $gtf_file: $!";
}
else{
	print "\nusing gtf file: $gtf_file\n";
	#print "----------------------------------\n";
	open (IN_GTF,$gtf_file) or die "Can't read $gtf_file: $!";
}

my $outfile = $input_file;
$outfile =~ s/.txt$//;
$outfile = $outfile.$outfile_suffix;

open(OUT, '>', $outfile) or die;

my $gene_id;
my $gene_name;

my %ids;

# print the header line
my $header = "gene_name\tscore\tensembl_id";
print OUT "$header\n"; 
 
 
print "\nprocessing gtf file...\n";

# load all the gene ids and names into a hash
GTF_LOOP: while (my $line = <IN_GTF>){

	# skip the header lines
	if($line =~ /^#/){
		next;
	}	
	
	my @line_info = split(/\t/, $line);	
	
	if($line_info[2] eq "gene"){
	
		my @gene_info = split(/;/, $line_info[8]);
		
		$gene_id = $gene_info[0];
		$gene_id =~ s/gene_id "//;
		$gene_id =~ s/\"//g;
		
		$gene_name = $gene_info[2];
		$gene_name =~ s/gene_name "//;
		$gene_name =~ s/\"//g;
		$gene_name =~ s/^\s+//g;
	
		$ids{$gene_id} = $gene_name;
	}
}	

close IN_GTF;
 
#Check if the read file is compressed and open accordingly
if($input_file =~ /\.gz$/){
	print "\nusing file: $input_file\n";
	print "----------------------------------\n";
	open (IN,"zcat \'$input_file\' |") or die "Can't read $input_file: $!";
}
else{
	print "\nusing file: $input_file\n";
	#print "----------------------------------\n";
	open (IN,$input_file) or die "Can't read $input_file: $!";
}

while (my $input_line = <IN>){

	my @genes_to_annotate = split(/\t/, $input_line);
	my $ensembl_id = $genes_to_annotate[0];
	my $score = $genes_to_annotate[1];
	chomp $ensembl_id;
	chomp $score;

	if (exists $ids{$ensembl_id}) {

		my $name = $ids{$ensembl_id};
		print OUT "$name\t$score\t$ensembl_id\n";
	}
}

close IN;
close OUT;
