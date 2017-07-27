#!/usr/bin/perl

# Filters a gene info file (created using the script create_gene_info_from_file.pl)
# Various filter options can be set.

# Usage: perl filter_gene_info.pl [options] file_to_filter.genfo 

use warnings;
use strict;
$|++;
use Getopt::Long;
use Cwd;
use feature 'say';

my $file_to_filter;
my $filtered_file_name;

my $chr_to_include;
my $chr_to_exclude;
my $strand_to_include;
my $min_length;
my $max_length;
my $min_GC;
my $max_GC;
my $min_transcripts;
my $max_transcripts;
my $exclude_biotypes;
my $include_biotypes;
my $list_biotypes;
my %gene_biotypes;
my $no_of_genes;
my @filtered_set;
my $dry_run;

my $help;
my $outfile_suffix = '.genfo';

my $config_result = GetOptions(
"chromosome=s" => \$chr_to_include,
"exclude_chr=s" => \$chr_to_exclude,
"strand=s" => \$strand_to_include,
"min_length=s" => \$min_length,
"max_length=s" => \$max_length,
"min_GC=s" => \$min_GC,
"max_GC=s" => \$max_GC,
"min_transcripts=s" => \$min_transcripts,
"max_transcripts=s" => \$max_transcripts,
"include_biotypes=s" => \$include_biotypes,
"exclude_biotypes=s" => \$exclude_biotypes,
"list_biotypes" => \$list_biotypes,
"output_file=s" => \$filtered_file_name,
"number_of_genes=s" => \$no_of_genes,
"dry_run" => \$dry_run,
"help" => \$help
);
die "Could not parse options" unless ($config_result);


if ($help) {
    print while (<DATA>);
    exit;
}

# some sanity checks on the options

if($max_GC && $min_GC){
	if($max_GC < $min_GC){
		print "\n max GC is smaller than min GC, no results will be generated.\n";
		exit;
	}
}

if($max_length && $min_length){
	if($max_length < $min_length){
		print "\n max length is smaller than min length, no results will be generated.\n";
		exit;
	}
}

if($max_transcripts && $min_transcripts){
	if($max_transcripts < $min_transcripts){
		print "\n max transcripts is smaller than min transcripts, no results will be generated.\n";
		exit;
	}
}


$file_to_filter = shift @ARGV;

unless($file_to_filter){

	print "gene info file is required \n";
	exit;
}	

#Check if the read file is compressed and open accordingly
if($file_to_filter =~ /\.gz$/){
	#print "\nusing file: $file_to_filter\n";
	#print "----------------------------------\n";
	open (FILE_IN,"zcat \'$file_to_filter\' |") or die "Can't read $file_to_filter: $!";
}
else{
	#print "\nusing file: $file_to_filter\n";
	#print "----------------------------------\n";
	open (FILE_IN,$file_to_filter) or die "Can't read $file_to_filter: $!";
}

if($list_biotypes){
	summarise_biotypes();
	close FILE_IN;
	exit;
}


if($filtered_file_name){
    open(OUT, '>', $filtered_file_name) or die;
}

else{  
    $filtered_file_name = $file_to_filter;
    $filtered_file_name =~ s/.txt$|.genfo$/_filtered.genfo/;
	open(OUT, '>', $filtered_file_name) or die;
}


my $line_counter = 0;

# print the header line
my $header = "gene_id\tgene_name\tchromosome\tstart\tend\tstrand\tbiotype\tbiotype_family\tlength\tGC_content\tno_of_transcripts";
print OUT "$header\n"; 

LINE_LOOP: while (my $line = <FILE_IN>){

	$line_counter++;

	# skip the header lines
	if($line =~ /^gene_id/){
		next;
	}		
	
	my @line_data = split(/\t/, $line);
	
	# all the filter options to run through
	
	if($chr_to_include){
		
		my $filt = 0;
		
		my $chr = $line_data[2];
		
		my @chrs = split(/ /,$chr_to_include);
		
		CHR_LOOP: foreach my $chrs_to_include(@chrs){
	
			if($chrs_to_include eq $chr){
				
				$filt = 1;
				# we don't want to have to cycle through each time if we've found a match
				last CHR_LOOP;
			}
		}
		if($filt==0){
			next;
		}	
	}
	
	if($chr_to_exclude){
		
	my $chr = $line_data[2];
		
		my @chrs = split(/ /,$chr_to_exclude);
		
		foreach my $chrs_to_exclude(@chrs){
		
			if($chrs_to_exclude eq $chr){
				
				next;		
			}
		}
	}
	
	if($strand_to_include){
	
		my $strand = $line_data[5];
		
		if($strand_to_include ne $strand){
		
			next;
		}
	}	
		
	if($min_length){
	
		my $length = $line_data[8];
		
		if($length < $min_length){
			next;
		}
	}
	
	if($max_length){

		my $length = $line_data[8];
	
		if($length > $max_length){
			next;
		}
	}
	
	if($min_GC){
		
		my $GC = $line_data[9];
		
		if($GC < $min_GC){
			next;
		}
	}
	
	if($max_GC){
	
		my $GC = $line_data[9];
		
		if($GC > $max_GC){
			next;
		}
	}
	
	if($min_transcripts){
		
		my $transcripts = $line_data[10];
		
		if($transcripts < $min_transcripts){
			next;
		}
	}
	
	if($max_transcripts){
		
		my $transcripts = $line_data[10];
		
		if($transcripts > $max_transcripts){
			next;
		}
	}
	
	if($exclude_biotypes){
		
	my $biotype = $line_data[6];
		
		my @biotypes = split(/ /,$exclude_biotypes);
		
		foreach my $biotypes_to_exclude(@biotypes){
		
			if($biotypes_to_exclude eq $biotype){
				
				next;		
			}
		}
	}
	
	
	if($include_biotypes){
		
		my $filt = 0;
		
		my $biotype = $line_data[6];
		
		my @biotypes = split(/ /,$include_biotypes);
		
		BIOTYPE_LOOP: foreach my $biotypes_to_include(@biotypes){
		
			if($biotypes_to_include eq $biotype){
				
				$filt = 1;
				# we don't want to have to cycle through each time if we've found a match
				last BIOTYPE_LOOP;
			}
		}
		# if there are no matches, skip this line
		if($filt==0){
			next;
		}	
	}
	
	push (@filtered_set, $line);
	
}

close FILE_IN;

if($dry_run){

	my $counter = 0;
	
	foreach my $filtered(@filtered_set){
		$counter++;
	}
	
	print "$counter genes passed the filter\n";
	
close OUT;
exit;	
}
	
	
my @chosen_set;

if($no_of_genes){
	
	# randomly select n genes from the filtered set
	for	(1 .. $no_of_genes){
		push @chosen_set, splice @filtered_set, rand @filtered_set, 1;
	}
}
else {

	@chosen_set = @filtered_set;
}		

foreach my $chosen(@chosen_set){
	print OUT $chosen;
}	
	


close OUT;	


sub summarise_biotypes{
	print "\ngene biotypes in file\n";
	print "------------------------------\n";
	
	while (my $line = <FILE_IN>){
		
		# skip the header lines
		if($line =~ /^#/){
			next;
		}	
		
		my @line_info = split(/\t/, $line);
		
		my $biotype = $line_info[6];
			
		if (exists $gene_biotypes{$biotype}){
			$gene_biotypes{$biotype} ++;
		}	
		else{
			$gene_biotypes{$biotype} = 1;
		}
	}
	
	# print the ordered hash
	foreach my $gene_biotype(sort {$gene_biotypes{$b} <=> $gene_biotypes{$a}} keys %gene_biotypes){

		print $gene_biotype, "\t", $gene_biotypes{$gene_biotype},"\n";
	}	
	exit;
}	




__DATA__

========================================================================================================================
 Perl script that filters a gene info file. Can be used for creating in silico datasets for testing purposes.
 Default is not to filter on any categories.
 
 Usage: filter_gene_info.pl [options] file_to_filter.genfo 
 
========================================================================================================================

Filter options:

 
  --chromosome          include the specified chromosomes separated by a space "1 3 X"     
  --exclude_chr         exclude the specified chromosomes
  --strand              '+', '-', filter by forward or reverse strand
  --min_length          minimum length of gene
  --max_length          maximum length of gene 
  --min_GC              minimum GC content
  --max_GC              maximum GC content
  --min_transcripts     minimum number of transcripts per gene
  --max_transcripts     maximum number of transcripts per gene
  --number_of_genes     number of genes that the output file should contain. Filtering of the file takes place first, 
                            then the specified number of genes are randomly selected. Default - all genes. 
  --exclude_biotypes    gene biotypes to exclude - string of biotypes enclosed within quotes and separated by whitespace
                            e.g. --exclude_biotypes "processed_pseudogene another_type protein_coding"
  --include_biotypes    gene biotypes to include - only these biotypes will be included. 
                            format is the same as for exclude_biotypes

Other options:

  --dry_run             doesn't create an output file, but provides a count for the number of genes that pass the filter options
  --output_file         output file name, defaults to appending input file name with _filtered.genfo
  --list_biotypes       lists all the gene biotypes along with the number of genes and exits  
  --help                display help and exit
========================================================================================================================

