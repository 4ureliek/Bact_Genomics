#!/usr/bin/perl -w
#######################################################
# Author  :  Aurelie Kapusta
# email   :  4urelie.k@gmail.com
#######################################################
use warnings;
use strict;
use Carp;
use Getopt::Long;
use Text::CSV;
use Bio::Seq;
#use Data::Dumper;

my $VERSION = "2.2";
my $SCRIPTNAME = "Prokka_Roary_extract.pl";
my $CHANGELOG = "
#  v1.0 = 27 Sep 2018
#  v1.1 = 20 Nov 2018
#		    Skip the header line of the csv
#  v1.2 = 28 Feb 2019
#         Allow prefix
#  v2.0 = 02 Jul 2019
#         Modify strategy to make this much faster
#  v2.1 = 18 Jul 2019
#         Load descriptions in lowercases
#         Add strand
#         Remove the flanking option here, could be too messy
#  v2.2 = 19 Sept 2019
#         Bug fix in group printing
\n";

my $USAGE = "\nUsage [$VERSION]: 
    Usage:
    perl $SCRIPTNAME -g <path_to_prokka_gff> [-r <gene_presence_absence.csv>]
    [-x <X>] [-d] [-p] [-a] [-m <path>] [-l] [-h] [-v]
	
    Reads the roary gene_presence_absence.csv and extract sequences
    Will skip any lines starting with #
	    
    MANDATORY ARGUMENTS:	
    -g,--gff   (STRING) location of the prokka outputs (path to folder)       
                            
    OPTIONAL ARGUMENTS:
    -r,--roary (STRING) gene_presence_absence.csv file from roary
    -x,--x     (STRING) if prefix in front of the sample IDs (e.g. HFLUA_XX-XXX)
    -d,--desc    (BOOL) To concatenate sequences based on description 
                        (skips hypotheticals)
    -p,--prot    (BOOL) To translate the sequence
    -a,--aln     (BOOL) to align the sequences
    -m,--m     (STRING) path to muscle alignment software
    -l,--log     (BOOL) Print the change log (updates)
    -v           (BOOL) Verbose mode, make the script talks to you
    -v           (BOOL) Print version if only option
    -h,--help    (BOOL) Print this usage\n\n";        
        
#-------------------------------------------------------------------------------
#------------------------------ LOAD AND CHECK ---------------------------------
#-------------------------------------------------------------------------------
my ($GFF,$DESC,$ALIGN,$PROT,$HELP,$V,$CHLOG);
my $ROARY;
my $PREFIX = "";
my $MUSCLE = "/home/akapusta/muscle3.8.31/muscle3.8.31";
#for staph:
#my $PATTERN = "Luk.faa";
#for Spneumo:
my $PATTERN = "unique.annot.fa";
GetOptions ('g=s'     => \$GFF, 
            'r=s'     => \$ROARY,
            'd'       => \$DESC,
            'm=s'     => \$MUSCLE,
            'p'       => \$PROT,
            'x=s'     => \$PREFIX,
            'a'       => \$ALIGN, 
            'l'       => \$CHLOG, 
            'h'       => \$HELP, 
            'v'       => \$V);

#check step to see if mandatory argument is provided + if help/changelog
die "\n Script $SCRIPTNAME version $VERSION\n\n" if (! $GFF && ! $ROARY && ! $HELP && ! $CHLOG && $V);
die $CHANGELOG if ($CHLOG);
die $USAGE if ((! $GFF && ! $ROARY) || $HELP);
die "\n -g $GFF does not exist?\n\n" if (! -e $GFF);
die "\n -r $ROARY does not exist?\n\n" if ($ROARY && ! -e $ROARY);
die "\n -m $MUSCLE does not exist?\n\n" if ($ALIGN && ! -e $MUSCLE);
$GFF =~ s/\/$//;

#-------------------------------------------------------------------------------
#----------------------------------- MAIN --------------------------------------
#-------------------------------------------------------------------------------
print STDERR "\n --- Script $SCRIPTNAME started (v$VERSION)\n" if ($V);

#First extract the gene sequences for prokka
print STDERR " --- Unless they exist, extracting gene sequences from the gff files in dir: $GFF\n" if ($V);
my @FILES = ();
my $GOUT = $GFF.".fa";
my $EOUT = $ROARY.".extract";
mkdir($GOUT) unless (-e $GOUT);
get_files();
extract_gff_gene_seqs();
print STDERR "     ...done\n" if ($V);

my %SEQS = ();
my %ROARY_NAME = ();
my %TOALN = ();
if ($ROARY) {
	print STDERR " --- loading groups from $ROARY\n" if ($V);
	load_roary_csv();
	print STDERR "     ...done\n" if ($V);

	print STDERR " --- Now grouping sequences based on group IDs\n" if ($V);
	print STDERR "     And descriptions\n" if ($DESC && $V);
	`rm -Rf $EOUT` if (-e $EOUT);
	mkdir($EOUT);
	extract_sequences();
	print STDERR "     ...done\n" if ($V);

	if ($ALIGN) {
		print STDERR " --- Aligning sequences\n" if ($V);
		align_sequences();
	}
}

print STDERR " --- Script done\n\n" if ($V);
exit;

#-------------------------------------------------------------------------------
#------------------------------- SUBROUTINES -----------------------------------
#-------------------------------------------------------------------------------
sub get_files {
	opendir (my $dir, $GFF) or confess "     \nERROR (get_files): could not open to read the directory $GFF $!\n";
	@FILES = grep { /\.gff$/ && -f "$GFF/$_" } readdir($dir);
	closedir $dir;
	
	if (! $FILES[0]) {
		print STDERR "     \nERROR: no file in $GFF?\n\n" if ($V);
		exit;
	}
	return 1;
}

#-------------------------------------------------------------------------------
sub get_fa_file {
	my $gff = shift;
	my $ln = `grep -n "##FASTA" $GFF/$gff`;
	chomp $ln;
	$ln =~ s/:.*$//;
	my $tot = `wc -l $GFF/$gff`;
	chomp $tot;
	$tot =~ s/([0-9]+?)\s.*/$1/;
	$tot =~ s/\s+//;
	my $faln = $tot - $ln +1;
	`tail -n $faln $GFF/$gff | perl -pe '/^>/ ? print "\n" : chomp' > $GOUT/$gff.fa`;
	#this makes a blank line at the top, no big deal, but ultimately wanna fix it
	return 1;
}

#-------------------------------------------------------------------------------
sub extract_gff_gene_seqs {
	foreach my $gff (@FILES) {
		get_fa_file($gff) unless (-e "$GOUT/$gff.fa");
		my $gfa = $GOUT."/".$gff.".genes.fa";
		$gfa = $gfa."a" if ($PROT);
		next if (-e $gfa);
		open (my $fhi, "<", $GFF."/".$gff) or confess "ERROR (sub extract_gff_gene_seqs): could not open to read $GFF/$gff $!\n";
		open (my $fho, ">", $gfa) or confess "ERROR (sub extract_gff_gene_seqs): could not open to write $gfa $!\n";
		while(defined(my $l = <$fhi>)) {
			chomp($l);
			last if ($l =~ /^##FASTA/);
			next if (substr($l,0,1) eq "#");
			my @l = split(/\t/,$l);
			next if ($l[2] ne "CDS");
			my ($chr,$st,$en) = ($l[0],$l[3],$l[4]);
			my ($strand,$info) = ($l[6],$l[8]);
			my @info = split(";",$info);
			my %desc = ();
			foreach my $vals (@info) {
				my @v = split("=",$vals);
				$desc{$v[0]}=$v[1];
			}
				#ID=10754X1_00001
				#inference=ab initio prediction:Prodigal:2.6,similar to AA sequence:GenBank_GenPept.20181108.unique.annot.fa:ADM90119
				#locus_tag=10754X1_00001
				#product=competence-specific global transcription modulator
				
			$desc{'inference'} =~ s/^.*$PATTERN:(.+?)/$1/; 
			$desc{'product'} = "na" unless ($desc{'product'});
			my $header = ">$desc{'ID'}\tstr=$strand;$desc{'product'};$desc{'inference'}";
			#now extract
			my $seq = `grep -A 1 "$chr" $GOUT/$gff.fa | tail -n 1`;
			chomp($seq);		
			my $totlen = length($seq);
			$st = $st-1; #-1 because substr is base 0
			my $len = $en-$st; #NOT +1 here - it was adding an extra nt... issue for antisense!
			my $subseq=substr($seq,$st,$len);	
			my $prot_obj;
			if ($strand eq "-" ) {
				my $seqio = Bio::Seq->new(-alphabet => 'dna', -seq => $subseq);
				my $revio = $seqio->revcom;
				$subseq = $revio->seq;	
			}
			if ($PROT) {
				my $seqio = Bio::Seq->new(-alphabet => 'dna', -seq => $subseq);
				$prot_obj = $seqio->translate;
				$subseq = $prot_obj->seq;
			}		
			print $fho "$header\n$subseq\n";		
		}
		close $fhi;
		close $fho;
	}	
	return 1;
}

#-------------------------------------------------------------------------------
sub load_roary_csv {
	#load the csv into array
	my $csv = Text::CSV->new ( { binary => 1 } )  # should set binary attribute.
						 or confess "ERROR (sub load_roary_csv): cannot use CSV: ".Text::CSV->error_diag ();
	
	open my $fh, "<:encoding(utf8)", $ROARY or confess "$ROARY: $!";
	#now loop through roaray csv - first 14 cols are details about the genes:
	# 	0		1						2			3				4				5							6					
	# 	Gene	Non-unique Gene name	Annotation	No. isolates	No. sequences	Avg sequences per isolate	Genome Fragment	
	#   7						8					9								10	11					12					13			
	#   Order within Fragment	Accessory Fragment	Accessory Order with Fragment	QC	Min group size nuc	Max group size nuc	Avg group size nuc
	while ( my $row = $csv->getline( $fh ) ) {
		#check if relevant line
		my $grp = $row->[0];
		next if $grp eq "Gene";
		my $desc = $row->[2];	
		my $name = $row->[1];
		$name = "na" if (! defined $name || $name eq "");
		
		#now get the exact gene ID for each sample so I can filter after
		my @row = @{$row};
		my $nb = scalar(@row);
		my @s = @row[14 .. $nb];
		my $i = 0;
		for ($i = 0; $i < scalar(@s)-1; $i++) {
			next if ($s[$i] !~ /\w/);
			my @ids = split(/\s+/,$s[$i]);
			foreach my $id (@ids) {
#				my $id = $1 if ($s[$i] =~ /^(.+?)_.*/);
				my $rdesc = lc($desc);
				if ($SEQS{$rdesc}{$grp}) {
					$SEQS{$rdesc}{$grp}.=",$id";
				} else {	
					$SEQS{$rdesc}{$grp}=$s[$i];
				}
				$ROARY_NAME{$id}=$name;
			}	
		}
	}
	$csv->eof or $csv->error_diag();
	close $fh;	
	return 1;
}

#-------------------------------------------------------------------------------
sub extract_sequences {
	DESC: foreach my $desc (keys %SEQS) {
		my $out;
		my $fh;
		if ($DESC) {
			next DESC if ($desc =~ /[Hh]ypothetical/);
			my $fdesc = $desc;
			$fdesc =~ s/\s+/_/g;
			$fdesc =~ s/\(/_/g;
			$fdesc =~ s/\)/_/g;
			$fdesc =~ s/\//_/g;
			$out = $out = $EOUT."/".$fdesc.".fa";
			open ($fh, ">", $out) or confess "ERROR (sub extract_sequences): could not open to write $out $!\n";
		}
		foreach my $grp (keys %{$SEQS{$desc}}) {
			my @ids = split(",",$SEQS{$desc}{$grp});
			for my $id (@ids) {
				my $sample = $1 if ($id =~ /^($PREFIX.+?)_.*/);
#				$sample =~ s/\.[1-2]//;
				my $gff = $GOUT."/".$sample.".gff.genes.fa";
				if (! $DESC) { 
					$out = $EOUT."/".$grp.".fa";
					open ($fh, ">>", $out) or confess "ERROR (sub extract_sequences): could not open to write $out $!\n";
				}	
				my @list = ();
				if ($id =~ /\s/) {
					@list = split(/\s+/,$id);
				} else {
					$list[0]=$id;
				}				
				for (my $i=0;$i<scalar(@list);$i++) {
					my $header = `grep $list[$i] $gff`;
					chomp $header;
					my ($iid,$idesc) = split(/\t/,$header);
 					my $seq = `grep -A 1 $list[$i] $gff | grep -v ">"`; 
					chomp $seq;
					print $fh ">$grp##$list[$i]\tPROKKA_DESC=$idesc\tROARY_DESC=$desc\tROARY_NAME=$ROARY_NAME{$list[$i]}\n$seq\n";
 				}
 				close $fh if (! $DESC);
 				$TOALN{$out}=1 if ($ALIGN);
			}
		}
		close $fh if ($DESC);
	}
	return 1;
}

#-------------------------------------------------------------------------------
sub align_sequences {
	foreach my $fa (keys %TOALN) {
		my $aln = $fa;
		$aln =~ s/\.fa/\.aln.fa/;
		`$MUSCLE -in $fa -out $aln -log $aln.log -quiet -verbose`;
	}
	return 1;
}





