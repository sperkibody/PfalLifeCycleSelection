use warnings;
use strict;

# Assumed fasta heading structure:
# >PF3D7_0100100.1 | organism=Plasmodium_falciparum_3D7 | product=erythrocyte membrane protein 1, PfEMP1 | location=Pf3D7_01_v3:29510-37126(+) | length=7617 | sequence_SO=chromosome | SO=protein_coding_gene

my $input = "<fasta>"; ## 1 entry per transcript (coding only) with above format
my $output = "<output.txt>";

## Input coding fasta file
open OUT, ">$output";
print OUT "GENE\tTRANS\tNAME\tNS\tSYN\tFFD\tTOTAL_CODING_LENGTH\tPROP_NS\tPROP_SYN\tPROP_FFD\tSTOP_CODONS\tFULL_LENGTH\tCOORD\n";
open CODING, "$input" or die;
my $transcript;
my @full_sequence;
my $length;
my $gene;
my $name;
my $coord;
my $codon_count;

while (<CODING>) {
  chomp;
  $_ =~ s/^\s*(.*?)\s*$/$1/;
  if ($_ =~ ">.*") {
    if (defined($full_sequence[0])) {
##### Count SYN sites
      my $syn_sites = 0;
      my $codon_count = 0;
      my $ffd_sites = 0;
      my $stop_codons = 0;
      while (defined($full_sequence[0])) {
      	++ $codon_count;
      	my $c1 = shift(@full_sequence);
		my $c2 = shift(@full_sequence);
		my $c3 = shift(@full_sequence);
		my $codon = "$c1"."$c2"."$c3";
		my $new_syn = &countSyn($codon);
		my $new_ffd = &countFFD($codon);
		if (($codon eq 'TAG') || ($codon eq 'TAA') || ($codon eq 'TGA')) {
			print "Stop codon ($codon) at $codon_count.\n";
			++$stop_codons
		}
		$syn_sites = ($syn_sites + $new_syn);
		$ffd_sites = ($ffd_sites + $new_ffd);
      }
      my $codon_length= ($codon_count * 3);
      my $non_syn = ($codon_length-$syn_sites);
      if ($codon_length != $length){
      	print "fasta length of $length but $codon_count codons ($codon_length nt)\n";
      }
      my $prop_ns = $non_syn/$codon_length;
      my $prop_syn = $syn_sites/$codon_length;
      my $prop_ffd = $ffd_sites/$codon_length;
      print OUT "$gene\t$transcript\t$name\t$non_syn\t$syn_sites\t$ffd_sites\t$codon_length\t$prop_ns\t$prop_syn\t$prop_ffd\t$stop_codons\t$length\t$coord\n";
    }
    @full_sequence = ();
    $length = 0;
    my @info = split(/ \| /, $_);
    my @coord_info;
    my @pos = ();
    my $complement=0;
    $coord = $info[3];
    if ($coord =~ '\(\-\)') {
      $complement = 1;
    }
    $transcript = $info[0];
    $transcript =~ s/>//g;
    print "Starting next transcript: $transcript. On strand $complement\n";
    $name = $info[2];
    $name =~ s/product=//g;
    ($gene, my $transID) = split(/\./, $transcript);
    $length = $info[4];
    $length =~ s/length=//g;   
    next;
  }
  my @dl = split("",$_);
  push(@full_sequence, @dl);
}
close CODING;

###For Final Transcript####

## Count SYN sites
my $syn_sites = 0;
my $ffd_sites = 0;
my $stop_codons = 0;
while (defined($full_sequence[0])) {
  ++ $codon_count;
  my $c1 = shift(@full_sequence);
  my $c2 = shift(@full_sequence);
  my $c3 = shift(@full_sequence);
  my $codon = "$c1"."$c2"."$c3";
  my $new_syn = &countSyn($codon);
  my $new_ffd = &countFFD($codon);
  if (($codon eq 'TAG') || ($codon eq 'TAA') || ($codon eq 'TGA')) {
    print "Stop codon ($codon) at $codon_count.\n";
    ++$stop_codons;
  }
  $syn_sites = ($syn_sites + $new_syn);
  $ffd_sites = ($ffd_sites + $new_ffd);
}
my $codon_length= ($codon_count * 3);
my $non_syn = ($codon_length-$syn_sites);
if ($codon_length != $length){
	print "fasta length of $length but $codon_count codons ($codon_length nt)\n";
}
my $prop_ns = $non_syn/$codon_length;
my $prop_syn = $syn_sites/$codon_length;
my $prop_ffd = $ffd_sites/$codon_length;
print OUT "$gene\t$transcript\t$name\t$non_syn\t$syn_sites\t$ffd_sites\t$codon_length\t$prop_ns\t$prop_syn\t$prop_ffd\t$stop_codons\t$length\t$coord\n";

close OUT;

## Subroutine based on Polymorphorama
sub countSyn {
  my($codon_look) = $_[0];
#  print "\n$codon_look";
  $codon_look = uc $codon_look;
  if (($codon_look eq 'TAG') || ($codon_look eq 'TAA') || ($codon_look eq 'TGA')) {
    print "Reached stop codon: $codon_look\n";
    return 0;
  }
  
  my(%genetic_code) = (
		       
		       'TCA' => 1,    # Serine
		       'TCC' => 1,    # Serine
		       'TCG' => 1,    # Serine
		       'TCT' => 1,    # Serine
		       'TTC' => (1/3),    # Phenylalanine
		       'TTT' => (1/3),    # Phenylalanine
		       'TTA' => (2/3),    # Leucine
		       'TTG' => (2/3),    # Leucine
		       'TAC' => 1,    # Tyrosine
		       'TAT' => 1,    # Tyrosine
		       #'TAA' => '_',    # Stop
		       #'TAG' => '_',    # Stop
		       'TGC' => 0.5,    # Cysteine
		       'TGT' => 0.5,    # Cysteine
		       #'TGA' => '_',    # Stop
		       'TGG' => 0,    # Tryptophan
		       'CTA' => (4/3),    # Leucine
		       'CTC' => 1,    # Leucine
		       'CTG' => (4/3),    # Leucine
		       'CTT' => 1,    # Leucine
		       'CCA' => 1,    # Proline
		       'CCC' => 1,    # Proline
		       'CCG' => 1,    # Proline
		       'CCT' => 1,    # Proline
		       'CAC' => (1/3),    # Histidine
		       'CAT' => (1/3),    # Histidine
		       'CAA' => (1/3),    # Glutamine
		       'CAG' => (1/3),    # Glutamine
		       'CGA' => (4/3),    # Arginine
		       'CGC' => 1,			# Arginine
		       'CGG' => (4/3),    # Arginine
		       'CGT' => 1,		# Arginine
		       'ATA' => (2/3),    # Isoleucine
		       'ATC' => (2/3),    # Isoleucine
		       'ATT' => (2/3),    # Isoleucine
		       'ATG' => 0,    # Methionine
		       'ACA' => 1,    # Threonine
		       'ACC' => 1,    # Threonine
		       'ACG' => 1,    # Threonine
		       'ACT' => 1,    # Threonine
		       'AAC' => (1/3),    # Asparagine
		       'AAT' => (1/3),    # Asparagine
		       'AAA' => (1/3),    # Lysine
		       'AAG' => (1/3),    # Lysine
		       'AGC' => (1/3),    # Serine
		       'AGT' => (1/3),    # Serine
		       'AGA' => (5/6),    # Arginine
		       'AGG' => (2/3),    # Arginine
		       'GTA' => 1,    # Valine
		       'GTC' => 1,    # Valine
		       'GTG' => 1,    # Valine
		       'GTT' => 1,    # Valine
		       'GCA' => 1,    # Alanine
		       'GCC' => 1,    # Alanine
		       'GCG' => 1,    # Alanine
		       'GCT' => 1,    # Alanine
		       'GAC' => (1/3),    # Aspartic Acid
		       'GAT' => (1/3),    # Aspartic Acid
		       'GAA' => (1/3),    # Glutamic Acid
		       'GAG' => (1/3),    # Glutamic Acid
		       'GGA' => 1,    # Glycine
		       'GGC' => 1,    # Glycine
		       'GGG' => 1,    # Glycine
		       'GGT' => 1,    # Glycine
		      );
  
  if(exists $genetic_code{$codon_look}) {
    return $genetic_code{$codon_look};
  } else {
    print "ERROR with codon: $codon_look\n"
  }
}

sub countFFD {
  my($codon_look) = $_[0];
#  print "\n$codon_look";
  $codon_look = uc $codon_look;
  if (($codon_look eq 'TAG') || ($codon_look eq 'TAA') || ($codon_look eq 'TGA')) {
    print "Reached stop codon: $codon_look\n";
    return 0;
  }
  
  my(%four_fold_deg) = (
		       
		       'TCA' => 1,    # Serine
		       'TCC' => 1,    # Serine
		       'TCG' => 1,    # Serine
		       'TCT' => 1,    # Serine
		       'TTC' => 0,    # Phenylalanine
		       'TTT' => 0,    # Phenylalanine
		       'TTA' => 0,    # Leucine
		       'TTG' => 0,    # Leucine
		       'TAC' => 0,    # Tyrosine
		       'TAT' => 0,    # Tyrosine
		       #'TAA' => '_', # Stop
		       #'TAG' => '_', # Stop
		       'TGC' => 0,    # Cysteine
		       'TGT' => 0,    # Cysteine
		       #'TGA' => '_', # Stop
		       'TGG' => 0,    # Tryptophan
		       'CTA' => 1,    # Leucine
		       'CTC' => 1,    # Leucine
		       'CTG' => 1,    # Leucine
		       'CTT' => 1,    # Leucine
		       'CCA' => 1,    # Proline
		       'CCC' => 1,    # Proline
		       'CCG' => 1,    # Proline
		       'CCT' => 1,    # Proline
		       'CAC' => 0,    # Histidine
		       'CAT' => 0,    # Histidine
		       'CAA' => 0,    # Glutamine
		       'CAG' => 0,    # Glutamine
		       'CGA' => 1,    # Arginine
		       'CGC' => 1,    # Arginine
		       'CGG' => 1,    # Arginine
		       'CGT' => 1,    # Arginine
		       'ATA' => 0,    # Isoleucine
		       'ATC' => 0,    # Isoleucine
		       'ATT' => 0,    # Isoleucine
		       'ATG' => 0,    # Methionine
		       'ACA' => 1,    # Threonine
		       'ACC' => 1,    # Threonine
		       'ACG' => 1,    # Threonine
		       'ACT' => 1,    # Threonine
		       'AAC' => 0,    # Asparagine
		       'AAT' => 0,    # Asparagine
		       'AAA' => 0,    # Lysine
		       'AAG' => 0,    # Lysine
		       'AGC' => 0,    # Serine
		       'AGT' => 0,    # Serine
		       'AGA' => 0,    # Arginine
		       'AGG' => 0,    # Arginine
		       'GTA' => 1,    # Valine
		       'GTC' => 1,    # Valine
		       'GTG' => 1,    # Valine
		       'GTT' => 1,    # Valine
		       'GCA' => 1,    # Alanine
		       'GCC' => 1,    # Alanine
		       'GCG' => 1,    # Alanine
		       'GCT' => 1,    # Alanine
		       'GAC' => 0,    # Aspartic Acid
		       'GAT' => 0,    # Aspartic Acid
		       'GAA' => 0,    # Glutamic Acid
		       'GAG' => 0,    # Glutamic Acid
		       'GGA' => 1,    # Glycine
		       'GGC' => 1,    # Glycine
		       'GGG' => 1,    # Glycine
		       'GGT' => 1,    # Glycine
		      );
  
  if(exists $four_fold_deg{$codon_look}) {
    return $four_fold_deg{$codon_look};
  } else {
    print "ERROR with codon: $codon_look\n"
  }
}
