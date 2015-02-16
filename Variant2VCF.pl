#!/bin/perl
# 2/6/2016 . Hugues Sicotte
# usage 
#      cat infile | perl Variant2VCF.pl [-fasta fasta] [-genome genomefile] [-bedtools full path to bedtools binary] [-annotfiles p16.gtf (can be multiple entries)] [-tmpprefix] [-out outfile.vcf] > logfile
#
# 
# THis program takes a files of CDKN2A annotation and produces a VCF
# Needs to load a transcript annotation file
#
# The input format will be
# GENE\tTRANSCRIPT\tOfficialNumbering\tOfficialProteinMutation
#
# 
# To deal with indels, 
#  need the path to bedtools,  human genome reference file (fasta) a
use strict;
print "WARNING: this is a partially incomplete program . ALL of the options are hardcoded!!!!!\n";
print "WARNING: This program does not support insertions at the splice boundary .. !!\n";
print "WARNING: Deletions that cross splice boundary will have whole intro deleted.\n";
print "WARNING: Insertions at the splice boundary will have alleles inserted in the intron\n";


my $bedtools="/projects/bsi/bictools/apps/misc/BEDTools/2.17.0/bin/bedtools";
my $fasta="/data2/bsi/reference/sequence/human/ncbi/37.1/allchr.fa";
my $hardcoded=1; # if user specifies -annotfiles, the @annotfile array should be erased and this set to $hardcoded==0, the first time a file is specified.
# Note CDS annotation does not include STOP codon..
my @annotfiles=("p14.gtf","p16.gtf");
my $tmpprefix="tmpfile"; # changing this prefix allows jobs to be parrallelized ... 

my $vcffile="p16_mutations.vcf";
my $tempbedfile="p16_mutations.temp.bed";
my $tempfasta="p16_mutations.temp.fasta.tab";
my $datestr="20150213";

for(my $i=0;$i<=$#ARGV;$i+=1) {
    if($ARGV[$i] =~ /^\-out/) {
	$i+=1;
	$vcffile=$ARGV[$i];
    }
}

open FBED,">$tempbedfile\n";


#  $bedtools  getfasta -tab -bed  bedfile  -fi $fasta -fo filename
# Transcript annotation file (tab separated bed file)
#chrom    start     stop    name score strand 
# chr15 73062668 73062770 ENST00000362710     0     -1 
#
# Need position of CDS start within exon .. and CDS end.

# gawk -F "\t" '$3=="CDS"{print $0}' p16.gtf
# Refseq file (note the geneid was modified for us to make it unique across chromosomes)
#    chr9    unknown CDS     21968231        21968241        .       -       2       gene_id "CDKN2A.chr9"; transcript_id "NM_000077"; gene_name "CDKN2A"; p_id "P8569"; tss_id "TSS12885";
# ENSEMBL file
# 
#    9       protein_coding  CDS     21974677        21974826        .       -       0       gene_id "ENSG00000147889"; transcript_id "ENST00000304494"; exon_number "1"; gene_name "CDKN2A"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "CDKN2A-001"; transcript_source "ensembl_havana"; tag "CCDS"; ccds_id "CCDS6510"; protein_id "ENSP00000307101";


my $transcripts_ref=&readGTF(@annotfiles);
my %TRANSCRIPTS=%$transcripts_ref;
my @allOK=();


my $l;

while($l=<STDIN>) {
#gene transcript cDNAmut AAmut
    chomp $l;
    my $OK=1;
    my @line=split(/\t/,$l);
    my $transcript=$line[1];
    my $cdsmut=$line[2];
    my $protmut=$line[3];
    my $relpos=0; # -1, before, 0 in CDS, 1, after CDS
    my $startpos=0;
    my $endpos=0;
    my $fromnuc="";
    my $tonuc="";
    my $type="";
    my $frame=0;
    my $goff1=0;
    my $goff2=0;
    $cdsmut =~ s/ //g;
    my $ocdsmut=$cdsmut;
    $Variant2VCF::CURRENT = $cdsmut;


    if($cdsmut =~ /^c\.(.*)/) {
	$cdsmut=$1;
    }
    if($cdsmut =~ /^\-/) {
	$relpos=-1;
    } elsif ($cdsmut =~ /^\*(.*)/) {
	$relpos=1;
	$cdsmut=$1;
    }
    if($cdsmut =~ /^([0-9,\+,\-]+)_([0-9,\+,\-]+)(.*)/) {
	# interval
	$startpos=$1;
	$endpos=$2;
	$cdsmut=$3;
    } elsif($cdsmut =~ /^([0-9,\+,\-]+)(.*)/) {
	# single point.
	$startpos=$1;
	$endpos=$1;
	$cdsmut=$2;
    }

    if($cdsmut =~ /^([A,C,G,T])>([A,C,G,T])$/) {
#CDKN2A(p14(ARF))        ENST00000304494 3G>A    M1I
#CDKN2A(p14(ARF))        ENST00000304494 283C>T  R95C
#CDKN2A(p16(INK4a))      NM_000077       329G>A  W110*
#CDKN2A(p16(INK4a))      NM_000077       -196C>T
#CDKN2A(p16(INK4a))      NM_000077       *29C>G
#CDKN2A(p16(INK4a))      NM_000077       *436A>G
	$fromnuc=$1;
	$tonuc=$2;
	$type="SNV";
	$frame=0;
    } elsif($cdsmut =~ /^del([A,C,G,T]*)$/) {
#CDKN2A(p16(INK4a))      NM_000077       238_251del      P81Cfs
#CDKN2A(p14(ARF))        ENST00000304494  404_417del      P136Lfs
#CDKN2A(p16(INK4a))      NM_000077       c.257_277delCCCGGGAGGGCTTCCTGGACA       R87_T93del
#CDKN2A(p14(ARF))        ENST00000304494  c.423_443delCCCGGGAGGGCTTCCTGGACA       P142_H148del
# the nucleotides after del are options.
#CDKN2A(p16(INK4a))      NM_000077       283del  V95Tfs
	$fromnuc=$1;
	$type="del";
	my $diffpos=&diffPos($endpos,$startpos);
	if($diffpos =~ /NA/) {
	    print "WARNING: unsupported variant on this line\n$l\n";
	    $OK=0;
	} else {
	    $frame=($diffpos+1) % 3;
	}
    } elsif($cdsmut =~ /^ins([A,C,G,T]+)$/) {
#CDKN2A(p16(INK4a))      NM_000077       280_281insC     L94Pfs
#CDKN2A(p14(ARF))        ENST00000304494  446_447insC     A149Afs
#CDKN2A(p16(INK4a))      NM_000077       *277_278insGA
	$tonuc=$1;
	$frame=length($tonuc) % 3;
	$type="ins";
    } elsif($cdsmut =~ /^del([A,C,G,T]*)ins([A,C,G,T]+)$/) {
# delins or del[A,C,G,T]ins[A,C,G,T]*
	$fromnuc=$1;
	$tonuc=$2;
	my $diffpos=&diffPos($endpos,$startpos);
	if($diffpos =~ /NA/) {
	    print "WARNING: unsupported variant on this line\n$l\n";
	    $OK=0;
	} else {
	    $frame=($diffpos+1+length($tonuc)) % 3;
	}
	$type="delins";
    } elsif($cdsmut =~ /^inv/) {
# inv  denotes inversion or inv22 (number of affected nucleotides)
	$frame=0;
	$type="inv";
    }  elsif($cdsmut =~ /^dup([A,C,G,T]*)$/) {
#CDKN2A(p16(INK4a))      NM_000077       19_22dupAGCA    S8Kfs
# the nucleotides after dup are optionsl
	$tonuc=$1;
	my $diffpos=&diffPos($endpos,$startpos);
	if($diffpos =~ /NA/) {
	    print "WARNING: unsupported variant on this line\n$l\n";
	    $OK=0;
	} else {
	    $frame=($diffpos+1) % 3;
	}
	$type="dup";
    } elsif(length($ocdsmut)>0) {
	print "WARNING: unsupported variant type ($ocdsmut) on this line\n$l\n";
	$OK=0;
    }
# $fromnuc is only for SNV.
# the $tonuc is only compulsory for SNV,ins, and delins.  For del, dup,inv we can get them from the genome.

#    print "$ocdsmut\t$startpos\t$endpos\t$type\t$relpos\t$frame\t$tonuc\n";
    if($OK==1) {
	my $bedline_ref=&mapVariant(\%TRANSCRIPTS,$startpos,$endpos,$relpos,$transcript,$ocdsmut);
	my @bedline=@$bedline_ref;
	if($bedline[2]<=$bedline[1]) {# Swap coordinates ... in proper 0-based format.
	    my $tmp=$bedline[2]-1;
	    $bedline[2]=$bedline[1]+1;
	    $bedline[1]=$tmp;
	}
	if($type eq "del") {
	    # Need an extra base for del (but not for delins)
	    $bedline[1]=$bedline[1]-1;
	}
	if($type eq "ins") {
	    # Do not need two bases for INS, just the one base before.
	    $bedline[2]=$bedline[2]-1;
	}
	# insertions already include the base before.
	print       join("\t",@bedline) . "\n";
	print  FBED join("\t",@bedline) . "\n";
	push(@bedline,$type);
	push(@bedline,$tonuc);
	push(@bedline,$protmut);
	push(@allOK,\@bedline);
    }

}

close FBED;


# unstranded bed file, with one extra base before for deletions
my $cmd="$bedtools  getfasta -tab -bed  $tempbedfile  -name -fi $fasta  -s -fo $tempfasta";
system($cmd);
open FIN,"<$tempfasta";
my %FROMNUC=();
while($l=<FIN>) {
    chomp $l;
    my @line=split(/\t/,$l);
    $FROMNUC{$line[0]}=$line[1];
}
close FIN;
open FVCF,">$vcffile";
print FVCF  "##fileformat=VCFv4.2\n";
print FVCF "##fileDate=$datestr\n";
print FVCF "##source=Variant2VCF.pl\n";
print FVCF "##reference=hg19\n";
print FVCF "##INFO=<ID=FPRED,Number=1,Type=String,Description=\"Annotated Functional Prediction\">\n";
print FVCF "##INFO=<ID=REFTRAN,Number=1,Type=String,Description=\"Reference Transcript\">\n";
print FVCF "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

for(my $i=0;$i<=$#allOK;$i+=1) {
    my $bedline_ref=$allOK[$i];
    my @bedline=@$bedline_ref;
    my $chr=$bedline[0];
    my $id=$bedline[3];
    my @idline=split(":",$id);
    my $transcript=$idline[0];
    my $fromnuc=$FROMNUC{$id};
    my $strand=$bedline[5];
    my $type=$bedline[6];
    my $tonuc=$bedline[7];
    my $fpred=$bedline[8];

    if($strand eq "-") {
	$fromnuc=&RevCompl($fromnuc);
	$tonuc=&RevCompl($tonuc);
    }

    if($type eq "ins") {
	# Need to add the base before to the "tonuc";
	if($strand eq "+") {
	    $tonuc = substr($fromnuc,0,1) . $tonuc;
	} else {
	    $tonuc = substr($fromnuc,0,1) . &RevCompl($tonuc);
	}
    } elsif($type eq "del") {
	    $tonuc = substr($fromnuc,0,1);	
    } elsif($type eq "dup") {
	$tonuc=$fromnuc . $fromnuc;
    } elsif($type eq "inv") {
	$tonuc=&RevCompl($fromnuc);
    }
# nothing to do for type=delins
    my $fpredstr="";
    if(length($fpred)>0) {
	$fpredstr=";FPRED=$fpred";
    }
    print FVCF "$chr\t" . ($bedline[1]+1) . "\t" . $bedline[3] . "\t$fromnuc\t$tonuc\t.\tPASS\tREFTRAN=\"$transcript\"" . "$fpredstr\n"; 

}
close FVCF;




    

# Note, not all known variant formats are accepted .. they are rare, just spey out a warning.


#protein consequence
# Deletion of Stop codon
#p.*110Glnext*17 (alternatively p.*110Qext*17)  p.*110Glnext*17 (alternatively p.*110Qext*17)
# in frame Del
# p.Lys2del (alternatively p.K2del) p.Lys2del (alternatively p.K2del)

# in frame Insertions
#p.*110Glnext*17 (alternatively p.*110Qext*17) p.*110Glnext*17 (alternatively p.*110Qext*17)

# frame-shifts
# p.Lys2del (alternatively p.K2del)




sub readGTF() {
# return a list with
#      HASH (per transcript) of transcript CDS exons
#             each element is a reference to a list consisting of 3 elements
#                     Strand.             
#                     reference to an exon start
#                     reference to an exon end
#                     reference to an exon number
#                     reference to a cumulative start coordinate in transcript coordinates
#                     reference to a cumulative end coordinate along transcript
#                     FASTA sequence
    my %TRANSCRIPTS=();
    my $fileid=1;
    while( my $file = shift @_) {
	$fileid+=1;
#	print "opening $file\n";
	open FIN,"<$file" or die "could not open $file\n";
	while (my $l=<FIN>){
	    chomp $l;
	    
#	chr9    unknown exon    21967751        21968241        .       -       .       gene_id "CDKN2A.chr9"; transcript_id "NM_000077"; gene_name "CDKN2A"; p_id "P8569"; tss_id "TSS12885"; exon_id "CDKN2A.NM_000077.E3";
	    my ($chr,$type,$exontype,$start,$stop,$id,$strand,$phase,$annot) = split(/\t/,$l);
	    if($exontype =~ /CDS/ || $exontype =~ /exon/) {
		my $gene="";
		my $geneid="";
		my $transcript="";
		my $exonnum=-1;
		my @elems=split(";",$annot);
		for(my $i=0;$i<=$#elems;$i+=1) {
		    my $elem=$elems[$i];
		    if($elem =~ /gene_id \"(.*)\"/) {
			$gene=$1;
		    }
		    if($elem =~ /transcript_id \"(.*)\"/) {
			$transcript=$1;
		    }
		    if($elem =~ /gene_name \"(.*)\"/) {
			$gene=$1;
		    }
#		    if($elem =~ /exon_id \"(.*)\"/) {
#			my $exonid=$1;
#			my @exonpieces=split(".",$exonid);
#			$exonnum = $exonpieces[$#exonpieces];
#		    }
		}
#		print "reading $transcript\n";
		my @tran;
		my $thisfileid=0;
		if(defined $TRANSCRIPTS{$transcript}) {
		    my $tr_ref=$TRANSCRIPTS{$transcript};
		    @tran=@$tr_ref;
		    $thisfileid=$tran[6];
		}
		if($thisfileid!=$fileid) {
		    if($thisfileid>0) {
			print "WARNING: in file $file, overwriting definition of transcript $transcript from previous file\n";
		    }
		    my @starts=();
		    my @ends=();
		    my @exonnum=();
		    my @cumstart=();
		    my @cumstop=();
		    my @phases=();
		    my @exstarts=();
		    my @exends=();
		    my $matchfirst=0;
		    my $matchlast=0;
		    my @anno = ($strand,\@starts,\@ends,\@exonnum,\@cumstart,\@cumstop,$fileid,\@phases,\@exstarts,\@exends,$transcript,$gene,$geneid,$chr,$matchfirst,$matchlast);
		    @tran=@anno;
		}
		if($exontype =~ /CDS/) {
		    my $ref=$tran[1];
		    push(@$ref,$start);
		    my $ref=$tran[2];
		    push(@$ref,$stop);
		 } elsif($exontype =~ /exon/) {
		    my $ref=$tran[8];
		    push(@$ref,$start);
		    my $ref=$tran[9];
		    push(@$ref,$stop);
		 }
#		print "read $transcript\n";
		$TRANSCRIPTS{$transcript}=\@tran;
	    }
	}
	close FIN;
	# Now reprocess each transcript and sort exons according to strand.
	foreach my $transcript (keys %TRANSCRIPTS) {
#	    print "processing $transcript\n";
	    my $tr_ref=$TRANSCRIPTS{$transcript};
	    my @tran=@$tr_ref;
	    # Only process newly loaded transcripts.
	    if($tran[6]==$fileid) {
		    # process transcripts
		    my @ordering=();
		    my @indexes=();
		    my $starts_ref=$tran[1];
		    my @starts=@$starts_ref;
		    my $ends_ref=$tran[2];
		    my @ends=@$ends_ref;
		    for(my $i=0;$i<=$#starts;$i+=1) {
			@indexes[$i]=$i;
		    }
		    if($tran[0] =~ /^\-/) {
		    @ordering=sort {$starts[$b] <=> $starts[$a]} @indexes;
		    } else {
		    @ordering=sort {$starts[$a] <=> $starts[$b]} @indexes;
		    }
		    my @exonnum=();
		    my $cumstart=1;
		    my $cumstop=0;
		    my @cumstarts=();
		    my @cumstops=();
		    my @cdsstarts=();
		    my @cdsends=();
		    my $phase=0;
		    my @phases=();
		    for(my $i=0;$i<=$#ordering;$i+=1) {
			push(@phases,$phase); # Phase Computed in previous round;
			my $len=1+$ends[$ordering[$i]]-$starts[$ordering[$i]];
			$cumstop+=$len;
			$phase = $cumstop-3*int($cumstop/3); # Will count for next round.
			push(@exonnum,$ordering[$i]+1);
			push(@cumstarts,$cumstart);
			push(@cumstops,$cumstop);
			push(@cdsstarts,$starts[$ordering[$i]]);
			push(@cdsends,$ends[$ordering[$i]]);
			$cumstart=$cumstop+1;
		    }
		    print "TRANSCRIPT $transcript CDS starts: " . join("\t",@cumstarts) . "\nTRANSCRIPT $transcript CDS ends:   " . join("\t",@cumstops) . "\n" ;

		    # Now process exons
		    my @exordering=();
		    my @exindexes=();
		    my $exstarts_ref=$tran[8];
		    my @exstarts=@$exstarts_ref;
		    my $exends_ref=$tran[9];
		    my @exends=@$exends_ref;
		    for(my $i=0;$i<=$#exstarts;$i+=1) {
			@exindexes[$i]=$i;
		    }
		    if($tran[0] =~ /^\-/) {
			@exordering=sort {$exstarts[$b] <=> $exstarts[$a]} @exindexes;
		    } else {
			@exordering=sort {$exstarts[$a] <=> $exstarts[$b]} @exindexes;
		    }
		    my @newexstarts=();
		    my @newexends=();
		    for(my $i=0;$i<=$#ordering;$i+=1) {
			push(@newexstarts,$exstarts[$exordering[$i]]);
			push(@newexends,$exends[$exordering[$i]]);
		    }

		    #
		    # To support exons before and/or after CDS, must align exons with CDS-exons
		    #   (one end will match) .. 
		    #

		    my @exonmap=();
		    my $matchfirst=-1;
		    my $matchlast=-1;
		    for(my $i=0;$i<=$#newexstarts && $matchfirst==-1;$i=$i+1) {
			if(($tran[0] eq "+" && $newexends[$i]==$cdsends[0]) ||  ($tran[0] eq "-" &&  $newexstarts[$i]==$cdsstarts[0])) {
			    $matchfirst=$i;
			}
		    }
		    for(my $i=$#newexstarts; $i>=0 && $matchlast==-1;$i=$i-1) {
			if(($tran[0] eq "+" && $newexstarts[$i]==$cdsstarts[$#cdsends]) ||  ($tran[0] eq "-" && $newexends[$i]==$cdsends[$#cdsends])) {
			    $matchlast=$i;
			}
		    }
		    die("ERROR: Error in GTF file, CDS exons do not line up with one end of exons entries for $transcript\n") unless ($matchlast!=-1 && $matchfirst!=-1);
		    
# Extend the CDS by 3 bases to account for STOP
		    if($tran[0] eq "+") {
			if($newexends[$matchlast]>= $cdsends[$#cdsends]+3){
			    $cdsends[$#cdsends]+=3;
			    $cumstops[$#cumstops]+=3;
			} else { # Add one CDS exon to extend CDS
			    my $extra=3-($newexends[$matchlast]-$cdsends[$#cdsends]);
			    $cdsends[$#cdsends]=$newexends[$matchlast];
			    $cumstops[$#cdsends]+=(3-$extra);
			    if($matchlast<$#newexends) {
				$matchlast+=1;
				my $additional_exon_start=$newexends[$matchlast];
				my $additional_exon_end= $additional_exon_start+$extra-1;
				push(@cdsstarts,$additional_exon_start);
				push(@cdsends,$additional_exon_end);
				my $cumstart=$cumstops[$#cumstops]+1;
				my $cumstop+=$cumstart+$extra-1;
				push(@cumstarts,$cumstart);
				push(@cumstops,$cumstop);
				my $phase = ($cumstart-1)-3*int(($cumstart-1)/3);
				push(@phases,$phase);
			    } else {
				print "WARNING: for transcript $transcript, had to extend (+strand) CDS by $extra bases beyond exon boundary (to account for STOP codon)\n";
				$cdsends[$#cdsends]+=3;
				$cumstops[$#cdsends]+=3;
				$newexends[$matchlast]+=3;
			    }
			}
		    }elsif($tran[0] eq "-") { 
			if($newexstarts[$matchlast]<= $cdsstarts[$#cdsstarts]-3){
			    $cdsstarts[$#cdsstarts]-=3;
			    $cumstops[$#cumstops]+=3;
			} else { # Add one CDS exon to extend CDS
			    my $extra=3-($cdsstarts[$#cdsstarts]-$newexstarts[$matchlast]);
			    $cdsstarts[$#cdsstarts]=$newexstarts[$matchlast];
			    $cumstops[$#cumstops]+=(3-$extra);
			    if($matchlast< $#newexstarts) {
				$matchlast+=1;
				my $additional_exon_end=$newexends[$matchlast];
				my $additional_exon_start= $additional_exon_end-($extra-1);
				push(@cdsstarts,$additional_exon_start);
				push(@cdsends,$additional_exon_end);
				my $cumstart=$cumstops[$#cumstops]+1;
				my $cumstop+=$cumstart+$extra-1;
				push(@cumstarts,$cumstart);
				push(@cumstops,$cumstop);
				my $phase = ($cumstart-1)-3*int(($cumstart-1)/3);
				push(@phases,$phase);
			    } else {
				print "WARNING: for trancript $transcript Must extend (-strand) CDS by $extra bases before exon boundary (to account for STOP codon)\n";
				$cdsstarts[$#cdsstarts]-=3;
				$cumstops[$#cdsends]+=3;
				$newexends[$matchlast]-=3;
			    }
			}
			
		    } else {
			die("ERROR: Invalid strand(" . $tran[0] . ") for $transcript\n");
		    }
			
#
# Pre-compute Cumulative position for "pre" (before START) and "post" (after STOP codon).
#
#
# first find out in which exon does the 5'UTR (pre) and 3'UTR (post) region are.

		    my $matchpost=$matchlast;
		    my $matchpre=$matchfirst;

		    if($tran[0] eq "+") {
			if($newexends[$matchlast]> $cdsends[$#cdsends]){
#			    $matchpost=$matchlast;
			} else {
			    if($matchlast<$#newexends) {
				$matchpost=$matchlast+1;
			    } else {
				print "WARNING: transcript $transcript ends after CDS Stop codon\n";
				$matchpost=-1;
			    }
			}
			if($newexstarts[$matchfirst]< $cdsstarts[0]){
#			    $matchpre=$matchfirst;
			} else {
			    if($matchfirst>0) {
				$matchpre=$matchfirst-1;
			    } else {
				print "WARNING: transcript $transcript has no 5' UTR\n";
				$matchpre=-1;
			    }
			}
			# Prepare looping .. for pre .. backwards ... for post .. forward
		    }elsif($tran[0] eq "-") { 
			if($newexstarts[$matchlast]< $cdsstarts[$#cdsstarts]){
#			    matchpost=$matchlast;
			} else {
			    if($matchlast<$#newexstarts) {
				$matchpost=$matchlast+1;
			    } else {
				print "WARNING: transcript $transcript has no 3'UTR\n";
				$matchpost=-1;
			    }
			}

			if($newexends[$matchfirst]>$cdsends[0]){
#              	            $matchpre=$matchfirst;
			} else { # Add one CDS exon to extend CDS
			    if($matchfirst>0) {
				$matchpre=$matchfirst-1;
			    } else {
				print "WARNING: transcript $transcript has no 5'UTR\n";
				$matchpre=-1;
			    }
			}
			# Prepare looping .. for pre .. backwards ... for post .. forward
		    }
		    
		    my $firstelem5=$matchpre;
		    my $firstelem3=$matchpost;
		    my $firstelemOff3=$cdsends[$#cdsends]-$cdsstarts[$#cdsends];
		    my $firstelemOff5=$cdsends[0]-$cdsstarts[0];


# Next, compute cumulative coordinates along the pre and post regions 

		    my @precumstarts=();
		    my @precumends=();
		    if($matchpre>=0){
			my $cumstart=1;
			my $cumstop=-$firstelemOff5;
			for(my $i=$firstelem5;$i>=0;$i=$i-1) {
			    my $len=$newexends[$i]-$newexstarts[$i];
			    $cumstop+=$len;
			    push(@precumstarts,$cumstart);
			    push(@precumends,$cumstop);
			    $cumstart=$cumstop+1;
			}
		    }

		    my @postcumstarts=();
		    my @postcumends=();
		    if($matchpost>=0){
			my $cumstart=1;
			my $cumstop=-$firstelemOff3;
			for(my $i=$firstelem3;$i<=$#newexstarts;$i=$i+1) {
			    my $len=$newexends[$i]-$newexstarts[$i];
			    $cumstop+=$len;
			    push(@postcumstarts,$cumstart);
			    push(@postcumends,$cumstop);
			    $cumstart=$cumstop+1;
			}
		    }









#		    my @anno = ($strand,\@starts,\@ends,\@exonnum,\@cumstart,\@cumstop,0,\@phases,\@exstarts,\@exends,$transcript,$gene,$geneid,$chr,$matchfirst,$matchlast);

		    @tran[1]=\@cdsstarts;
		    @tran[2]=\@cdsends;
		    @tran[3]=\@exonnum;
		    @tran[4]=\@cumstarts;
		    @tran[5]=\@cumstops;
		    @tran[6]=$fileid;
		    @tran[7]=\@phases;
		    @tran[8]=\@newexstarts;
		    @tran[9]=\@newexends;
		    @tran[14]=$matchfirst;
		    @tran[15]=$matchlast;
		    @tran[16]=\@precumstarts;
		    @tran[17]=\@precumends;
		    @tran[18]=\@postcumstarts;
		    @tran[19]=\@postcumends;
		    @tran[20]=$matchpre;
		    @tran[21]=$matchpost;
		    $TRANSCRIPTS{$transcript}=\@tran;
#		    print "processed $transcript\n";
	    }
	}
    }
    return \%TRANSCRIPTS;
}
# Do not support complex intervals overlapping the CDS and the introns.
sub diffPos(){
    my $endpos=shift @_;
    my $startpos = shift @_;
    if($startpos =~ /[\-,\+]+/ || $endpos =~ /[\-,\+]+/) {
	my $s1=$startpos;
	my $sign1=0;
	my $off1=0;
	if($startpos =~ /^([0-9]+)([\+,-])([0-9]+)/) {
	    $s1=$1;
	    $sign1=$2;
	    $off1=$3;
	}
	my $s2=$endpos;
	my $sign2=0;
	my $off2=0;
	if($endpos =~ /^([0-9]+)([\+,-])([0-9]+)/) {
	    $s2=$1;
	    $sign2=$2;
	    $off2=$3;
	}
	if($s1 == $s2 && $sign1 eq $sign2) {
	    return($off2-$off1);
	} else {
	    return "NA";
	}

    } else {
	return $endpos-$startpos;
    }
    return "NA";
}


#    print "$ocdsmut\t$startpos\t$endpos\t$type\t$relpos\t$frame\t$tonuc\n";
sub mapVariant() {
    my ($transcripts_ref,$startpos,$endpos,$relpos,$transcriptid,$cdsmut) = @_;
#    print "called mapVariant\n";
    my $tr_ref = $transcripts_ref->{$transcriptid};
    my $strand = $tr_ref->[0];
    my $chr = $tr_ref->[13];
    my $s1=$startpos;
    my $sign1=0;
    my $off1=0;
    if($startpos =~ /^([0-9]+)([\+,\-])([0-9]+)/) {
	$s1=$1;
	$sign1=$2;
	$off1=$3;
	if($sign1 eq "-") {
	    $sign1=-1;
	} else {
	    $sign1=1;
	}
#	print "Found startpos=$startpos with p=$s1,$sign1,$off1\n";
    }
#    print "calling getGenomicPos0Base with s1\n";
    my $gpos1_0base = &getGenomicPos0Base($transcripts_ref,$s1,$transcriptid,$relpos);
    my $dir1=$sign1;
    if($strand eq "-") {
	$dir1 = - $dir1;
    }
#    if($off1>0) {print "offset = " . ($dir1*$off1);}
    $gpos1_0base += ($dir1*$off1);
    my $gpos2_0base=$gpos1_0base;
    if(!($startpos eq $endpos)) {
	my $s2=$endpos;
	my $sign2=0;
	my $off2=0;
	if($endpos =~ /^([0-9]+)([\+,\-])([0-9]+)/) {
	    $s2=$1;
	    $sign2=$2;
	    $off2=$3;
	    if($sign2 eq "-") {
		$sign2=-1;
	    } else {
		$sign2=1;
	    }
#	    print "Found startpos=$startpos with p=$s1,$sign1,$off1\n";
	}
#	print "calling getGenomicPos0Base with s2\n";
	$gpos2_0base = &getGenomicPos0Base($transcripts_ref,$s2,$transcriptid,$relpos);
	my $dir2=$sign2;
	if($strand eq "-") {
	    $dir2 = - $dir2;
	}
	$gpos2_0base+=($dir2*$off2);
    }
    my @bed=($chr,$gpos1_0base,$gpos2_0base+1,$transcriptid . ":" . $cdsmut,0,$strand);
    return \@bed;

}

sub getGenomicPos0Base() {
    my $transcripts_ref=shift @_;
    my $startpos=shift @_;
    my $transcriptid = shift @_;
    my $relpos= shift @_;
    our $CURRENT;
    my $tran_ref=$transcripts_ref->{$transcriptid};
    my $strand = $tran_ref->[0];

#    print "calling getGenomicPos0Base\n";

#		    my @anno = ($strand,\@starts,\@ends,\@exonnum,\@cumstart,\@cumstop,0,\@phases,\@exstarts,\@exends,$transcript,$gene,$geneid,$chr,$matchfirst,$matchlast);
#		    @tran[1]=\@cdsstarts;
#		    @tran[2]=\@cdsends;
#		    @tran[3]=\@exonnum;
#		    @tran[4]=\@cumstarts;
#		    @tran[5]=\@cumstops;
#		    @tran[6]=$fileid;
#		    @tran[7]=\@phases;
#		    @tran[8]=\@newexstarts;
#		    @tran[9]=\@newexends;
#		    @tran[14]=$matchfirst;
#		    @tran[15]=$matchlast;
#		    @tran[16]=\@precumstarts;
#		    @tran[17]=\@precumends;
#		    @tran[18]=\@postcumstarts;
#		    @tran[19]=\@postumends;
#		    @tran[20]=$matchpre;
#		    @tran[21]=$matchpost;
#		    $TRANSCRIPTS{$transcript}=\@tran;




    my $pos;
    my $newexstarts_refs=$tran_ref->[8];
    my @newexstarts=@$newexstarts_refs;
    my $newexends_refs=$tran_ref->[9];
    my @newexends=@$newexends_refs;
    
    if($relpos==-1) {
	my $matchpre_index=$tran_ref->[20];
#	print "preindex=$matchpre_index\n";
	my $inexon=-1;
	my $i=0;
	my $cumstarts_ref=$tran_ref->[16];
	my @cumstarts=@$cumstarts_ref;
	my $cumends_ref=$tran_ref->[17];
	my @cumends=@$cumends_ref;
	my $newstartpos=-$startpos;
	if($matchpre_index>=0) {
	    while($inexon==-1 && $i<=$#cumstarts) {	    
#		print "PRE: Searching for ($startpos) in " . $cumstarts[$i] . " - " . $cumends[$i] . "\n";
		if($newstartpos>=$cumstarts[$i] && $newstartpos<=$cumends[$i]) {
#		print "PRE: FOUND         ($startpos) in " .$cumstarts[$i] . " - " . $cumends[$i]. "\n";
		    $inexon=$matchpre_index-$i;
		} else {
		    $i=$i+1;
		}
	    }
	}
	if($inexon>=0) {
# the "start" of an exon .. may still have some CDS in it, so have to count backward from end of exon.
	    if($strand eq "+") {
		my $endoffset=$cumends[$i]-$newstartpos;
		$pos=$newexstarts[$inexon]+$endoffset;
	    } elsif($strand eq "-") {
		my $endoffset=$cumends[$i]-$newstartpos;
		$pos=$newexends[$inexon]-$endoffset;
	    } else {
		die "Invalid strand: $strand\n";
	    }	    
	} else {
	    print "WARNING: For $CURRENT Position $startpos extends before transcript\n";
	    my $offset=$newstartpos-$cumends[$#cumends];
	    if($strand eq "+") {
		$pos=$newexstarts[0]-$offset;
	    } elsif($strand eq "-") {
		$pos=$newexends[0]+$offset;
	    } else {
		die "Invalid strand: $strand\n";
	    }

	}





    } elsif($relpos==1) {
	my $matchpost_index=$tran_ref->[21];
	my $inexon=-1;
	my $i=0;
	my $cumstarts_ref=$tran_ref->[18];
	my @cumstarts=@$cumstarts_ref;
	my $cumends_ref=$tran_ref->[19];
	my @cumends=@$cumends_ref;
	if($matchpost_index) {
	    while($inexon==-1 && $i<=$#cumstarts) {	    
		if($startpos>=$cumstarts[$i] && $startpos<=$cumends[$i]) {
		    $inexon=$matchpost_index-$i;
		} else {
		    $i=$i+1;
		}
	    }
	}
	if($inexon>=0) {
# the "start" of an exon .. may still have some CDS in it, so have to count backward from end of exon.
	    if($strand eq "+") {
		my $endoffset=$cumends[$i]-$startpos;
		$pos=$newexends[$inexon]-$endoffset;
	    } elsif($strand eq "-") {
		my $endoffset=$cumends[$i]-$startpos;
		$pos=$newexstarts[$inexon]+$endoffset;
	    } else {
		die "Invalid strand: $strand\n";
	    }	    
	} else {
	    print "WARNING: For $CURRENT Position $startpos extends 5' of transcript\n";
	    my $offset=$startpos-$cumends[$#cumends];
	    if($strand eq "+") {
		$pos=$newexends[$#newexends]+$offset;
	    } elsif($strand eq "-") {
		$pos=$newexstarts[$#newexstarts]-$offset;
	    } else {
		die "Invalid strand: $strand\n";
	    }

	}




    } elsif($relpos==0) {
	my $inexon=-1;
	my $i=0;
	my $cumstarts_ref=$tran_ref->[4];
	my @cumstarts=@$cumstarts_ref;
	my $cumends_ref=$tran_ref->[5];
	my @cumends=@$cumends_ref;
	while($inexon==-1 && $i<=$#cumstarts) {	    
	    # The code extends the CDS to include the STOP codon when reading the GTF.
	    if($startpos>=$cumstarts[$i] && $startpos<=$cumends[$i]) {
		$inexon=$i;
	    }
	    $i=$i+1;
	}
	if($inexon>=0) {
	    my $offset=$startpos-$cumstarts[$inexon];
	    if($strand eq "+") {
		$pos=$tran_ref->[1]->[$inexon]+$offset;
	    } elsif($strand eq "-") {
		$pos=$tran_ref->[2]->[$inexon]-$offset;
	    } else {
		die "Invalid strand: $strand\n";
	    }	    
	} else {
	    print "WARNING: For $CURRENT Position $startpos not in CDS for this reference transcript.\n"
	}
    } else {
	die "BUG: invalid relpos=$relpos at line\n$l\n";
    }

    # Since GFF positions are 1-based, must substract 1 to get a 0-based position
    return $pos-1;
}    





 sub Bed2VCFAddAlleles() {
     my ($bedline_ref,$type,$tonuc,$cdsmut) = @_;
# Get the alleles using bedtools
# return a pointer to a VCF line.
 }




#
## We probably won't need this in current version of code.
#



my %aa3to1 = ("Ala"=>"A",
"Arg"=>"R",
"Asn"=>"N",
"Asp"=>"D",
"Cys"=>"C",
"Glu"=>"E",
"Gln"=>"Q",
"Gly"=>"G",
"His"=>"H",
"Ile"=>"I",
"Leu"=>"L",
"Lys"=>"K",
"Met"=>"M",
"Phe"=>"F",
"Pro"=>"P",
"Ser"=>"S",
"Thr"=>"T",
"Trp"=>"W",
"Tyr"=>"Y",
"Val"=>"V");

my %aa1to3 = (
"A"=>"Ala",
"R"=>"Arg",
"N"=>"Asn",
"D"=>"Asp",
"C"=>"Cys",
"E"=>"Glu",
"Q"=>"Gln",
"G"=>"Gly",
"H"=>"His",
"I"=>"Ile",
"L"=>"Leu",
"K"=>"Lys",
"M"=>"Met",
"F"=>"Phe",
"P"=>"Pro",
"S"=>"Ser",
"T"=>"Thr",
"W"=>"Trp",
"Y"=>"Tyr",
    "V"=>"Val");

my %aacode = (
  "TTT" => "F", "TTC" => "F", "TTA" => "L", "TTG" => "L",
  "TCT" => "S", "TCC" => "S", "TCA" => "S", "TCG" => "S",
  "TAT" => "Y", "TAC" => "Y", "TAA" => "STOP", "TAG" => "STOP",
  "TGT" => "C", "TGC" => "C", "TGA" => "STOP", "TGG" => "W",
  "CTT" => "L", "CTC" => "L", "CTA" => "L", "CTG" => "L",
  "CCT" => "P", "CCC" => "P", "CCA" => "P", "CCG" => "P",
  "CAT" => "H", "CAC" => "H", "CAA" => "Q", "CAG" => "Q",
  "CGT" => "R", "CGC" => "R", "CGA" => "R", "CGG" => "R",
  "ATT" => "I", "ATC" => "I", "ATA" => "I", "ATG" => "M",
  "ACT" => "T", "ACC" => "T", "ACA" => "T", "ACG" => "T",
  "AAT" => "N", "AAC" => "N", "AAA" => "K", "AAG" => "K",
  "AGT" => "S", "AGC" => "S", "AGA" => "R", "AGG" => "R",
  "GTT" => "V", "GTC" => "V", "GTA" => "V", "GTG" => "V",
  "GCT" => "A", "GCC" => "A", "GCA" => "A", "GCG" => "A",
  "GAT" => "D", "GAC" => "D", "GAA" => "E", "GAG" => "E",
  "GGT" => "G", "GGC" => "G", "GGA" => "G", "GGG" => "G"
    ); 


sub RevCompl {
    my $geno = shift @_;
    my @bases = split(//,$geno);
    @bases = reverse(@bases);
    my $geno = join("",@bases);
    $geno =~ tr/A,C,T,G,N/T,G,A,C,N/;
    return $geno;
}
