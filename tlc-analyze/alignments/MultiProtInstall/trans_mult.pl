#!/usr/bin/perl -w

#Ver 1.5


use strict;
use FindBin;

my $home="$FindBin::Bin";
my $suffix=".Linux";

# #detect the machine type
# my $mtype=`uname -s`;

# if($mtype =~ /IRIX/){
#     $suffix= ".IRIX";
# }else{
#     if($mtype =~ /Linux/){
#         $suffix= ".Linux";
#     }else{
#         print "This type of machine: $mtype";
#         print "is not supported.\n";
#         print "The supported OS's are: Linux, IRIX and IRIX64\n";
#         print "You can check the type of the system by: uname -s\n";
#         exit(0);
#     }
# }


if( $#ARGV < 1) {
    print "Usage: trans_mult.pl <results.txt> <res_num> [destination_dir]\n";
    exit(0);
}

open RES,$ARGV[0];

 
my @mols=();

my @selected_mols=();


my @match_list=();

my $filename;
my $referenceMol;
my $b_proc=0;
my $boolStartReadMatchList=0;
my @tmp;
my $tmp;
my $molid;
my $dest_dir="";

if($#ARGV == 2){
  if(! (-d $ARGV[2])){
    `mkdir $ARGV[2]`;
  }
  $dest_dir=$ARGV[2]."/";
}


my $counter=0;
    
while(<RES>){
    chomp;
    if(/^Mol-/){
        ($tmp, $filename)=split(' : ',$_);
        $mols[$#mols+1]=$filename;
        next;
    }
    if( /^Solution Num :/){
	if(/^Solution Num : $ARGV[1]/){
	    $b_proc=1;
	}else{
	    $b_proc=0;
	} 
    }

    if($b_proc==1 and /^Reference Molecule/){
	@tmp=split(':',$_);
        $filename=$mols[ ($tmp[1]+0)];

	my $dest_filename;
	if($#ARGV == 2){
	  @tmp=split('/',$filename);
	  $dest_filename=$tmp[$#tmp];
	}else{
	  $dest_filename=$filename;
	}
	
        
	my $tline=$dest_dir."0$counter"."_".$dest_filename;
        system("cp $filename $tline");
        print "Creating file $tline\n";
	
	$counter++;
	next;
    }


    if($b_proc==1 and /^Molecule/){
	@tmp=split(':',$_);
        $filename=$mols[ ($tmp[1]+0)];

      	next;
    }

    if($b_proc==1 and /^Trans :/){
     	my $tline=$_;
     
        @tmp=split(':',$tline);
	my $params=$tmp[1]; 
	
	
#$params=$tmp[0]." ".$tmp[1]." ".$tmp[2]." ".$tmp[3]." ".$tmp[4]." ".$tmp[5];
        #system("$home/utils/getPDBheader.pl $filename > $tmp_file2");


	my $dest_filename;
	if($#ARGV == 2){
	  @tmp=split('/',$filename);
	  $dest_filename=$tmp[$#tmp];
	}else{
	  $dest_filename=$filename;
	}
	
	if($counter<10){
	  $tline=$dest_dir."0$counter"."_".$dest_filename;
        }else{
	  $tline=$dest_dir."$counter"."_".$dest_filename;
	}
        system("cat $filename | $home/utils/pdb_trans_all_atoms$suffix $params > $tline");
        print "Creating file $tline\n";

	$counter++;
	next;
    }
}
