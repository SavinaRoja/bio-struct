#!/usr/bin/perl

#ver 1.2

if($#ARGV<1){
    print("Usage:joinXpdb.pl pdbFile1 pdbFile2 ..\n");
    exit 0;
}

my $index=1;
foreach $file (@ARGV){
  print "MODEL     ";
  printf "%4d   $file\n",$index;
  $index++;

  open FILE,"$file";
  #my @array = grep !/MODEL/, <FILE>;
  while(<FILE>){
    if(/^MODEL/){
      print "MODEL     ";
      printf "%4d   $file",$index;
      my (undef, $model)=split /MODEL/;
      print " $model\n";

      $index++;
    }else{
      print $_;
    }
  }

  print @array;
  close FILE;
  print "ENDMDL \n"; 
}

