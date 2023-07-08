#!/usr/bin/perl

use strict;
use warnings;

#my $dir_to_open = "/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result_scheduler2";
my $dir_to_open = "./";
# Directory where the output files are

# my $dir_to_open_1 = "/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/scheduler/out";
my $dir_to_open_1 = "../scheduler/tmp/";
# Directory where you submit jobs

#my $dir_to_open_2 = "/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/scheduler/tmp";

#tmp/

opendir my $dh, $dir_to_open or die " Could not open $!\n";
# define a handle to open $dir_to_open

my $name   = $ARGV[0]; my $filenumber = $ARGV[1];
$name .= "_";
# $ARGV[0] is the first input argument: ./myPerl1_jobID.pl AB8F79ABBB628CC0B96256E4863BF059
my $JobID;

my @unsorted;
# creat an unsorted array 

# Read contents in the first directory
while(my $thing = readdir $dh)
{
    # Choose files which have the pattern "name"  
    if($thing =~ m/$name\w*/)
    {
	 my @aFiles   = split /_/, $thing;
	 #split the file name into several sections by "_"
	 #Example: AB8F79ABBB628CC0B96256E4863BF059_50.picoDst.result.root ==> @aFiles = {AB8F79ABBB628CC0B96256E4863BF059, 50.picoDst.result.root}
         my $vAppend  = $aFiles[1];
	 # 50.picoDst.result.root
               $JobID = $aFiles[0];
	 # AB8F79ABBB628CC0B96256E4863BF059
         my @aNumbers = split /\./, $vAppend;
	 #  50.picoDst.result.root ==> @aNumbers = {50, picoDst, result, root}
         push(@unsorted, $aNumbers[0]);
	 # Put 50 into array @unsorted

	 #print "$aNumbers[0]\n";
	 # print " $thing \n";
    }

}
closedir $dh;

my @sorted = sort{$a <=> $b} @unsorted;
# sort the array of numbers by order

my %missingnum;
# a hash of missing numbers
for(my $i=0; $i<$filenumber; $i++)
{
    $missingnum{$i} = 0;
}

foreach my $number (@sorted)
{
    $missingnum{$number} = 1;
    #print "$number\n";
}

chdir "${dir_to_open_1}";
# change directory to the one submitting jobs

my @missingDst;
# an array of missing Dst
foreach my $i (keys %missingnum)
{
    if($missingnum{$i}==0)
    {
      push(@missingDst, $i);

      # print ("star-submit -r $i $JobID.session.xml\n");

    }

}
# print  join(',',@missingDst);
my $jobString = join(',',@missingDst);
# make array of missingDst into a string separated by ","

if ($jobString ne "") {
    print ("star-submit -r $jobString $JobID.session.xml\n");
}
else {
    print ("All things were caught up! Good.\n");
}
# system("star-submit -r $jobString $JobID.session.xml\n");
