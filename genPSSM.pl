#!/gpfs/group/randen/home10/yoojin/programs/perl-5.8.9/perl

# Copyright (c) 2007, Gue Su Chang.
# All Rights Reserved.
#
# This module is free software; you can redistribute it
# and/or modify it under the same terms as Perl itself.
# This software is provided "as is" without warranty of any kind.
#
# Author:
#  Chang, Gue Su (gs(AT)medchem(DOT)info)
#     Graduate student
#     IBIOS Graduate Program Option in Bioinformatics and Genomics
#     Huck Institutes for Life Science
#     The Pennsylvania State University, UP
#
#     http://www.huck.psu.edu/IBIOS/

#  Modified by Yoojin Hong (yuh108@psu.edu) for AdaGB

# includes perl modules
use strict;
use warnings;
use Carp;

# add library path to @INC, for Unix and MSWin
use lib (
    'C:\Perl\shared',

    '/gpfs/group/randen/home10/yoojin/programs/perl-5.8.9/shared',
);

# for the necessary Perl modules
use FileHandle;
use DirHandle;
use File::Basename;
use File::Spec::Functions; 
use Date::Formatter;
use Cwd;

# BioPerl
use Bio::SeqIO;

# predeclares global variable names
use vars qw/$version $is_windows $scriptname $scriptpath/;

# determines if the operating system is Windows.
$is_windows = ($^O =~ /MSWin\d+/);


# TODO: adds the necessary modules
use FuncUtils;

use File::Copy;
use File::Path;
use File::Util;
use File::Basename;
use POSIX;
use Bio::SeqIO;
use GTool::Queries;
use GTool::Parallel;


## ------------ Setting global variables ------------------------------------- #

# predeclares global variable names
use vars qw/
            $seqpath $blastpgp $rpsblast $formatrpsdb $formatdb $copymat $makemat $proteindb $procno $localwater $evalue_h
            $startindex $dbname $dbsize $tempdir $iteration $dbdir $testpath
            $dbdir_cdd $dbdir_gt $dbdir_pssm $consensus $submatrix
            $seqs $dbdesc $scorescale $dbdir_check $threadno $maxtime
        /;

# the version
$version = '1.1.2';

# finds the perl script name and path
$scriptname = basename($0, '.pl') . '.pl';
$scriptpath = getcwd;       # dirname($0);



# input ------------------------------------------------------- #

# the path of the FASTA file containing a set of sequences
$seqpath = '';

# the path of the test FASTA file
# Note: the database checking process will be skipped if there is no test set
$testpath = '';

# the title and size of database 
$dbname = '';      # e.g. "virus"
$dbdesc = '';     # the database description


# the database size
# 5 means the maximum 10,000 entry number for the database
$dbsize = 5;

# the starting the index number (default: 1)
$startindex = 1;


# parameters ------------------------------------------------ #

# the e-value threshold for inclusion in multipass model (-h) (default: 0.01)
# Note: this threshold must be larger than or equal to 0.01;
$evalue_h = 1.0e-006;   # 1.0e-003, 1.0e-006, 1.0e-009

# the iteration number to construct each PSSM (default: 6)
$iteration = 2;

# the substitution matrix (default: BLOSUM62)
$submatrix = 'BLOSUM62';        # BLOSUM45


# the use of Smith-Waterman alignment algorithms (-s) (default: T)
$localwater = 'T';      # True or False

# the number of processors to use (-a) (default: 1)
# NOTE: Please don't change this option.
$procno = 1;

# the score scaling
$scorescale = 100;


# configuration ---------------------------------------------- #

# the protein sequence database for search in order to produce PSSMs
$proteindb = $is_windows ? 'D:\home\db\nr\nr' : '/gpfs/group/randen/home10/yoojin/DB/NR/nr';


# the path of the psi-BLAST executable file
$blastpgp = $is_windows ? 'D:\home\bin\blast\bin\blastpgp.exe' : '/gpfs/group/randen/home10/yoojin/programs/blast-2.2.18/bin/blastpgp';

# the path of the formatrpsdb executable file to produce the PSSM database
$formatrpsdb = $is_windows ? 'D:\home\bin\blast\bin\formatrpsdb.exe' : '/gpfs/group/randen/home10/yoojin/programs/blast-2.2.18/bin/formatrpsdb';

$makemat = $is_windows ? 'D:\home\bin\blast\bin\makemat.exe' : '/gpfs/group/randen/home10/yoojin/programs/blast-2.2.18/bin/makemat';
$copymat = $is_windows ? 'D:\home\bin\blast\bin\copymat.exe' : '/gpfs/group/randen/home10/yoojin/programs/blast-2.2.18/bin/copymat';
$formatdb = $is_windows ? 'D:\home\bin\blast\bin\formatdb.exe' : '/gpfs/group/randen/home10/yoojin/programs/blast-2.2.18/bin/formatdb';

# the path of the rpsblast executable file
$rpsblast = $is_windows ? 'D:\home\bin\blast\bin\rpsblast.exe' : '/gpfs/group/randen/home10/yoojin/programs/blast-2.2.18/bin/rpsblast';

# the number of child processes running
$threadno = $is_windows ? 4 : 12;

# the maximum running time in seconds (60 min = 3600)
$maxtime = 3600 * 2;


# the output directory for the PSSMs produced
# the three database will be stored together
#    "cdd": the PSSM database for the normal domain database
#    "gt": the PSSM database for the GTools
#    "pssm": the PSSMs as checkpoint and PSI-BLAST output
$dbdir = '';

# each directory name
$dbdir_cdd = 'cdd';
$dbdir_gt = 'gt';
$dbdir_pssm = 'pssm';       # used as a temporary directory together
$dbdir_check = 'test';

# the FASTA file name containing consensus sequences (or master sequences)
$consensus = '_masters.fa';



## ------------ Main subroutine ---------------------------------------------- #
# prints the about information
print About();

# main subroutine
Main();

# this program is terminated here
exit 1;



## ------------ Library of the main subroutine ----------------------------- #

sub Main {
 	  
    $seqpath = $ARGV[0];
    $testpath = $seqpath;
    $dbname = $ARGV[1];
    $dbdesc = $ARGV[2];
    $dbdir = $ARGV[3];
    #$dbdir = "./" . $dbdir_gt;
    $startindex = $ARGV[4];

    printf "seqpath = $seqpath\n";
    printf "testpath = $testpath\n";
    printf "dbname = $dbname\n";
    printf "dbdesc = $dbdesc\n";
    printf "dbdir = $dbdir\n"; 
    printf "startindex = $startindex\n";

    # creats the PSSM database
    # displays the parameters
    DisplayParameters();
    
    # sets up the path configuration
    SetupConfiguration();
    
    # indexes the sequences for PSSMs
    my $seqidx = IndexSequences();
    
    # creates the input files to construct PSSMs
    CreateInputfiles($seqidx);
    
    # constructs PSSMs of the sequences
    ConstructPSSMs($seqidx);

    # create profile ID list
    #CreateProfileIDList($seqidx);   
 
    # constructs the PSSM database
    ConstructDatabase($seqidx);
    
    # test the database files
    #CheckDatabase($seqidx);
    
    # program terminated
    printf "\nSuccessfully done.\n";
}



## ------------ Library of subroutines --------------------------------------- #


sub CheckDatabase {
    # test the database files
    my ($index) = @_;

    # local variables
    my $fu = File::Util->new();
    
    # the test sequence set
    if ($fu->existent($testpath))
    {
        printf "\nChecking databases: %s\n", $testpath;
    }
    else
    {
        printf "\nDB Checking Skipped: %s\n", $testpath;
        return;
    }
    
    # checks the individual PSSM
    # the parameters for formatrpsdb
    my $param = sprintf "-p T -e 0.01 -m 9 -F F -a %d", $procno;
    
    # creates the input file
    foreach my $id (sort keys %$index)
    {
        # the sequence index
        my $name = $index->{$id}->{name};
        #my $filename = $name . "_" . $id;
	my $filename = $name;
        my $seq = $index->{$id}->{seq};
        
        # the input file paths
        my $output = catfile $dbdir_check, $filename;
        my $dbpath = catfile $dbdir_gt, $name, $name;
        
        # the input FASTA file: $path.fa
        # the scoremat file: $filename.asnt
        
        # the command line for formatrpsdb
        my $cmd = sprintf "%s -i \"%s\" -d \"%s\" -o %s -l %s %s",
                                $rpsblast, $testpath, $dbpath, "$output.out", "$output.log", $param;

        printl $cmd;   # for debug
        system($cmd);
    }
    
    # the input file paths
    my $output = catfile $dbdir_check, $dbname;
    my $dbpath = catfile $dbdir_cdd, $dbname;
    
    # the command line for formatrpsdb
    my $cmd = sprintf "%s -i \"%s\" -d \"%s\" -o %s -l %s %s",
                            $rpsblast, $testpath, $dbpath, "$output.out", "$output.log", $param;
                            
    printl $cmd;   ## for debug
    system($cmd);
    
    printf "Checking databases finished: %s\n", $dbdir_check;
}


sub DisplayParameters {
    # displays the parameters

    # local variables

    # the e-value
    my $evalue = sprintf "%f", $evalue_h;

    print qq!
PSI-BLAST parameters:
\tE-value (for inclusion):\t$evalue
\tIterations:\t$iteration
\tSubstitution matrix:\t$submatrix
\tLocal Smith-Waterman:\t$localwater
 
FORMATDB parameters:
\tScale:\t$scorescale
!;
}


sub ConstructDatabase {
    # constructs the PSSM database
    my ($index) = @_;

    # local variables
    my $fu = File::Util->new();
    
    printf "\nConstructing databases:\n";
    
    # the logfile path
    my $log = catfile $dbdir, "formatdb.log";

    my $cmd = '';   

    chdir($dbdir_pssm);

    # creates the input file
    foreach my $id (sort keys %$index)
    {
        # the sequence index
        my $name = $index->{$id}->{name};
        #my $filename = $name . "_" . $id;
	my $filename = $name;
        my $seq = $index->{$id}->{seq};
        
        # the input file path
        my $path = catfile $dbdir_pssm, $filename;
        
        # the input FASTA file: $path.fa
        # the scoremat file: $filename.chk
      
	# runs the command for makemat
	$cmd = sprintf "%s -P %s",
				$makemat, $filename;
	printl $cmd;   # for debug
        system($cmd);

	# runs the command for copymat
	$cmd = sprintf "%s -P %s",
				$copymat, $filename;
	printl $cmd;   # for debug
        system($cmd);

	# runs the command for formatdb
	$cmd = sprintf "%s -i %s -l %s -o T -p T",
				$formatdb, $filename, $log;
	printl $cmd;   # for debug
        system($cmd);	

        # the PSSM database directory
        my $pssmdbdir = catfile $dbdir_gt, $name;
        mkpath($pssmdbdir);
        
        # copies the database files
        CopyDB($name, $pssmdbdir);

        printf "\t%s done and copied\n", $filename;
    }

    printf "\n%d databases constructed\n", $seqs->query_number;
    
    
    # constructs the total database
    my $path = catfile $dbdir_pssm, $dbname;
   
    # runs the command for makemat
    $cmd = sprintf "%s -P %s",
    			$makemat, $dbname;
    printl $cmd;   # for debug
    system($cmd);

    # runs the command for copymat
    $cmd = sprintf "%s -P %s",
                        $copymat, $dbname;
    printl $cmd;   # for debug
    system($cmd);

    # runs the command for formatdb
    $cmd = sprintf "%s -i %s -l %s -o T -p T",
                        $formatdb, $dbname, $log;
    printl $cmd;   # for debug
    system($cmd);
 
    
    # moves the database files
    CopyDB($dbname, $dbdir_cdd);
    
    printf "%s database of %d PSSMs constructed\n", $dbname, $seqs->query_number;
}


sub CopyDB {
    # copies the database files
    my ($name, $dir) = @_;

    # local variables
    
    # the list of the extensions of database files
    my @ext = ('aux', 'loo', 'phr', 'pin', 'psd', 'psi', 'psq', 'rps', 'mtx');

    # copies the database files
    foreach my $ext (@ext)
    {
        #printl "$name.$ext";        # only for debug
        move("$name.$ext", $dir);
    }
}


sub ConstructPSSMs {
    # makes the input files
    my ($index) = @_;

    # local variables
    my $fu = File::Util->new();
    
    printf "\nConstructing PSSMs:\n";

    # the PSI-BLAST parameter
    my $param = sprintf "-h %.1e -j %d -M %s -F F",
                            $evalue_h, $iteration, $submatrix;

    # creates the input file
    my @cmd;
    foreach my $id (sort keys %$index)
    {
        # the sequence index
        my $name = $index->{$id}->{name};
        #my $filename = $name . "_" . $id;
	my $filename = $name;
        my $seq = $index->{$id}->{seq};
        
        # the input file path
        my $path = catfile $dbdir_pssm, $filename;
        
        # the input FASTA file: $path.fa
        # the scoremat file: $filename.chk
        
        # the command line for PSI-BLAST
        my $cmd = sprintf "%s -d \"%s\" -i \"%s\" -o \"%s\" -C \"%s\" -e 0.01 -m 0 -s %s -a %d %s",
                                $blastpgp, $proteindb, "$path.fa", "$path.out", "$path.chk",
                                $localwater, $procno, $param;
        
        # runs the PSI-BLAST
        #printl $cmd;   ## for debug
        push @cmd, $cmd;
        
        #printf "\t%s done\n", $filename;
    }
    
    # the multiple process job
    my $pa = GTool::Parallel->new(
                        cmds => \@cmd,
                        
                        threadno => $threadno,
                        maxtime => $maxtime,
                        try => 3,
                        
                        logpath => $dbdir_pssm,
                    );
    
    # runs the given commands
    $pa->run;    

    printf "\n%d PSSMs constructed (parameter: %s)\n", $seqs->query_number, $param;
}


sub CreateInputfiles {
    # makes the input files
    my ($index) = @_;

    # local variables
    my $fu = File::Util->new();
    my $fu2 = File::Util->new();
 
    printf "\nCreating input files:\n";    
    
    # creates the consensus sequences
    my $cout = Bio::SeqIO->new(-file => "> $consensus", -format => 'fasta');
    
    # creates the input file
    my $list = '';
    my $list2 = '';
    foreach my $id (sort keys %$index)
    {
        # the sequence index
        my $name = $index->{$id}->{name};
        #my $filename = $name . "_" . $id;
	my $filename = $name;
        my $seq = $index->{$id}->{seq};
        
        # the input file path
        my $path = catfile $dbdir_pssm, $filename;
 
        # the input FASTA file
        my $out = Bio::SeqIO->new(-file => "> $path.fa", -format => 'fasta');
        $seq->display_id($name);
        $seq->desc("$name, $id, $dbdesc");
        $out->write_seq($seq);
        
        # for the consensus sequence file
        $cout->write_seq($seq);
        
        # creates the list file
        $fu->write_file(
           'file' => "$path.pn",
           'content' => "$path.chk",
           'bitmask' => 0644
        );
	$fu2->write_file(
           'file' => "$path.sn",
           'content' => "$path.fa",
           'bitmask' => 0644
        );

        
        # the total list file
        $list .= "$path.chk\n";
       	$list2 .= "$path.fa\n";
 
        printf "\t$path done\n"
    }

    # creates the total list file
    my $path = catfile $dbdir_pssm, $dbname;

    $fu->write_file(
       'file' => "$path.pn",
       'content' => $list,
       'bitmask' => 0644
    );
    $fu2->write_file(
       'file' => "$path.sn",
       'content' => $list2,
       'bitmask' => 0644
    );

   printf "\n\t$path.pn and $path.sn created\n";
}

sub CreateProfileIDList {
   my ($index) = @_;

   # local variables
   my $fu = File::Util->new();

   printf "\nCreating Profile ID List:\n"; 

   # creates the list file
   my $list = '';

   foreach my $id (sort keys %$index)
   {
	# the sequence index
        my $name = $index->{$id}->{name};

	# creates the list
	$list .= "$name\n";
   }
 
   # creates the total list file
    my $path = catfile $dbdir_pssm, $dbname;
 
    $fu->write_file(
       'file' => "$path.list",
       'content' => $list,
       'bitmask' => 0644
    );

}


sub IndexSequences {
    # indexes the sequences for PSSMs

    # local variables
    my %index;

    printf "\nMaking the PSSM index:\n";
    
    # the database title
    printf "Database title: $dbname\n";
    printf "Description: $dbdesc\n";
    printf "The index list:\n";
    
    # the sequence indice {seq=>, name=>}
    my $idx = $startindex;
    foreach my $id (sort $seqs->pids)
    {
        # the PSSM name like "cd00012"
        my $format = "%s%." . $dbsize . "d";
        my $name = sprintf $format, $dbname, $idx;
        
        #$index{$id} = {seq=>$seqs->queries->{$id}, name=>$name, idx=>$idx};
        $index{$id} = {seq=>$seqs->queries->{$id}, name=>$name};
        
        printf "\t$name =>\t$id\n";
        
        $idx ++;
    }
    
    # the result message
    printf "\n%d PSSM indices created\n", scalar(keys %index);
    
    return \%index;
}


sub SetupConfiguration {
    # sets up the path configuration

    # local variables
    my $fu = File::Util->new();
    
    printf "\nSetting the configuration:\n", $seqpath;

    # verifies the file paths
    croak "Cannot find " . $seqpath unless $fu->existent($seqpath);
    croak "Cannot find " . $blastpgp unless $fu->existent($blastpgp);
    croak "Cannot find " . $formatrpsdb unless $fu->existent($formatrpsdb);
    
    my ($name, $path) = fileparse($proteindb);
    croak "Cannot find " . $path unless $fu->existent($path);
    croak "Cannot find " . $proteindb unless $fu->existent($proteindb . ".pin");

    # opens the fasta file containing sequences
    $seqs = GTool::Queries->new(
                        path => $seqpath,
                        name => $dbname,
                    );
    
    if ($seqs->open())
    {
        printf "Reading %d sequences...\n", $seqs->query_number;
        printf "  from $seqpath\n"
    }
    else
    {
        croak $seqs->get_error
    }
    
    # checks the database size
    #printl POSIX::pow(10, $dbsize); # for debug
    croak "Please increase the database size: $dbsize" if $seqs->query_number >= POSIX::pow(10, $dbsize);

    # sets the working directories
    $dbdir_cdd = catfile $dbdir, $dbdir_cdd;
    $dbdir_gt = catfile $dbdir, $dbdir_gt;
    $dbdir_pssm = catfile $dbdir, $dbdir_pssm;
    $dbdir_check = catfile $dbdir, $dbdir_check;
    $consensus = catfile $dbdir, $dbname . $consensus;
 
    # creates the working direcotories
    printf "Creating working directories...\n";
    mkpath([$dbdir_cdd, $dbdir_gt, $dbdir_pssm, $dbdir_check], 1);

    # checks the output directories
    croak "Cannot find " . $dbdir unless $fu->existent($dbdir);
    croak "Cannot find " . $dbdir_cdd unless $fu->existent($dbdir_cdd);
    croak "Cannot find " . $dbdir_gt unless $fu->existent($dbdir_gt);
    croak "Cannot find " . $dbdir_pssm unless $fu->existent($dbdir_pssm);
    
    printf "\nThe job successfully set up in $dbdir\n";
    
    return;
}


## ------------ Library of basic subroutines --------------------------------- #

sub About {
    # param: none, return: string
    # TODO: adds the description of this program
    return qq!
$scriptname $version   (Created by Gue Su Chang)
  The DB manager for PSSM (Position-specific scoring matrix)
!;
}


