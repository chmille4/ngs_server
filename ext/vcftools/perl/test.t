#!/usr/bin/env perl
#
# Author: petr.danecek@sanger
#
# Usage: test.t [-d]
#

use strict;
use warnings;
use Carp;
use IPC::Open2;
use FindBin;
use lib "$FindBin::Bin";
use Vcf;

BEGIN {
    use Test::Most tests => 34;
}


my $path = $FindBin::RealBin;

my $debug = ($ARGV[0] && $ARGV[0] eq '-d') ? 1 : 0;

test_bgzip_and_tabix("$path/../examples/merge-test-a.vcf");
test_validator($path,"$path/../examples/valid-3.3.vcf");
test_validator($path,"$path/../examples/valid-4.0.vcf");
test_validator($path,"$path/../examples/valid-4.1.vcf");
test_validator($path,"$path/../examples/floats.vcf");
test_format_validation($path,'3.3');
test_format_validation($path,'4.0');
test_format_validation($path,'4.1');
test_vcf_stats($path,"$path/../examples/valid-4.0.vcf");
test_empty_cols($path,'4.0');
test_merge($path,'merge-test.vcf.out','merge-test-a.vcf','merge-test-b.vcf','merge-test-c.vcf');
test_compare($path,'cmp-test-a.vcf','cmp-test-b.vcf','cmp-test.out');
test_isec($path,'-n +2','isec-n2-test.vcf.out','merge-test-a.vcf','merge-test-b.vcf','merge-test-c.vcf');
test_query_vcf("$path/../examples/",'cmp-test-a.vcf','query-test.out','%CHROM:%POS\tref=%REF\talt=%ALT\tqual=%QUAL\t%INFO/DP[\t%SAMPLE=%GT]\n');
test_shuffle("$path/../examples/",'cmp-test-a.vcf','shuffle-test.vcf');
test_concat("$path/../examples/",'concat.out','concat-a.vcf','concat-b.vcf','concat-c.vcf');
test_annotate("$path/../examples/",'-c FROM,TO,CHROM,-,-,-,INFO/HM2,INFO/GN,INFO/DP -d key=INFO,ID=HM2,Number=0,Type=Flag,Description="HapMap2 membership" -d key=INFO,ID=GN,Number=1,Type=String,Description="Gene Name" -d key=INFO,ID=DP,Number=0,Type=Integer,Description="Depth,etc"','annotate.out','concat-a.vcf','annotate.txt');
test_annotate("$path/../examples/",'-c FROM,TO,CHROM,ID,REF,ALT,INFO/HM2,INFO/GN,INFO/DP -d key=INFO,ID=HM2,Number=0,Type=Flag,Description="HapMap2 membership" -d key=INFO,ID=GN,Number=1,Type=String,Description="Gene Name" -d key=INFO,ID=DP,Number=0,Type=Integer,Description="Depth,etc"','annotate3.out','concat-a.vcf','annotate.txt');
test_annotate("$path/../examples/",'-f +/D=34/c=2,3','annotate2.out','annotate-test.vcf');
test_fill_an_ac("$path/../examples/",'fill-an-ac.out','concat-a.vcf');
test_api_event_type([qw(A C),'s 1 C'],[qw(A ACGT),'i 3 CGT'],[qw(ACGT A),'i -3 CGT'],[qw(ACGT ACT),'i -1 G'],
    [qw(ACGT AAA),'o 3 AAA'],[qw(A .),'r 0 A'],[qw(A <ID>),'u 0 <ID>'],[qw(ACG AGC),'s 2 AGC'], [qw(A .A),'b'], [qw(A A.),'b']);

exit;

#--------------------------------------

sub test_bgzip_and_tabix
{
    my ($file) = @_;
    my $cmd;

    $cmd = "cat $file | bgzip -c > $file.gz";
    system($cmd);
    is(0,$?,"Is bgzip OK? .. $cmd");

    $cmd = "tabix $file.gz";
    system($cmd);
    is(0,$?,"Is tabix OK? .. $cmd");
}

sub test_validator
{
    my ($path,$fname) = @_;
    
    my $cmd = "perl -I$path -MVcf -e validate $fname";
    my @out = `$cmd 2>&1`;
    my @exp = ();
    is_deeply(\@out,\@exp,"Testing validator .. $cmd");
}

sub test_format_validation
{
    my ($path,$version) = @_;

    my ($chld_in,$chld_out);
    my $cmd = "perl -I$path -MVcf -e validate 2>&1";
    my $pid = open2($chld_out, $chld_in, $cmd);

    my $vcf = Vcf->new(version=>$version);
    $vcf->recalc_ac_an(2);
    $vcf->add_header_line({key=>'INFO', ID=>'AC',Number=>-1,Type=>'Integer',Description=>'Allele count in genotypes'});
    $vcf->add_header_line({key=>'INFO', ID=>'AN',Number=>1,Type=>'Integer',Description=>'Total number of alleles in called genotypes'});
    $vcf->add_header_line({key=>'FORMAT', ID=>'GT',Number=>1,Type=>'String',Description=>'Genotype'});
    if ( $version >= 4.0 )
    {
        $vcf->add_header_line({key=>'ALT',ID=>'DEL:ME:ALU', Description=>'Deletion of ALU element'});
    }
    if ( $version >= 4.1 )
    {
        $vcf->add_header_line({key=>'reference',value=>'file:/some/file.fa'});
        $vcf->add_header_line({key=>'contig',ID=>'1',length=>12345,md5=>'f126cdf8a6e0c7f379d618ff66beb2da',assembly=>'E.T.'});
    }
    $vcf->add_columns('NA0001','NA0002');
    print $vcf->format_header() unless !$debug;
    print $chld_in $vcf->format_header();

    my %rec = ( CHROM=>1, POS=>1, REF=>'A', QUAL=>$$vcf{defaults}{QUAL}, FORMAT=>['GT'] );
    $rec{gtypes}{NA0001}{GT} = 'A/A';
    $rec{gtypes}{NA0002}{GT} = $$vcf{defaults}{GT};
    $vcf->format_genotype_strings(\%rec);
    print $vcf->format_line(\%rec) unless !$debug;
    print $chld_in $vcf->format_line(\%rec);

    $rec{POS} = 2;
    $rec{gtypes}{NA0002}{GT} = 'IA|D1';
    if ( $version >= 4.0 ) 
    { 
        $rec{REF} = 'AC';
        $rec{gtypes}{NA0002}{GT} = 'ATC|<DEL:ME:ALU>'; 
    }
    $vcf->format_genotype_strings(\%rec);
    print $vcf->format_line(\%rec) unless !$debug;
    print $chld_in $vcf->format_line(\%rec);
    close($chld_in);

    my @exp = ();
    my @out = ();
    while (my $line=<$chld_out>)
    {
        chomp($line);
        push @out,$line;
    }
    close($chld_out);
    waitpid $pid, 0;

    if ( !is_deeply(\@out,\@exp,"Testing formatting followed by validation .. $cmd") )
    {
        print STDERR @out;
    }
}

sub test_vcf_stats
{
    my ($path,$file) = @_;
    my $cmd = "perl -I$path -MVcf $path/vcf-stats $file";
    my @out = `$cmd 2>&1`;
    open(my $fh,'<',"$file.stats") or confess("$file.stats: $!");
    my @exp = <$fh>;
    close($fh);

    is_deeply(\@out,\@exp,"Testing vcf-stats .. $cmd");
}

sub test_empty_cols
{
    my ($path,$version) = @_;

    my ($header,$vcf,@out,$exp);

    $vcf = Vcf->new(version=>$version);
    $vcf->add_header_line({key=>'FORMAT', ID=>'GT',Number=>1,Type=>'String',Description=>'Genotype'});
    $vcf->add_columns(qw(CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA0001));
    $header = $vcf->format_header();
    @out = split(/\n/,$header);
    $exp = join("\t",qw(CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA0001));
    is_deeply($out[-1],'#'.$exp,"Testing add_columns with genotypes full, $version.");

    $vcf = Vcf->new(version=>$version);
    $vcf->add_header_line({key=>'FORMAT', ID=>'GT',Number=>1,Type=>'String',Description=>'Genotype'});
    $vcf->add_columns('NA0001');
    $header = $vcf->format_header();
    @out = split(/\n/,$header);
    $exp = join("\t",qw(CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA0001));
    is_deeply($out[-1],'#'.$exp,"Testing add_columns with genotypes brief, $version.");

    $vcf = Vcf->new(version=>$version);
    $vcf->add_header_line({key=>'FORMAT', ID=>'GT',Number=>1,Type=>'String',Description=>'Genotype'});
    $vcf->add_columns();
    $header = $vcf->format_header();
    @out = split(/\n/,$header);
    $exp = join("\t",qw(CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO));
    is_deeply($out[-1],'#'.$exp,"Testing add_columns brief, $version.");

    $vcf = Vcf->new(version=>$version);
    $vcf->add_header_line({key=>'FORMAT', ID=>'GT',Number=>1,Type=>'String',Description=>'Genotype'});
    $vcf->add_columns('FORMAT');
    $header = $vcf->format_header();
    @out = split(/\n/,$header);
    $exp = join("\t",qw(CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO FORMAT));
    is_deeply($out[-1],'#'.$exp,"Testing add_columns no gtypes, $version.");
}

sub test_compare
{
    my ($path,$a,$b,$expected) = @_;

    my $curdir = `pwd`;
    chomp($curdir);
    chdir("$path/../examples");

    for my $file ($a,$b)
    {
        `cat $file | bgzip -c > $file.gz`;
        `tabix -p vcf -f $file.gz`;
    }
    
    my $cmd = "perl -I../perl/ -MVcf ../perl/vcf-compare -g $a.gz $b.gz";
    my @out = `$cmd 2>&1`;
    open(my $fh,'<',"$expected") or confess("$expected: $!");
    my @exp = <$fh>;
    close($fh);

    chdir($curdir);

    is_deeply(\@out,\@exp,"Testing vcf-compare .. $cmd");
}

sub test_merge
{
    my ($path,$expected,@files) = @_;

    my $curdir = `pwd`;
    chomp($curdir);
    chdir("$path/../examples");

    my $cmd = "perl -I../perl/ -MVcf ../perl/vcf-merge";
    for my $file (@files)
    {
        `cat $file | bgzip -c > $file.gz; tabix -f -p vcf $file.gz`;
        $cmd .= " $file.gz";
    }
    my @out = `$cmd 2>&1 | grep -v ^##source`;
    open(my $fh,'<',$expected) or confess("$expected: $!");
    my @exp = <$fh>;
    close($fh);

    chdir($curdir);
    is_deeply(\@out,\@exp,"Testing vcf-merge .. $cmd");
}

sub test_isec
{
    my ($path,$opts,$expected,@files) = @_;

    my $curdir = `pwd`;
    chomp($curdir);
    chdir("$path/../examples");

    my $cmd = "perl -I../perl/ -MVcf ../perl/vcf-isec -f $opts";
    for my $file (@files)
    {
        `cat $file | bgzip -c > $file.gz; tabix -f -p vcf $file.gz`;
        $cmd .= " $file.gz";
    }
    my @out = `$cmd 2>&1 | grep -v ^##source`;
    open(my $fh,'<',$expected) or confess("$expected: $!");
    my @exp = <$fh>;
    close($fh);

    chdir($curdir);
    is_deeply(\@out,\@exp,"Testing vcf-isec .. $cmd");
}


sub test_query_vcf
{
    my ($path,$file,$expected,$query) = @_;

    my $curdir = `pwd`;
    chomp($curdir);
    chdir("$path/../examples");

    my $cmd = "perl -I../perl/ -MVcf ../perl/vcf-query -f '$query' $file";
    my @out = `$cmd 2>&1`;
    open(my $fh,'<',$expected) or confess("$expected: $!");
    my @exp = <$fh>;
    close($fh);

    chdir($curdir);
    is_deeply(\@out,\@exp,"Testing vcf-query .. $cmd");
}


sub test_shuffle
{
    my ($path,$template,$file) = @_;

    my $curdir = `pwd`;
    chomp($curdir);
    chdir("$path/../examples");

    my $cmd = "perl -I../perl/ -MVcf ../perl/vcf-shuffle-cols -t $template $file";
    my @out = `$cmd 2>&1`;
    open(my $fh,'<',$template) or confess("$template: $!");
    my @exp = <$fh>;
    close($fh);

    chdir($curdir);
    is_deeply(\@out,\@exp,"Testing vcf-shuffle-cols .. $cmd");
}

sub test_concat
{
    my ($path,$out,@files) = @_;

    my $curdir = `pwd`;
    chomp($curdir);
    chdir("$path/../examples");

    my $cmd = "perl -I../perl/ -MVcf ../perl/vcf-concat -s 3";
    for my $file (@files)
    {
        `cat $file | bgzip -c > $file.gz`;
        `tabix -p vcf -f $file.gz`;
        $cmd .= " $file.gz";
    }
 
    my @out = `$cmd 2>&1`;
    open(my $fh,'<',$out) or confess("$out: $!");
    my @exp = <$fh>;
    close($fh);

    chdir($curdir);
    is_deeply(\@out,\@exp,"Testing vcf-concat .. $cmd");
}


sub test_annotate
{
    my ($path,$args,$out,$vcf,$annot) = @_;

    my $curdir = `pwd`;
    chomp($curdir);
    chdir("$path/../examples");

    my $cmd = "perl -I../perl/ -MVcf ../perl/vcf-annotate $args $vcf";

    if ( defined $annot )
    {
        `cat $annot | bgzip -c > $annot.gz`;
        `tabix -s 3 -b 1 -e 2 -f $annot.gz`;
        $cmd .= " -a $annot.gz";
    }

    my @out = `$cmd 2>&1 | grep -v ^##source`;
    open(my $fh,'<',$out) or confess("$out: $!");
    my @exp = <$fh>;
    close($fh);

    chdir($curdir);
    is_deeply(\@out,\@exp,"Testing vcf-annotate .. $cmd");
}

sub test_fill_an_ac
{
    my ($path,$out,$vcf) = @_;

    my $curdir = `pwd`;
    chomp($curdir);
    chdir("$path/../examples");

    my $cmd = "perl -I../perl/ -MVcf ../perl/fill-an-ac $vcf";
    my @out = `$cmd 2>&1`;
    open(my $fh,'<',$out) or confess("$out: $!");
    my @exp = <$fh>;
    close($fh);

    chdir($curdir);
    is_deeply(\@out,\@exp,"Testing fill-an-ac .. $cmd");
}

sub test_api_event_type
{
    my (@subs) = @_;
    my $vcf = Vcf->new();
    for my $mut (@subs)
    {
        my $exp = join(' ', $vcf->event_type($$mut[0],$$mut[1]));
        is_deeply($$mut[2],$exp,"Testing API event_type($$mut[0],$$mut[1]) .. $exp");
    }
}
