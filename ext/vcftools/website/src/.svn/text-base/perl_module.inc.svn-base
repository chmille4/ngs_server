
<div class="item">
<h1>The Perl modules and scripts</h1>
<p>VCFtools contains a Perl API (<a href="#Vcf.pm">Vcf.pm</a>) and a number
of Perl scripts that can be used to perform common tasks with VCF files such as
file validation, file merging, intersecting, complements, etc. The Perl tools
support all versions of the VCF specification (3.2, 3.3, 4.0 and 4.1), 
nevertheless, the users are encouraged to use the latest versions VCFv4.0 or VCFv4.1.
The VCFtools in general have been used mainly with diploid data, but the Perl tools
aim to support polyploid data as well.
</p>

<p>
Run any of the Perl scripts with the --help switch to obtain more help.
Note that the <b><span class="cmd">PERL5LIB</span> environment variable must contain 
the path to your VCFtools installation</b> in order for the scripts to work.
</p>
<p class="codebox">export PERL5LIB=/path/to/your/installation</p>

<p>
Many of the <b>Perl scripts require that the VCF files are compressed by <span class="cmd">bgzip</span>
and indexed by <span class="cmd">tabix</span></b> (both tools are part of the tabix package,
available for <a href="https://sourceforge.net/projects/samtools/files/tabix/">download here</a>).
The VCF files can be compressed and indexed using the following commands
<p class="codebox">bgzip my_file.vcf<br>
tabix -p vcf my_file.vcf.gz
</p>

<h3>The tools</h3>
<ul class="faq">
<li> <a href="#fill-aa">fill-aa</a>
<li> <a href="#fill-an-ac">fill-an-ac</a>
<li> <a href="#fill-ref-md5">fill-ref-md5</a>
<li> <a href="#fill-rsIDs">fill-rsIDs</a>
<li> <a href="#vcf-annotate">vcf-annotate</a>
<li> <a href="#vcf-compare">vcf-compare</a>
<li> <a href="#vcf-concat">vcf-concat</a>
<li> <a href="#vcf-convert">vcf-convert</a>
<li> <a href="#vcf-isec">vcf-isec</a>
<li> <a href="#vcf-merge">vcf-merge</a>
<li> <a href="#vcf-query">vcf-query</a>
<li> <a href="#vcf-shuffle-cols">vcf-shuffle-cols</a>
<li> <a href="#vcf-sort">vcf-sort</a>
<li> <a href="#vcf-stats">vcf-stats</a>
<li> <a href="#vcf-subset">vcf-subset</a>
<li> <a href="#vcf-to-tab">vcf-to-tab</a>
<li> <a href="#vcf-validator">vcf-validator</a>
<li> <a href="#Vcf.pm">Vcf.pm</a>
</ul>
</div>


<div class="item">
<a name="fill-aa" class="Q">fill-aa</a>
<p> Fill in ancestral alleles.
</p>
<p class="codebox">zcat file.vcf.gz | fill-aa -a ancestral-alleles.fa.gz | bgzip -c &gt; out.vcf.gz</p>
</div>


<div class="item">
<a name="fill-an-ac" class="Q">fill-an-ac</a>
<p> Fill or recalculate AN and AC INFO fields.
</p>
<p class="codebox">zcat file.vcf.gz | fill-an-ac | bgzip -c &gt; out.vcf.gz</p>
</div>


<div class="item">
<a name="fill-ref-md5" class="Q">fill-ref-md5</a>
<p>Fill missing reference info and sequence MD5s into VCF header.
</p>
<p class="codebox">fill-ref-md5 -i "SP:Homo\ Sapiens" -r ref.fasta in.vcf.gz -d ref.dict out.vcf.gz</p>
</div>


<div class="item">
<a name="fill-rsIDs" class="Q">fill-rsIDs</a>
<p>Fill missing rsIDs. This script has been discontinued, please use <a href="#vcf-annotate">vcf-annotate</a> instead.
</p>
</div>


<div class="item">
<a name="vcf-annotate" class="Q">vcf-annotate</a>
<p> Annotates VCF files with filtering information, IDs or custom annotations.
To add custom annotations to VCF files, create TAB delimited file with annotations
such as
<p class="codebox" style="white-space:pre;font-family:monospace;">#CHR     FROM   TO      ANNOTATION 
1        12345  22345   gene1 
1        67890  77890   gene2 
</p>
<p>
Compress the file (using <span class="cmd">bgzip annotations</span>), index (using <span class="cmd">tabix -s 1 -b 2 -e 3 annotations.gz</span>) and
run
</p>
<p class="codebox">cat in.vcf | vcf-annotate -a annotations.gz \ <br>
&nbsp;&nbsp; -d key=INFO,ID=ANN,Number=1,Type=Integer,Description='My custom annotation' \ <br>
&nbsp;&nbsp; -c CHROM,FROM,TO,INFO/ANN &gt; out.vcf </p>
<p style='padding-bottom:0px'>The script is also routinely used to apply filters. There are a number of predefined filters and 
custom filters can be easily added, see <span class="cmd">vcf-annotate -h</span> for examples. 
Some of the predefined filters take advantage of tags added by bcftools, the descriptions of the most frequently asked
ones follow:</p>
    <div style='margin-left:1em;margin-top:0.2em;'>
        <b>Strand Bias</b> 
        .. Tests if variant bases tend to come from one strand. Fisher's exact test
        for 2x2 contingency table where the row variable is being the reference allele
        or not and the column variable is strand. Two-tail P-value is used.
    </div>
    <div style='margin-left:1em;margin-top:0.2em;'>
        <b>End Distance Bias</b> 
        .. Tests if variant bases tend to occur
        towards the end of reads, which is usually an indication of
        misalignment. (T-test)
    </div>
    <div style='margin-left:1em;margin-top:0.2em;'>
        <b>Base Quality Bias</b> 
        .. Tests if variant bases tend to occur with a quality bias (T-test). This filter is by default 
        effectively disabled as it is set to 0.
    </div>
</div>

<div class="item">
<a name="vcf-compare" class="Q">vcf-compare</a>
<p> Compares positions in two or more VCF files and outputs the numbers of
positions contained in one but not the other files; two but not the other files, etc, which
comes handy when generating Venn diagrams. The script also computes numbers
such as nonreference discordance rates (including multiallelic sites), compares
actual sequence (useful when comparing indels), etc.
<br>
</p>
<p class="codebox">vcf-compare -H A.vcf.gz B.vcf.gz C.vcf.gz</p>
</div>




<div class="item">
<a name="vcf-concat" class="Q">vcf-concat</a>
<p>Very simple script for concatenating multiple VCF files into one. Checks that the count and
column names match.
</p>
<p class="codebox">vcf-concat A.vcf.gz B.vcf.gz C.vcf.gz | gzip -c &gt; out.vcf.gz</p>
</div>



<div class="item">
<a name="vcf-convert" class="Q">vcf-convert</a>
<p>Convert between VCF versions, currently from VCFv3.3 to VCFv4.0.
</p>
<p class="codebox">zcat file.vcf.gz | vcf-convert -r reference.fa &gt; out.vcf</p>
</div>


<div class="item">
<a name="vcf-isec" class="Q">vcf-isec</a>
<p> Creates intersections and complements of two or more VCF files. Given multiple
    VCF files, it can output the list of positions which are shared by at least N files,
    at most N files, exactly N files, etc. The first example below outputs positions shared by
    at least two files and the second outputs positions present in the files A but absent
    from files B and C.
</p>
<p class="codebox">
vcf-isec -n +2 A.vcf.gz B.vcf.gz | bgzip -c &gt; out.vcf.gz <br>
vcf-isec -c A.vcf.gz B.vcf.gz C.vcf.gz | bgzip -c &gt; out.vcf.gz
</p>
</div>


<div class="item">
<a name="vcf-merge" class="Q">vcf-merge</a>
<p> Merges two or more VCF files into one so that, for example, if two source files
    had one column each, on output will be printed a file with two columns.
</p>
<p class="codebox">vcf-merge A.vcf.gz B.vcf.gz C.vcf.gz | bgzip -c &gt; out.vcf.gz</p>
<p>Note that this script is not intended for concatenating VCF files. For this, use
<a href="perl_module.html#vcf-concat">vcf-concat</a> instead.
</p>
</div>


<div class="item">
<a name="vcf-query" class="Q">vcf-query</a>
<p> Powerful tool for converting VCF files into user-definable formats.
Supports retrieval of subsets of positions, columns and fields.
</p>
<p class="codebox">vcf-query file.vcf.gz 1:10327-10330</p>
</div>



<div class="item">
<a name="shuffle" class="Q">vcf-shuffle-cols</a>
<p>Reorder columns</p>
<p class="codebox"><a href="perl_module.html#vcf-shuffle-cols">vcf-shuffle-cols</a> -t template.vcf.gz file.vcf.gz > out.vcf
</p>
</div>


<div class="item">
<a name="vcf-sort" class="Q">vcf-sort</a>
<p>Sort a VCF file.
</p>
<p class="codebox">vcf-sort file.vcf.gz</p>
</div>


<div class="item">
<a name="vcf-stats" class="Q">vcf-stats</a>
<p>Outputs some basic statistics: the number of SNPs, indels, etc.
</p>
<p class="codebox">vcf-stats file.vcf.gz</p>
</div>


<div class="item">
<a name="vcf-subset" class="Q">vcf-subset</a>
<p>Remove some columns from the VCF file.
</p>
<p class="codebox">vcf-subset -c NA0001,NA0002 file.vcf.gz | bgzip -c &gt; out.vcf.gz</p>
</div>


<div class="item">
<a name="vcf-to-tab" class="Q">vcf-to-tab</a>
<p>A simple script which converts the VCF file into a tab-delimited text file listing the
    actual variants instead of ALT indexes.
</p>
<p class="codebox">zcat file.vcf.gz | vcf-to-tab &gt; out.tab</p>
</div>


<div class="item">
<a name="vcf-validator" class="Q">vcf-validator</a>
<p>
</p>
<p class="codebox">vcf-validator file.vcf.gz</p>
</div>


<div class="item">
<a name="Vcf.pm" class="Q">Vcf.pm</a>
<p>For examples how to use the Perl API, it is best to look at some of the simpler scripts, 
for example <span class="cmd">vcf-to-tab</span>. The detailed documentation can be obtained by running</p>
<p class="codebox">perldoc Vcf.pm</p>
</div>


