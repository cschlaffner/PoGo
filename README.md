# PoGo

## Overview
In proteogenomic analyses it is essential to know the loci giving rise to peptides in order to improve genomic annotation and the functional characterization of protein products in their biological context. With next-generation sequencing of DNA and RNA for each sample studied by proteomic mass spectrometry integration and visualisation in a common coordinate system, i.e. the genome, is vital for systems biology. Advances in technology in mass spectrometry now allow almost complete quantification of the sample proteome. With research moving to protein quantitative trait loci (pQTL) to identify genomic alterations with functional effects on the proteome and the high complexity of combinations thereof integration and visualisation of protein and peptide quantification on genomic loci is paramount for this type of analysis. Furthermore, moving towards more personal multi-omics studies comparative visualisation of proteomic data on a genome has been lacking. Not only genomic variation affecting proteins have come into focus of functional integration studies but also post-translational modifications (PTM), the effect of single nucleotide variants and other alterations on PTMs and alternative modification loci, and the effects of alternative PTMs on protein abundance have become more a centre of attention for researchers. To facilitate this type of integration not only the genomic locations of modified peptides but specifically the genomic loci of associated with these modifications is required. Here, we provide a mapping tool, PoGo, to quickly and efficiently identify genomic loci of peptides and post-translational modifications and couple these mappings with associated quantitative values over multiple samples. Using reference gene annotation and an associated transcript translations our tool identifies the genomic loci of peptides given as input and generates output in different formats borrowed from genomics and transcriptomics which can be loaded in various genome browsers such as UCSC Genome Browser, Ensembl Genome Browser, BioDalliance, and the Integrative Genomics Viewer.

## Download and Instiallation
PoGo source code is written to support compilation on Windows, Mac and Linux systems. Executables can be downloaded here: ftp://ftp.sanger.ac.uk/pub/teams/17/software/PoGo/. 

## Learn and Support
PoGo uses transcript translations and reference gene annotations to identify the genomic loci of peptides and post-translational modifications. Multiple occurrences of peptides in the input data resulting in the same genomic loci will be collapsed as a single occurrence in the output.

### Input format
The input format required by PoGo is a tab delimited file with four columns.

<table border="0" width="100%"><thead><tr><th scope="col">Column</th><th scope="col">Column header</th><th scope="col">Description</th></tr></thead><tbody><tr><td>1</td><td>Sample</td><td>Name of sample or experiment</td></tr><tr><td>2</td><td>Peptide</td><td>Peptide sequence with PSI-MS nodification names in round brackets following the mpdified amino acid, e.g. PEPT(Phopsho)IDE for a phosphorylated threonine</td></tr><tr><td>3</td><td>PSMs</td><td>Number of peptide-spectrum matches (PSMs) for the given peptide</td></tr><tr><td>4</td><td>Quant</td><td>Quantitative value for the given peptide in the given sample</td></tr></tbody></table>

### Output formats

#### BED
This format contains the genomic loci for peptides, the exon-structure, the peptide sequence, as well as a colour code for uniqueness of peptides within the genome.

<table align="left" border="0" width="100%">
	<thead>
		<tr>
			<th scope="col" width="20%">Colour</th>
			<th scope="col" width="80%">Description</th>
		</tr>
	</thead>
	<tbody>
		<tr>
			<td bgcolor="#F00000"><img src="http://ngs.sanger.ac.uk/production/proteogenomics/trackhubs_files/color/uniquetranscript.svg" height="25px"></td>
			<td>Peptide is unique to single gene AND single transcript</td>
		</tr>
		<tr></tr>
		<tr>
			<td bgcolor="#000000"><img src="http://ngs.sanger.ac.uk/production/proteogenomics/trackhubs_files/color/uniquegene.svg" height="25px"></td>
			<td>Peptide is unique to single gene BUT shared between multiple transcripts</td>
		</tr>
		<tr></tr>
		<tr>
			<td bgcolor="#808080"><img src="http://ngs.sanger.ac.uk/production/proteogenomics/trackhubs_files/color/notunique.svg" height="25px"></td>
			<td>Peptide is shared between multiple genes</td>
		</tr>
	</tbody>
</table>

#### PTMBED
Like BED but containing the location of the post-translational modification on the genome. Thick parts of the peptide blocks indicate the position of the post-translational modification on a single amino acid (short thick block) while longer blocks indicate the occurrence of the first and last post-translational modification and residues in between. In the PTMBED the colour code is changed to indicate the type of modification.

<table border="0" width="100%">
	<thead>
		<tr>
			<th scope="col" width="20%">Colour</th>
			<th scope="col" width="80%">Post-translational Modification</th>
		</tr>
	</thead>
	<tbody>
		<tr>
			<td bgcolor="#FF3333"><img src="http://ngs.sanger.ac.uk/production/proteogenomics/trackhubs_files/color/phospho.svg" height="25px"></td>
			<td>Phosphorylation (phospho)</td>
		</tr>
		<tr></tr>
		<tr>
			<td bgcolor="#CC6600"><img src="http://ngs.sanger.ac.uk/production/proteogenomics/trackhubs_files/color/acetyl.svg" height="25px"></td>
			<td>Acetylation (acetyl)</td>
		</tr>
		<tr></tr>
		<tr>
			<td bgcolor="#FF9933"><img src="http://ngs.sanger.ac.uk/production/proteogenomics/trackhubs_files/color/amidated.svg" height="25px"></td>
			<td>Amidation (amidated)</td>
		</tr>
		<tr></tr>
		<tr>
			<td bgcolor="#CCCC00"><img src="http://ngs.sanger.ac.uk/production/proteogenomics/trackhubs_files/color/oxidation.svg" height="25px"></td>
			<td>Oxidation (oxidation)</td>
		</tr>
		<tr></tr>
		<tr>
			<td bgcolor="#00CC00"><img src="http://ngs.sanger.ac.uk/production/proteogenomics/trackhubs_files/color/methyl.svg" height="25px"></td>
			<td>Methylation (methyl)</td>
		</tr>
		<tr></tr>
		<tr>
			<td bgcolor="#33FF33"><img src="http://ngs.sanger.ac.uk/production/proteogenomics/trackhubs_files/color/glygly.svg" height="25px"></td>
			<td>Ubiquitinylation (glygly; gg)</td>
		</tr>
		<tr></tr>
		<tr>
			<td bgcolor="#00CCCC"><img src="http://ngs.sanger.ac.uk/production/proteogenomics/trackhubs_files/color/sulfo.svg" height="25px"></td>
			<td>Sulfation (sulfo)</td>
		</tr>
		<tr></tr>
		<tr>
			<td bgcolor="#3399FF"><img src="http://ngs.sanger.ac.uk/production/proteogenomics/trackhubs_files/color/palmitoyl.svg" height="25px"></td>
			<td>Palmitoylation (palmitoyl)</td>
		</tr>
		<tr></tr>
		<tr>
			<td bgcolor="#0000CC"><img src="http://ngs.sanger.ac.uk/production/proteogenomics/trackhubs_files/color/formyl.svg" height="25px"></td>
			<td>Formylation (formyl)</td>
		</tr>
		<tr></tr>
		<tr>
			<td bgcolor="#3333FF"><img src="http://ngs.sanger.ac.uk/production/proteogenomics/trackhubs_files/color/deamidated.svg" height="25px"></td>
			<td>Deamidation (deamidated)</td>
		</tr>
		<tr></tr>
		<tr>
			<td bgcolor="#FF3399"><img src="http://ngs.sanger.ac.uk/production/proteogenomics/trackhubs_files/color/any.svg" height="25px"></td>
			<td>Any other post-translational modification</td>
		</tr>
	</tbody>
</table>

#### GTF
This output format contains besides the genomic loci the annotated information for the genes giving rise to each peptide sequence including status and biotype. For each mapped peptide the sample, number of peptide-spectrum matches and associated quantitative value as tags.

#### GCT
In this format the peptide sequences are combines with the Ensembl gene identifier. It contains the genomic loci for each peptide as well as the quantitative values for each peptide in different samples as a matrix.

### Usage
PoGo is a readily compiled command line tool written in C++, basic usage to generate map peptides to corresponding genomic loci with default parameters:

`Windows:	PoGo.exe -fasta annotation.translation.fasta -gtf annotation.gtf -in example.tsv`
`Linux/Unix/Mac:	PoGo -fasta annotation.translation.fasta -gtf annotation.gtf -in example.tsv`

Full usage:
`PoGo/PoGo.exe -fasta TRANSL -gtf ANNO -in *.tsv[,*.tsv] [-format OUTF] [-merge TRUE/FALSE] [-source SRC] [-mm NUM] [-mmmode TRUE/FALSE]`

Required arguments:
<table border="0" widht="100%"><tbody><tr><td width="20%">
<pre>-fasta TRANSL</pre>
</td><td>Filepath for file containing protein sequences in FASTA format</td></tr><tr><td>
<pre>-gtf ANNO</pre>
</td><td width="80%">Gene annotation with coding sequences (CDS) in GTF format</td></tr><tr><td>
<pre>-in *.tsv</pre>
</td><td>Path to single input file or comma separated list of paths to input files containing peptides to be mapped with associated number of peptide to spectrum matches, sample name and quantitative value (see input file format)</td></tr></tbody></table>

Optional arguments:
<table border="0" width="100%"><tbody><tr><td width="20%">
<pre>-format OUTF</pre>
</td><td width="80%">Set output format GTF, GCT, BED, PTMBED or ALL. Comma separated combination possible. Default = ALL</td></tr><tr><td>
<pre>-merge TRUE/FALSE</pre>
</td><td>Set TRUE to merge output of multiple input files (output will be named after last input file *_merged). Default = FALSE</td></tr><tr><td><pre>-source SRC</pre></td><td>Set TRUE to merge output of multiple input files (output will be named after last input file *_merged). Default = FALSE</td></tr><tr><td>
<pre>-mm NUM</pre>
</td><td>Number of mismatches allowed in mapping (0, 1 or 2). DEFAULT = 0</td></tr><tr><td>
<pre>-mmmode TRUE/FALSE</pre>
</td><td>Set TRUE to restrict number of mismatch in kmer to 1. DEFAULT = FALSE</td></tr><tr><td>
<pre>-species SPECIES</pre>
</td><td>Set species using common or scientific name or taxonomy ID. Default is Human (Homo sapiens, 9606).</td></tr></tbody></table>

Table of supported species
<table border="0" width="100%"><thead>
<tr><th>Common name</th><th>Scientific name</th><th>Taxon ID</th></tr></thead><tbody>
<tr><td>C.intestinalis</td><td>Ciona intestinalis</td><td>7719</td></tr>
<tr><td>Cat</td><td>Felis catus</td><td>9685</td></tr>
<tr><td>Chicken</td><td>Gallus gallus</td><td>9031</td></tr>
<tr><td>Chimpanzee</td><td>Pan troglodytes</td><td>9598</td></tr>
<tr><td>Cow</td><td>Bos taurus</td><td>9913</td></tr>
<tr><td>Dog</td><td>Canis lupus familiaris</td><td>9615</td></tr>
<tr><td>Gorilla</td><td>Gorilla gorilla gorilla</td><td>9595</td></tr>
<tr><td>Horse</td><td>Equus caballus</td><td>9796</td></tr>
<tr><td>Human</td><td>Homo sapiens</td><td>9606</td></tr>
<tr><td>Macaque</td><td>Macaca mulatta</td><td>9544</td></tr>
<tr><td>Marmoset</td><td>Callithrix jacchus</td><td>9483</td></tr>
<tr><td>Medaka</td><td>Oryzias latipes</td><td>8090</td></tr>
<tr><td>Mouse</td><td>Mus musculus</td><td>10090</td></tr>
<tr><td>Olive baboon</td><td>Papio anubis</td><td>9555</td></tr>
<tr><td>Opossum</td><td>Monodelphis domestica</td><td>13616</td></tr>
<tr><td>Orangutan</td><td>Pongo abelii</td><td>9601</td></tr>
<tr><td>Pig</td><td>Sus scrofa</td><td>9823</td></tr>
<tr><td>Platypus</td><td>Ornithorhynchus anatinus</td><td>9258</td></tr>
<tr><td>Rabbit</td><td>Oryctolagus cuniculus</td><td>9986</td></tr>
<tr><td>Rat</td><td>Rattus norvegicus</td><td>10116</td></tr>
<tr><td>Sheep</td><td>Ovis aries</td><td>9940</td></tr>
<tr><td>Tetraodon</td><td>Tetraodon nigroviridis</td><td>99883</td></tr>
<tr><td>Turkey</td><td>Meleagris gallopavo</td><td>9103</td></tr>
<tr><td>Vervet-AGM</td><td>Chlorocebus sabaeus</td><td>60711</td></tr>
<tr><td>Zebra Finch</td><td>Taeniopygia guttata</td><td>59729</td></tr></tbody>
<tr><td>Zebrafish</td><td>Danio rerio</td><td>7955</td></tr></tbody>
<tr><td>Yeast</td><td>Saccheromyces cerevisiae</td><td>4932</td></tr></tbody>
</table>

### Step by step
<ol><li>Download annotation and translated sequences for human from GENCODE, e.g. release 25. Go to www.gencodegenes.org/release/25.html and download the GTF file containing 'Comprehensive gene annotation' and the 'Protein-coding transcript translation sequences' as Fasta file. Store and unzip both files into a folder, e.g. ${POGO_DIR}/input/</li><li>Navigate to the folder that contains the PoGo executable (cd ${POGO})</li><li>Execute the following command to generate gtf, gct, bed and ptmbed output for of the your input file referred to as ${Peptides.txt}<br><newline>Linux/Unix</newline>

<pre>./PoGo -fasta ./input/gencode.v25.pc_translations.fa -gtf ./input/gencode.v25.annotation.gtf -in /PATH/TO/${Peptides.txt}</pre>
<newline>Windows</newline>

<pre>.\PoGo.exe -fasta .\input\gencode.v25.pc_translations.fa -gtf .\input\gencode.v25.annotation.gtf -in \PATH\TO\${Peptides.txt}</pre>
</li><li>You can load the generated BED and/or GTF files into a genome browser or create an web accessible track hub for your data through the TrackHub Generator (http://www.sanger.ac.uk/science/tools/trackhub-generator). Here the example is shown for visualisation in the UCSC genome browser.<br><newline>To load the data into the browser please follow these steps:</newline>
<ol type='a'><li>Go to https://genome.ucsc.edu and navigate to 'My Data' -&gt; 'Custom Tracks'.</li><li>After clicking 'Choose File' select the file you want to upload and submit via the 'Submit' button.</li><li>You will be redirected to the 'Manage Custom Tracks' webpage.</li><li>Proceed from the 'Manage Custom Tracks' page by selecting 'Genome Browser' and confirm ('go').</li><li>Now you can browse the peptides mapped to their genomic loci on the reference genome.</li></ol></li></ol>

### Runtime and Memory Estimation
Runtime and required memory (RAM) for PoGo execution across different settings for inclusion of mismatches and depending on number of peptides in the input file.
<table border="0" width="100%"><thead><tr><th scope="col">Mismatch Parameter Setting</th><th scope="col">~1,500 peptides in Input</th><th scope="col">~250,000 peptides in Input</th></tr></thead><tbody><tr><td>
<pre>-mm 0</pre>
</td><td>&lt; 5 min / &lt; 4 GB</td><td>&lt;5 min / min 10 GB</td></tr><tr><td>
<pre>-mm 1</pre>
</td><td>~5 min / &lt; 4 GB</td><td>~5 min / min 16 GB</td></tr><tr><td>
<pre>-mm 2 -mmode true</pre>
</td><td>&lt; 15 min / min 6 GB</td><td>&lt; 20 min / min 32 GB</td></tr><tr><td>
<pre>-mm 2 -mmode false</pre>
</td><td>~ 1.5 h / min 64 GB</td><td>~ 2 h / min 160 GB</td></tr></tbody></table>

### Test Examples
Test examples, requirement specifications and time estimations are available here: ftp://ftp.sanger.ac.uk/pub/teams/17/software/PoGo/PoGo_Testprocedures.zip.

## Contact
Christoph Schlaffner (christoph.schlaffner@childrens.harvard.edu)
