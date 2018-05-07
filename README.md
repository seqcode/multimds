# MultiMDS

MultiMDS is a tool for locus-specific structural comparisons of two Hi-C datasets. It jointly infers and aligns 3D structures from two datasets, such as different cell types. The output is aligned 3D structure files (which can be plotted, see below) and a list of loci that significantly relocalize between the datasets (measured using Euclidean distance). These may represent A/B compartment changes, changes in enhancer localization, or other structural changes. 

## Installation

Requirements:
* python 2.7
* Python dependencies can be installed using
``pip install -r requirements.txt``
* The following optional dependencies can be installed manually:
    * [mayavi](http://docs.enthought.com/mayavi/mayavi/) (for plotting)
    * [ImageMagick](https://www.imagemagick.org/script/index.php) (for creating gifs)

## TLDR

``python relocalization_peaks.py [Hi-C BED path 1] [Hi-C BED path 2]``

## Usage

### Input files

MultiMDS uses intrachromosomal BED files as input. Data must be normalized prior to use (for example, using [HiC-Pro](http://nservant.github.io/HiC-Pro/)). 

Format:

>chrom	bin1\_start	bin1\_end	chrom	bin2\_start	bin2\_end	normalized\_contact\_frequency

Example - chr22 data at 10-Kbp resolution:

>chr22	16050000	16060000	chr22	16050000	16060000	12441.5189291
> 
>...

Important: BED file 1 and BED file 2 must be the same species, chromosome, and resolution!

### Running the program

To view help:

``python relocalization_peaks.py -h``

To run with default parameters:

``python relocalization_peaks.py [BED path 1] [BED path 2]``

For example:

``python relocalization_peaks.py GM12878_combined_22_100kb.bed K562_22_100kb.bed``

#### Parameters (optional)

##### Number of partitions

Partitioning is used in the structural inference step for greater efficiency and accuracy. By default 4 partitions are used. This can be controlled with the -N parameter: 

``python relocalization_peaks.py -N 2 GM12878_combined_22_100kb.bed K562_22_100kb.bed``

##### Centromere

Better results are achieved if the chromosome is partitioned at the centromere in the partitioning step. The genomic coordinate of the centromere can be entered with the -m parameter

``python relocalization_peaks.py -m 28000000 GM12878_combined_20_100kb.bed K562_220_100kb.bed``

### Output files

### Relocalization peaks

Genomic regions that significantly relocalize between the cell types are saved to a BED file, with the format [PREFIX1]_[PREFIX2]_peaks.bed

For example the output of

``python relocalization_peaks.py GM12878_combined_22_100kb.bed K562_22_100kb.bed``

is GM12878_combined_22_100kb_K562_22_100kb_peaks.bed

#### Structure files
Aligned structures are saved to tsv files, which can be used for plotting (see below). The header contains the name of the chromosome, the resolution, and the starting genomic coordinate. Each line in the file contains the genomic bin number followed by the 3D coordinates (with "nan" for missing data). 

Example - chr22 at 10-Kbp resolution:

>chr22
> 
>10000
> 
>16050000
> 
>0	0.589878298045	0.200029092422	0.182515056542
> 
>1	0.592088232028	0.213915817254	0.186657230841
> 
>2	nan	nan	nan
> 
>...

0 corresponds to the bin 16050000-16060000, 1 corresponds to the bin 16060000-16070000, etc. 

### Plotting

Read a structure:

    import data_tools
    structure = data_tools.structure_from_file("GM12878_combined_22_100kb_structure.tsv")``

Create an interactive 3D plot in Mayavi. (Mayavi allows you to rotate the image and save a view.)

    import plotting
    plotting.plot_structure_interactive(structure, color=(0,0.5,0.7), radius=0.01, enrichments=my_enrichments)``

If _radius_ is not selected, the to-scale radius of heterochromatin is used. 

_enrichments_ is a vector with a numerical value for each bin in the structure (i.e. bins that do not have a nan coordinate). For example, this could represent ChIP-seq enrichments for each bin. This option overrides _color_ and will use a rainbow colormap, with blue representing low values and red representing high values. 

Multiple structures can be plotted simultaneously:

    chroms = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, X)
    structures = [data_tools.structure_from_file("GM12878_combined_{}_100kb_structure.tsv".format(chrom) for chrom in chroms)]
    plotting.plot_structures_interactive(structures)

plotting.py has 23 built-in colors designed to be maximally different to the human eye. By default, these colors are used when plotting multiple structures. You can also specify a list of colors:

    chroms = (1, 2)
    structures = [data_tools.structure_from_file("GM12878_combined_{}_100kb_structure.tsv".format(chrom) for chrom in chroms)]
    plotting.plot_structures_interactive(structures, colors=[(1,0,0), (0,0,1)])

_all_enrichments_ is a list of enrichments, e.g. 
     
     plotting.plot_structures_interactive(structures, all_enrichments=[enrichments1, enrichments2])

The radius can also be specified, as above. 

The option _cut_ creates a cross-section of the plot. For example, this is useful for viewing the interior of the nucleus.

    chroms = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, X)
    structures = [data_tools.structure_from_file("GM12878_combined_{}_100kb_structure.tsv".format(chrom) for chrom in chroms)]
    plotting.plot_structures_interactive(structures, cut=True)

A plot can be saved as a gif:

``plotting.plot_structure_gif(structure, struct, color=(1,0,0), radius=None, increment=10)``

will create struct.gif

A smaller value of _increment_ will lead to a smoother gif. Increments must be a factor of 360. 

Multiple structures can also be plotted in a single gif:

``plotting.plot_structures_gif(structures, struct, colors=default_colors, radius=None, increment=10)``
