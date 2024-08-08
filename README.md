# Molecular Mimicry as a Mechanism of Viral Immune Evasion and Autoimmunity
<a href="https://doi.org/10.5281/zenodo.13272863"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.13272863.svg" alt="DOI"></a>

This is the repository for the paper "Molecular Mimicry as a Mechanism of Viral Immune Evasion and Autoimmunity" currently pre-printed at BioRxiv. Below is a brief summary of the repository structure.

All analysis outside of the alignment including Figure assembly and statistical testing is in the "Analysis" folder. All k-mer alignment is in the epitope_alignment folder. Within both of these directories, there were several instances of files being too large to be housed on github. In these instances we have included a link to where the files are housed on Box if you desired the exact file used in this manuscript. The only instances were we did not do this are for some files in the "Databases" folder under "Analysis" which were from other sources/manuscripts. There we have reported where to retrieve the file from.

Many of the figure's quarto documents have intermediate data files uploaded on the github (file size permitting, otherwise a link was provided to where you can retrieve the large file) which will enable you to run the figure code directly out of the box where permitting. The alternative is to go into the respective quarto document and set eval: true for many of the longer running chunks that generate these intermediate files (e.g. the chunks that save data to Key_Data). So also pay attention to which chunks have eval set to false, as this was usually done because the file they generate is provided and running that chunk everytime greatly slows down the speed of rendering the quarto.

On 08-08-2024 changes made for revisions to the manuscript in response to reviewer's comments were finalized (v1.1.1). 

If you use any of the code or novel data (so excluding data in the Databases folder) from this github, please cite the manuscript at <a href="https://zenodo.org/doi/10.5281/zenodo.11411891">https://doi.org/10.1101/2024.03.08.583134</a>.
