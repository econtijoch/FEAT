<h2>FMT Efficacy Analysis Toolkit</h2>
<p>A distributable package that includes a Shiny app and command-line functions to analyze the efficacy of FMT experiments.</p>
<p>The app will perform most all useful analysis downstream of generating an OTU table (for which you can use QIIME or any other method)</p>
<p>This app requires two packages that are not available on CRAN, but must be downloaded from bioconductor. To do this, simply type the following into the R command line:</p>

	source("https://bioconductor.org/biocLite.R")
	biocLite("rhdf5")
	biocLite("phyloseq")


<p>To launch the Shiny app, simply load the package and run the launchFEAT() command.</p>

<p>Alternatively, you can run it directly from GitHub:</p>


	require(shiny)
	runGitHub('econtijoch/FEAT/', subdir = 'inst/shiny-FEAT/')


