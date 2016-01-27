<h2>FMT Efficacy Analysis Toolkit</h2>
<p>A shiny app that provides an interactive UI to help analyze 16S data for FMT experiments.</p>
<p>The app will perform most all useful analysis downstream of generating an OTU table (for which you can use QIIME or NINJA-OPS)</p>

<p>The app has several dependencies listed below. It is also important to have the <a href=https://cran.r-project.org>latest version of R</a>. The app depends on the following packages:</p>
<ul>
	<li> shiny </li>
	<li> dplyr</li>
	<li> biom </li>
	<li> ggplot2 </li>
	<li> DT </li>
</ul>

<p> If you are starting from scratch, you can add all of these to your installation of R:</p>
<pre>
	install.packages(c('shiny', 'dplyr', 'biom', 'ggplot2', 'DT'))
</pre>

<p> When you have all of the dependencies installed, you can then run the app:</p>

<pre>
	require(shiny)
	runApp('FMT-Efficacy-Analysis-Toolkit')
</pre>

<p>It is important to note that the runApp() command's argument is a directory. For this app, it should be the name of the folder containing the contents of this GitHub repository.</p>