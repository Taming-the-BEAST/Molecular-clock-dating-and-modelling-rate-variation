---
author: Simon Ho,Jacqueline Nguyen
level: Intermediate
title: Molecular clock dating and modelling rate variation
subtitle: Accounting for evolutionary rate variation in passerine birds
beastversion: 2.7.7
---




# Background

Evolutionary rate variation is an intrinsic feature of biological data, being driven by a range of life-history, environmental, and biochemical factors. Rates of molecular evolution can vary across lineages, across nucleotide sites, and across regions of the genome. These forms of variation can interact to produce complex patterns of heterogeneity in the data, particularly at the genomic scale {% cite Ho2014 --file Molecular-clock-dating-and-modelling-rate-variation/master-refs.bib %}. Thus, accounting for evolutionary rate variation represents an important aspect of phylogenetic analysis and molecular dating {% cite dosReis2016 --file Molecular-clock-dating-and-modelling-rate-variation/master-refs.bib %}. A range of models are available for this purpose. 

This tutorial will provide an introduction to using multiple relaxed clock models to account for rate variation across lineages and across subsets of the data set. The tutorial will use mitochondrial genome data from passerine birds. The exercise will guide you through setting up two contrasting models of among-lineage rate variation to estimate phylogenetic relationships, evolutionary rates, and node times using the program BEAST v{{ page.beastversion }}.


----

# Programs used in this Exercise 



----

# Practical: Using multiple clock models

This tutorial will walk you through an analysis of mitochondrial genome data from 20 passerine birds, representing the three suborders of this highly diverse order of vertebrates. Specifically, the data set comprises the first and second codon sites of the 13 protein-coding genes that are found in the mitochondrial genome. The analyses will include two different configurations of clock models to account for evolutionary rate variation in the data. 



## The data 

Passerines (Order Passeriformes), also known as perching birds, comprise nearly 60% of all living bird species. They have diversified into a broad range of ecological roles and can be found in almost all habitats across the world ([Figure 1](#passerines)). Mitochondrial genomes have been sequenced for hundreds of passerine species. 


<figure>
	<a id="passerines"></a>
	<img src="figures/passerines.png" alt="passerines">
	<figcaption>Figure 1: Rifleman (<i>Acanthisitta chloris</i>; <i>Acanthisitti</i>), Fairy Pitta (<i>Pitta nympha</i>, <i>Tyranni</i>), and American Robin (<i>Turdus migratorius</i>; <i>Passeri</i>), representing the three major lineages of passerine birds. Creative Commons photographs by Brian Ralphs, Jason Thompson, Rhododendrites.</figcaption>
</figure>
<br>


The data files for this tutorial can be downloaded from the left-hand panel. The data set in this tutorial includes publicly available nucleotide sequences that were assembled for a study of evolutionary rate variation in mitochondrial genomes from passerine birds {% cite Nguyen2016 --file Molecular-clock-dating-and-modelling-rate-variation/master-refs.bib %}.

The data set comprises two separate mitochondrial sequence alignments from the 20 passerine taxa. The files contain the first codon sites (`passerines_pc1.nex`) and second codon sites (`passerines_pc2.nex`) of the 13 protein-coding genes in the mitochondrial genome. Together, these data sets comprise 7232 aligned nucleotide sites. The third codon sites have been excluded because they show signs of substitution saturation – they evolve rapidly and have undergone such a large amount of change that the evolutionary signal has been eroded. 


## Creating the analysis file with BEAUti

We will use BEAUti2 to select the priors and starting values for our analysis and save these settings into a BEAST2 XML file.

> Begin by starting **BEAUti2**.


### Installing BEAST2 packages

Next, we need to install a BEAST2 package that will be used in this analysis. The package is called **bModelTest**.

> Open the **BEAST 2 Package Manager** by navigating to **File > Manage Packages**. Install the **bModelTest** package by selecting it and clicking the **Install/Upgrade** button ([Figure 2]()). 
>
>

<figure>
	<a id="packagemanager"></a>
	<img src="figures/packagemanager.png" alt="BEAST2 Package Manager">
	<figcaption>Figure 2: Installing bModelTest in the BEAST 2 Package Manager. </figcaption>
</figure>
<br>


After the installation of a package, the program is on your computer, but BEAUti2 is unable to load the template files for the newly installed model unless it is restarted. So, you will now need to restart BEAUti2 so that bModelTest is available. 

> Close the **BEAST 2 Package Manager** and restart **BEAUti2** to load the **bModelTest** package.




<figure>
	<a id=""></a>
	<img src="figures/.png" alt="">
	<figcaption>Figure 2: </figcaption>
</figure>
<br>




## Figures


<figure>
	<a id="fig:example1"></a>
	<img style="width:25%;" src="figures/Logo_bw.png" alt="">
	<figcaption>Figure 1: This figure is 25% of the page width.</figcaption>
</figure>


<figure>
	<a id="fig:example2"></a>
	<img style="width:10%;" src="figures/Logo_bw.png" alt="">
	<figcaption>Figure 2: This figure is only 10% of the page width.</figcaption>
</figure>



# Code

A bit of inline monospaced font can be made `like this`. Larger code blocks can be made by using the code environment:

Java:

```java
public class HelloWorld {

    public static void main(String[] args) {
        // Prints "Hello, World" to the terminal window.
        System.out.println("Hello, World");
    }

}
```

XML:

```xml
	<BirthDeathSkylineModel spec="BirthDeathSkylineModel" id="birthDeath" tree="@tree" contemp="true">
	      <parameter name="origin" id="origin" value ="100" lower="0."/>    
	      <parameter name="R0" id="R0" value="2" lower="0." dimension ="10"/>
	      <parameter name="becomeUninfectiousRate" id="becomeUninfectiousRate" value="1" lower="0." dimension ="10"/>
	      <parameter name="samplingProportion" id="samplingProportion" value="0."/>
	      <parameter name="rho" id="rho" value="1e-6" lower="0." upper="1."/>
	</BirthDeathSkylineModel>
```

R:

```R
	> myString <- "Hello, World!"
	> print (myString)
	[1] "Hello, World!"
```

# Equations

Inline equations: {% eqinline \dot{x} = \sigma(y-x) %}

Displayed equations: 
{% eq \left( \sum_{k=1}^n a_k b_k \right)^2 \leq \left( \sum_{k=1}^n a_k^2 \right) \left( \sum_{k=1}^n b_k^2 \right) %}



## Instruction boxes

Use block-quotes for step-by-step instruction that the user should perform (this will produce a framed box on the website):

> The data we have is not the data we want, and the data we need is not the data we have.
> 
> We can input **any** formatted text in here:
>
> - Even
> - Lists
>
> or equations:
>
> {% eq (x_1, \ldots, x_n) \left( \begin{array}{ccc}
      \phi(e_1, e_1) & \cdots & \phi(e_1, e_n) \\
      \vdots & \ddots & \vdots \\
      \phi(e_n, e_1) & \cdots & \phi(e_n, e_n)
    \end{array} \right)
  \left( \begin{array}{c}
      y_1 \\
      \vdots \\
      y_n
    \end{array} \right) %}






# Hyperlinks

Add links to figures like this: 

- [Figure 1](#fig:example1) is 25% of the page width.
- [Figure 2](#fig:example2) is 10% of the page width. 

Add links to external URLs like [this](http://www.google.com). 

Links to equations or different sections within the same document are a little buggy.


----

# Useful Links

- [Bayesian Evolutionary Analysis with BEAST 2](http://www.beast2.org/book.html) {% cite BEAST2book2014 --file Molecular-clock-dating-and-modelling-rate-variation/master-refs.bib %}
- BEAST 2 website and documentation: [http://www.beast2.org/](http://www.beast2.org/)
- BEAST 1 website and documentation: [http://beast.bio.ed.ac.uk](http://beast.bio.ed.ac.uk)
- Join the BEAST user discussion: [http://groups.google.com/group/beast-users](http://groups.google.com/group/beast-users) 

----

# Relevant References

{% bibliography --cited --file Molecular-clock-dating-and-modelling-rate-variation/master-refs.bib %}
