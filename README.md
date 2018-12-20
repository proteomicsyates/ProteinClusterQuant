# ProteinClusterQuant

*Protein Cluster Quant* PCQ is a Java software for the analysis of complex proteomics samples (quantitative or not). It helps to reduce redundancy of peptide-to-protein relationship and to visualize the results in a bipartite network (Cytoscape).  

## Use cases:
 - Working with NGS databases results in a very redundant MS dataset in which each peptide is mapped to numerous protein sequence entries (sequencing products). PCQ will help you to 1) visualize all the data in a bipartite network, removing all the redundancy of the peptide-to-protein associations, but not loosing any piece of information as other protein inference approaches do; 2) try to align peptides (it is actually an option in the program) in order to reconstruct the protein sequences even with some SNPs or sequencing errors.
 - If you want to quantitatively distinguish isoforms (difficult problem when proteins share many peptides), or simply if you are worried about how other software handles quant information of shared peptides.
 - If you want to apply a proper statistical aggregation of ratios ([SanXot](https://www.ncbi.nlm.nih.gov/pubmed/24512137)).
 - Even having a non-quant experiment, if you have a lot of redundancy (peptides shared by lots of proteins).
 - If you are analyzing samples from multiple organisms, you can use PCQ to reduction redundancy and use the peptide alignment option to connect and compare ortholog proteins.
 - You can now use PCQ to analyze quantitative data at site level, that is, PTM quantitative data or surface accesibility data. In this way, all the ratios assigned to a particular site in a protein will be collapsed and analyzed together.

## What does PCQ exactly do:
 - Can work with quant data (any, even isobaric isotopologues) or just spectral counts.
 - Represents peptide-to-protein associations in a bipartite graph (peptides nodes and protein nodes). Proteins sharing same set of peptides are collapsed in a single node. Peptides corresponding to the same set of proteins are collapsed in a single node.
 - Helps to reduce redundancy and to visualize data easily, without losing any piece of information (as other protein inference algorithms do).
 - Using quantitative information:
    - Analyzes protein pairs (protein nodes sharing a peptide node) to detect interesting inconsistencies.
    - It can optionally aggregate ratios using SanXot* algorithm (handles source of errors properly, from individual ratio measurements in a PSM to a peptide node and having technical and/or biological replicates, using a weight average of the ratios). It calculates a quantitative FDR.
 - When used with quantified sites, such as PTMs as phosphorilations or surface accesibility measurements, PCQ helps you to aggregate (average) the ratios obtained in different peptide sequences but containing the same site of interest (the same phosphorilation site of a protein, for example).
 - In contains a sub-module in order to compare quantified surface accessibility sites from multiple samples. [surface_accessibility](go to Comparison of surface accessibility in multiple samples).

## Detailed instructions:  
For a detailed description of the software, input files, intput parameters and outputs go to our **[wiki page](https://github.com/proteomicsyates/ProteinClusterQuant/wiki)**.

## How to get ProteinClusterQuant:

**Download latest build from [here](http://sealion.scripps.edu/PCQ)**  


**Using maven:**   
Add the dependency:  
```
<dependency>  
  <groupId>edu.scripps.yates</groupId>   
  <artifactId>ProtClusterQuant</artifactId>    
  <version>$version</version>  
</dependency>
```  

Add the repositories:  
```
<repository>    
  <id>internal</id>  
  <url>http://sealion.scripps.edu/archiva/repository/internal/</url>  
</repository>  
<snapshotRepository>  
  <id>snapshots</id>  
  <url>http://sealion.scripps.edu/archiva/repository/snapshots/</url>  
</snapshotRepository>
``` 

## Getting help
If you have any question, you want to report a bug or a suggestion, please contact to Salvador Martinez-Bartolome (salvador at scripps.edu)

