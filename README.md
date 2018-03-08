# ProteinClusterQuant

*Protein Cluster Quant* PCQ is a Java software for the analysis of complex proteomics samples (quantitative or not). It helps to reduce redundancy of peptide-to-protein relationship and to visualize the results in a bipartite network (Cytoscape).  

## Use cases:
 - Working with NGS databases results in a very redundant MS dataset in which each peptide is mapped to numerous protein sequence entries (sequencing products). PCQ will help you to 1) visualize all the data in a bipartite network, removing all the redundancy of the peptide-to-protein associations, but not loosing any piece of information as other protein inference approaches do; 2) try to  
  
For a detailed description of the software, input files, intput parameters and outputs go to our **[wiki page](https://github.com/proteomicsyates/ProteinClusterQuant/wiki)**.

### How to get ProteinClusterQuant:

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



