## Segmentduplication Sequence Annotation Tool for Sniffles  
#### Version 1
A simple segmentduplication sequence annotation tool for sniffles SV result.  


### Usage 
>python3 segmentduplication_anno.py `<`vcf`>` `<`outfile`>` `<`bp`>` `<`dataset`>`  
>usage:segmentduplication_anno.py [options]  

>Segmentduplication sequence annotation for sniffles vcf  

>optional arguments:   
> `<`vcf`>`         	&nbsp;sniffles vcf file  
> `<`outfile`>`     &nbsp;output filename including path  
> `<`bp`>`          &nbsp;set the range of sequence  
> `<`dataset`>`     &nbsp;the segmentduplication sequence that you refer

### Work Flow

1. Input the sniffles SV result and range of the chrmosome position 

3. If the range of the two breakpoint of sniffles result is covered by dataset, and the strand direction is the same as the dataset, the result will be annotated   

3. Output the annotated result 



### Reference Dataset  
dataset was downloaded form [http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz ](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz)




  
