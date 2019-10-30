#!/usr/bin/awk -f

######################################################################
##        This is a awk code written by Sujan Mamidi                ##
##                                                                  ## 
##               This code is licensed under the                    ##
##                GNU General Public License v3.0                   ##
##                                                                  ##
## This code reads a vcf file (SNP called from population Bulks) to ##
##                      determine the type of site                  ##
## Fixed site - Allele frequency is 100% for each of the pop, but   ##
## 					Alternate alleles
## Shared site - Both alleles resent in both the populations        ##
## Unique_pop - Allele frequency of one allele is 100% in only      ##
##   one population, but other population has both alleles          ##
##                                                                  ##
## The output file can be used to look for frequency distribution   ##
##                       based on site or delta                     ##
######################################################################

## Usage:
## Uncompressed file : cat my.vcf | awk -f identify_snpType.awk > my.vcf.sitetype
## bz2 compressed file : bzip2 -dck my.vcf.bz2 | awk -f identify_snpType.awk > my.vcf.sitetype


{ !/^#/ ;

split($10,a,":");
split($11,b,":");
split(a[2],c,",");
split(b[2],d,",");

pop1Ref =c[1];
pop1Var =c[2];

pop2Ref =d[1];
pop2Var =d[2];

if((pop1Ref + pop1Var) !=0) pop1_reffreq = pop1Ref/(pop1Ref + pop1Var);
if((pop2Ref + pop2Var) !=0) pop2_reffreq = pop2Ref/(pop2Ref + pop2Var);

delta=pop1_reffreq-pop2_reffreq;
(delta<0)?-delta:delta;

if(delta==0) site="NotAnSNP";
else if(delta==1) site="fixed";
else if(pop1_reffreq=0 && pop2_reffreq>0) site = "unique_pop2";
else if(pop2_reffreq=0 && pop1_reffreq>0) site = "unique_pop1";
else if(pop1_reffreq=1 && pop2_reffreq>0) site = "unique_pop2";
else if(pop2_reffreq=1 && pop1_reffreq>0) site = "unique_pop1";
else site="shared";
if((c[1]+c[2]) !=0 && (d[1]+d[2]) !=0) print $1,$2,$4,$5,delta,site
}
