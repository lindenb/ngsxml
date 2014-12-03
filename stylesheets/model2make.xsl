<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet
	version="1.0"
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
	>


<xsl:output method='text' />

<xsl:template match='/'>
  <xsl:apply-templates select="model"/>
</xsl:template>

<xsl:template match='model'>
# Name
# 	<xsl:value-of select="@name"/>
# Description:
# 	<xsl:value-of select="@description"/>
# 
include config.mk
OUTDIR=<xsl:value-of select="@directory"/>

.PHONY=all clean all_bams all_vcfs

all: all_vcfs


all_vcfs: <xsl:for-each select="project"> \
	<xsl:apply-templates select="." mode="vcf.final"/></xsl:for-each>


all_bams: <xsl:for-each select="project/sample"> \
	<xsl:apply-templates select="." mode="bam.final"/></xsl:for-each>


<xsl:apply-templates select="project"/>

$(addsuffix .fai,${REFERENCE}): ${REFERENCE}
	${samtools.exe} faidx $lt;

$(addsuffix .bwt,${REFERENCE}): ${REFERENCE}
	${bwa.exe} index $lt;
	

clean:

</xsl:template>


<xsl:template match='project'>

#
# VCF for project '<xsl:value-of select="@name"/>'
# 
<xsl:apply-templates select="." mode="vcf.final"/> : <xsl:for-each select="sample"> \
	<xsl:apply-templates select="." mode="bam.final"/></xsl:for-each>
	mkdir -p $(dir $@) &amp;&amp; \
	samtools mpilup &amp;&amp; \
	bgzip $(basename$@) &amp;&amp; \
	tabix -p vcf $@
	

<xsl:apply-templates select="sample"/>

</xsl:template>


<xsl:template match='sample'>

#
# prepare final BAM for Sample '<xsl:value-of select="@name"/>'
# 
<xsl:apply-templates select="." mode="bam.final"/><xsl:text> : </xsl:text><xsl:apply-templates select="." mode="bam.rmdup"/>
	mkdir -p $(dir $@) &amp;&amp; \
	picard rmdup IN=$&lt; -O $@	

<xsl:apply-templates select="." mode="bam.rmdup"/><xsl:text> : </xsl:text><xsl:apply-templates select="." mode="bam.merged"/>
	mkdir -p $(dir $@) &amp;&amp; \
	picard rmdup IN=$&lt; -O $@


<xsl:if test="count(fastq) &gt; 1">
#
# merge BAMs 
#
<xsl:apply-templates select="." mode="bam.merged"/><xsl:text> : </xsl:text><xsl:for-each select="fastq"> \
 	<xsl:apply-templates select="." mode="bam.sorted"/> </xsl:for-each>
	mkdir -p $(dir $@) &amp;&amp; \
 	picard merge $(foreach B,$(filter %.bam,$^), I=$B ) -O $@
  
</xsl:if>


<xsl:apply-templates select="fastq"/>

</xsl:template>


<xsl:template match='fastq'>


	
#
# Index BAM <xsl:apply-templates select="." mode="bam.sorted"/>
#
<xsl:text>$(addsuffix .bai,</xsl:text><xsl:apply-templates select="." mode="bam.sorted"/><xsl:text> ): </xsl:text><xsl:apply-templates select="." mode="bam.sorted"/>
	${samtools} index $&lt;

#
# Align <xsl:value-of select="for"/> and <xsl:value-of select="rev"/>
#
<xsl:apply-templates select="." mode="bam.sorted"/> : \
	<xsl:value-of select="for"/>  \
	<xsl:value-of select="rev"/> \
	$(addsuffix .bwt,${REFERENCE})
	mkdir -p $(dir $@) &amp;&amp; \
	${bwa.exe} mem -R '@RG\tID:<xsl:value-of select="generate-id(.)"/>\tSM:<xsl:value-of select="../@name"/>' ${REFERENCE} $(filter %q.gz $^) |\
	${samtools.exe} view -uS - |\
	${samtools.exe} sort -o $(basename $@) -


</xsl:template>

<!-- BAM mapped to reference and sorted -->
<xsl:template match='fastq' mode="bam.sorted">
<xsl:apply-templates select="." mode="bam.dir"/>
<xsl:value-of select="concat(count(preceding-sibling::fastq) + 1,'_sorted.bam')"/>
</xsl:template>

<!-- directory for the BAMs -->
<xsl:template match='fastq' mode="bam.dir">
<xsl:apply-templates select=".." mode="bam.dir"/>
</xsl:template>

<!-- directory for the BAMs -->
<xsl:template match='sample' mode="bam.dir">
<xsl:apply-templates select="." mode="dir"/>
<xsl:text>/BAM/</xsl:text>
</xsl:template>

<!--Base directory for a sample -->
<xsl:template match='sample' mode="dir">
<xsl:apply-templates select=".." mode="dir"/>
<xsl:text>/Samples/</xsl:text>
<xsl:value-of select="@name"/>
</xsl:template>



<!-- Print the name BAM with merged sorted BAM -->
<xsl:template match='sample' mode="bam.merged">
<xsl:choose>
  <!-- only one pair of fastq . No need to merge -->
  <xsl:when test="count(fastq)=1">
  	<xsl:apply-templates select="." mode="bam.sorted"/>
  </xsl:when>
  <!-- merge BAMS -->
  <xsl:otherwise>
    <xsl:apply-templates select="." mode="bam.dir"/>
    <xsl:text>${tmp.prefix}merged.bam</xsl:text>
  </xsl:otherwise>
</xsl:choose>
</xsl:template>

<!-- Print the name BAM with marked duplicates for a sample -->
<xsl:template match='sample' mode="bam.rmdup">
  <xsl:apply-templates select="." mode="bam.dir"/>
  <xsl:text>${tmp.prefix}rmdup.bam</xsl:text>
</xsl:template>

<!-- Print the name of final BAM for a sample -->
<xsl:template match='sample' mode="bam.final">
  <xsl:apply-templates select="." mode="bam.dir"/>
  <xsl:value-of select="concat(../@name,'_',@name,'.bam')"/>
</xsl:template>

<!-- Print the name of the bgzipped VCF for a project. All samples merged -->
<xsl:template match='project' mode="vcf.final">
  <xsl:apply-templates select="." mode="dir"/>
  <xsl:text>/VCF/</xsl:text>
  <xsl:value-of select="concat(@name,'.vcf.gz')"/>
</xsl:template>

<!-- Print the PROJECT directory -->
<xsl:template match='project' mode="dir">
<xsl:text>$(OUTDIR)/Projects/</xsl:text>
<xsl:value-of select="@name"/>
</xsl:template>


</xsl:stylesheet>
