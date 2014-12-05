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
BINDIR=$(abspath ${OUTDIR})/bin


# if tools are undefined
bwa.exe ?=${BINDIR}/bwa-0.7.10/bwa
samtools.exe ?=${BINDIR}/samtools-0.1.19/samtools
bcftools.exe ?=${BINDIR}/samtools-0.1.19/bcftools/bcftools
tabix.exe ?=${BINDIR}/tabix-0.2.6/tabix
bgzip.exe ?=${BINDIR}/tabix-0.2.6/bgzip
java.exe ?= java
picard.version?=1.119
picard.dir=${BINDIR}/picard-tools-${picard.version}

.PHONY= all clean all_bams all_vcfs

all: all_vcfs


all_vcfs: <xsl:for-each select="project"> \
	<xsl:apply-templates select="." mode="vcf.final"/></xsl:for-each>


all_bams: <xsl:for-each select="project/sample"> \
	<xsl:apply-templates select="." mode="bam.final"/></xsl:for-each>


<xsl:apply-templates select="project"/>

$(addsuffix .fai,${REFERENCE}): ${REFERENCE} ${samtools.exe}
	${samtools.exe} faidx $&lt;

$(addsuffix .bwt,${REFERENCE}): ${REFERENCE} ${bwa.exe}
	${bwa.exe} index $&lt;


${BINDIR}/bwa-0.7.10/bwa :
	rm -rf $(BINDIR)/bwa-0.7.10/ &amp;&amp; \
	mkdir -p $(BINDIR) &amp;&amp; \
	curl -o $(BINDIR)/bwa-0.7.10.tar.bz2 -L "http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.10.tar.bz2/download?use_mirror=freefr" &amp;&amp; \
	tar xvfj $(BINDIR)/bwa-0.7.10.tar.bz2 -C $(OUTDIR)/bin  &amp;&amp; \
	rm $(BINDIR)/bwa-0.7.10.tar.bz2 &amp;&amp; \
	make -C $(dir $@)

${BINDIR}/samtools-0.1.19/bcftools/bcftools: ${BINDIR}/samtools-0.1.19/samtools

${BINDIR}/samtools-0.1.19/samtools  :	
	rm -rf $(BINDIR)/samtools-0.1.19/ &amp;&amp; \
	mkdir -p $(BINDIR) &amp;&amp; \
	curl -o $(BINDIR)/samtools-0.1.19.tar.bz2 -L "http://sourceforge.net/projects/samtools/files/samtools-0.1.19.tar.bz2/download?use_mirror=freefr" &amp;&amp; \
	tar xvfj $(BINDIR)/samtools-0.1.19.tar.bz2 -C $(BINDIR)  &amp;&amp; \
	rm $(BINDIR)/samtools-0.1.19.tar.bz2 &amp;&amp; \
	make -C $(dir $@)


${BINDIR}/tabix-0.2.6/bgzip : ${BINDIR}/tabix-0.2.6/tabix


${BINDIR}/tabix-0.2.6/tabix  :	
	rm -rf $(BINDIR)/tabix-0.2.6/ &amp;&amp; \
	mkdir -p $(BINDIR) &amp;&amp; \
	curl -o $(BINDIR)/tabix-0.2.6.tar.bz2 -L "http://sourceforge.net/projects/samtools/files/tabix-0.2.6.tar.bz2/download?use_mirror=freefr" &amp;&amp; \
	tar xvfj $(BINDIR)/tabix-0.2.6.tar.bz2 -C $(BINDIR)  &amp;&amp; \
	rm $(BINDIR)/tabix-0.2.6.tar.bz2 &amp;&amp; \
	make -C $(dir $@) tabix bgzip

${BINDIR}/picard-tools-1.119/MergeSamFiles.jar : ${BINDIR}/picard-tools-1.119/MarkDuplicates.jar
${BINDIR}/picard-tools-1.119/MarkDuplicates.jar :
	rm -rf $(dir $@) &amp;&amp; \
	mkdir -p $(BINDIR) &amp;&amp; \
	curl -L -k -o ${BINDIR}/picard-tools-1.119.zip -L "http://downloads.sourceforge.net/project/picard/picard-tools/1.119/picard-tools-1.119.zip?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fpicard%2Ffiles%2Fpicard-tools%2F1.119%2F&amp;ts=1417770270&amp;use_mirror=skylink" &amp;&amp; \
	unzip ${BINDIR}/picard-tools-1.119.zip -d ${BINDIR} &amp;&amp; \
	rm $(BINDIR)/picard-tools-1.119.zip



clean:
	rm -rf ${BINDIR}

</xsl:template>


<xsl:template match='project'>

#
# VCF for project '<xsl:value-of select="@name"/>'
# 
<xsl:apply-templates select="." mode="vcf.final"/> : $(addsuffix .bai,<xsl:for-each select="sample"><xsl:text> </xsl:text><xsl:apply-templates select="." mode="bam.final"/> </xsl:for-each>) \
	$(addsuffix .fai,${REFERENCE}) ${samtools.exe}  ${bgzip.exe} ${tabix.exe} ${bcftools.exe}
	mkdir -p $(dir $@) &amp;&amp; \
	${samtools.exe} mpileup -uf ${REFERENCE} $(basename $(filter %.bai,$^)) | \
	${bcftools.exe} view -vcg - > $(basename $@)  &amp;&amp; \
	${bgzip.exe} -f $(basename $@) &amp;&amp; \
	${tabix.exe} -f -p vcf $@ 
	
	   

<xsl:apply-templates select="sample"/>

</xsl:template>


<xsl:template match='sample'>

#
# index final BAM for Sample '<xsl:value-of select="@name"/>'
# 
$(addsuffix .bai, <xsl:apply-templates select="." mode="bam.final"/>): <xsl:apply-templates select="." mode="bam.final"/> ${samtools.exe}
	mkdir -p $(dir $@) &amp;&amp; \
	${samtools.exe} index  $&lt;
#
# prepare final BAM for Sample '<xsl:value-of select="@name"/>'
# 
<xsl:apply-templates select="." mode="bam.final"/><xsl:text> : </xsl:text><xsl:apply-templates select="." mode="bam.rmdup"/>
	mkdir -p $(dir $@) &amp;&amp; \
	cp  $&lt; $@


<xsl:apply-templates select="." mode="bam.rmdup"/><xsl:text> : </xsl:text><xsl:apply-templates select="." mode="bam.merged"/> ${picard.dir}/MarkDuplicates.jar
	mkdir -p $(dir $@) &amp;&amp; \
	${java.exe} -jar $(filter %.jar,$^) I=$&lt; O=$@ M=$(addsuffix .metrics,$@) AS=true VALIDATION_STRINGENCY=SILENT


<xsl:if test="count(fastq) &gt; 1">
#
# merge BAMs 
#
<xsl:apply-templates select="." mode="bam.merged"/><xsl:text> : </xsl:text><xsl:for-each select="fastq"> \
 	<xsl:apply-templates select="." mode="bam.sorted"/> </xsl:for-each> ${picard.dir}/MergeSamFiles.jar
	mkdir -p $(dir $@) &amp;&amp; \
 	${java.exe} -jar $(filter %.jar,$^)  $(foreach B,$(filter %.bam,$^), I=${B} ) O=$@ AS=true VALIDATION_STRINGENCY=SILENT
  
</xsl:if>


<xsl:apply-templates select="fastq"/>

</xsl:template>


<xsl:template match='fastq'>


	
#
# Index BAM <xsl:apply-templates select="." mode="bam.sorted"/>
#
<xsl:text>$(addsuffix .bai,</xsl:text><xsl:apply-templates select="." mode="bam.sorted"/><xsl:text> ): </xsl:text><xsl:apply-templates select="." mode="bam.sorted"/> ${samtools.exe}
	${samtools} index $&lt;

#
# Align <xsl:value-of select="for"/> and <xsl:value-of select="rev"/>
#
<xsl:apply-templates select="." mode="bam.sorted"/> : \
	<xsl:apply-templates select="for"/>  \
	<xsl:apply-templates select="rev"/> \
	$(addsuffix .bwt,${REFERENCE}) \
	${bwa.exe} ${samtools.exe}
	mkdir -p $(dir $@) &amp;&amp; \
	${bwa.exe} mem -R '@RG\tID:<xsl:apply-templates select="." mode="fastq.id"/>\tSM:<xsl:value-of select="../@name"/>\tLB:<xsl:choose>
		<xsl:when test="@library"><xsl:value-of select="@library"/></xsl:when>
		<xsl:otherwise><xsl:value-of select="../@name"/></xsl:otherwise>
		</xsl:choose>\tPL:<xsl:choose>
		<xsl:when test="@platform"><xsl:value-of select="@platform"/></xsl:when>
		<xsl:otherwise><xsl:text>ILLUMINA</xsl:text></xsl:otherwise>
		</xsl:choose>\tPU:<xsl:choose>
		<xsl:when test="@lane"><xsl:value-of select="@lane"/></xsl:when>
		<xsl:otherwise><xsl:text>1</xsl:text></xsl:otherwise>
		</xsl:choose><xsl:choose>
		<xsl:when test="@median-size">\tPI:<xsl:value-of select="@median-size"/></xsl:when>
		</xsl:choose>' \
		${REFERENCE} \
		<xsl:apply-templates select="for"/> \
		<xsl:apply-templates select="rev"/> |\
	${samtools.exe} view -uS - |\
	${samtools.exe} sort - $(basename $@) 


</xsl:template>

<!-- BAM mapped to reference and sorted -->
<xsl:template match='for|rev'>

<xsl:value-of select="normalize-space(./text())"/>

</xsl:template>



<!-- BAM mapped to reference and sorted -->
<xsl:template match='fastq' mode="bam.sorted">
<xsl:apply-templates select="." mode="bam.dir"/>
<xsl:value-of select="concat('${tmp.prefix}',count(preceding-sibling::fastq) + 1,'_sorted.bam')"/>
</xsl:template>

<!-- directory for the BAMs -->
<xsl:template match='fastq' mode="bam.dir">
<xsl:apply-templates select=".." mode="bam.dir"/>
</xsl:template>

<!-- fastq-id will be used for read-groups in BAM -->
<xsl:template match='fastq' mode="fastq.id">
<xsl:choose>
<xsl:when test="@id"><xsl:value-of select="@id"/></xsl:when>
<xsl:otherwise><xsl:value-of select="generate-id(.)"/></xsl:otherwise>
</xsl:choose>
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
  <xsl:when test="count(fastq)=1">
  	<xsl:apply-templates select="fastq" mode="bam.sorted"/>
  </xsl:when>
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
