## Motivation

build Makefile-based workflows for Next Generation Sequencing with XML model and XSLT transformations


## Example
model of data
```xml
<?xml version="1.0" encoding="UTF-8"?>
<model name="myProject" description="my project" directory="OUT">
  <project name="Proj1">
    <sample name="Sample1">
      <fastq>
        <for>sample_1_01_R1.fastq.gz</for>
        <rev>sample_1_01_R2.fastq.gz</rev>
      </fastq>
      <fastq>
        <for>sample_1_02_R1.fastq.gz</for>
        <rev>sample_1_02_R2.fastq.gz</rev>
      </fastq>
    </sample>
    <sample name="Sample2">
      <fastq>
        <for>sample_2_01_R1.fastq.gz</for>
        <rev>sample_2_01_R2.fastq.gz</rev>
      </fastq>
    </sample>
  </project>
</model>
```

process the model with the xslt-stylesheet and create a **Makefile**

```bash
$ xsltproc --output makefile stylesheets/model2make.xsl test/model01.xml
```




## Author
Author: Pierre Lindenbaum

twitter @yokofakun
