#! /usr/bin/env python
"""

"""
import sys
import subprocess
import os
import types
import gzip
import yaml
import ConfigParser
import io
from datetime import datetime

class snp():

    def __init__(self, input, outdir, reference, name, paired,input2, verbose, argString):
        self.name               = name
        self.fOut               = outdir
        self.flog               = "QC"
        self.input              = input
        self.outdir             = self.fOut + "/tmp"
        self.tmp                = self.outdir + "/tmp"
        self.trimmomatic        = self.fOut + "/trimmomatic"
        self.paired             = paired
        self.input2             = input2
        self.verbose            = verbose
        self.reference          = reference
        self.__finalVCF         = ''
        self.__annotation       = ''
        self.__final_annotation = ''
        self.__unclear          = ''
        self.__mixed            = ''
        self.__low              = ''
        self.__exception        = ''

        # Create the output directory, and start the log file.
        self.__logged = False
        self.__CallCommand('mkdir', ['mkdir', self.fOut])

        if not os.path.isfile(self.flog):
           self.__CallCommand('mkdir', ['mkdir', self.flog])

        self.__CallCommand('mkdir', ['mkdir', '-p', self.tmp])
        self.__CallCommand('mkdir', ['mkdir', '-p', self.trimmomatic])
        self.__log     = self.fOut + "/" + self.name + ".log"
        cwd = os.getcwd()
        with open(os.path.join(os.path.dirname(__file__), "config.yml"), 'r') as ymlfile:
             cfg       = yaml.load(ymlfile)
        self.__lineage = self.fOut + "/" + self.name + ".lineage_report.txt"
        self.__logFH   = open(self.__log, 'w')
        self.__logFH.write(argString + "\n\n")
        self.__mlog    = self.flog + "/" + "master.log"
        self.__qlog    = self.flog + "/" + "low_quals.txt"
        self.__logFH2  = open(self.__mlog, 'a')
        self.__logged  = True
				
	# Format Validation
        self.__pigz               = cfg['tools']['pigz']
        self.__unpigz             = cfg['tools']['unpigz']
        self.__trimmomatic        = cfg['tools']['trimmomatic']
        # Mapping
        self.__bwa                = cfg['tools']['bwa']
        self.__samtools           = cfg['tools']['samtools']
        # Picard-Tools
        self.__picard             = cfg['tools']['picard']
        # SNP / InDel Calling
        self.__gatk               = cfg['tools']['gatk']
        # Other
        self.__bcftools           = cfg['tools']['bcftools']
        self.__bedtools           = cfg['tools']['bedtools']
        self.__vcfannotate        = cfg['tools']['vcfannotate']
        self.__vcftools           = cfg['tools']['vcftools']
        self.__vcfutils           = cfg['tools']['vcfutils'] 
        self.__annotator          = cfg['tools']['annotator'] 
        self.__parser             = cfg['scripts']['parser']
        self.__lineage_parser     = cfg['scripts']['lineage_parser']
        self.__vcf_parser         = cfg['scripts']['vcf_parser']
        self.__merger             = cfg['scripts']['merger']
        self.__lineages           = cfg['scripts']['lineages']
        self.__excluded           = cfg['scripts']['excluded']
        self.__coverage_estimator = cfg['scripts']['coverage_estimator']
        self.__bedlist_one        = cfg['scripts']['bedlist_one']
        self.__bedlist_two        = cfg['scripts']['bedlist_two']
        self.__bedlist_merge      = cfg['scripts']['bedlist_merge']     
        self.__resis_parser       = cfg['scripts']['resis_parser']
        self.mutationloci         = cfg['scripts']['mutationloci']
        self.snplist              = cfg['scripts']['snplist']
        self.__threads            = cfg['other']['threads']

    """ Shell Execution Functions """
    def __CallCommand(self, program, command):
        """ Allows execution of a simple command. """
        out = ""
        err = ""
        p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out,err = p.communicate()

        if (type(program) is list):
            o = open(program[1], 'w')
            o.write(out)
            o.close()
            out = ""
            program = program[0]
        
        if (self.__logged):
            self.__logFH.write('---[ '+ program +' ]---\n')
            self.__logFH.write('Command: \n' + ' '.join(command) + '\n\n')
            if out:
                self.__logFH.write('Standard Output: \n' + out + '\n\n')
            if err:
                self.__logFH.write('Standard Error: \n' + err + '\n\n')
        return 1


    """ QC Trimmomatic """
    def runTrimmomatic(self):
        self.__ifVerbose("Performing trimmomatic trimming.")
        if self.paired:
           self.__CallCommand('trimmomatic', ['java', '-jar', self.__trimmomatic, 'PE', '-threads', self.__threads, '-trimlog', self.trimmomatic + "/" + 'trimLog.txt', self.input, self.input2,
                              self.trimmomatic + "/" + self.name + '_paired_1.fastq.gz', self.trimmomatic + "/" + self.name + '_unpaired_1.fastq.gz',
                              self.trimmomatic + "/" + self.name + '_paired_2.fastq.gz', self.trimmomatic + "/" + self.name + '_unpaired_2.fastq.gz',
                              'LEADING:3', 'TRAILING:3', 'SLIDINGWINDOW:4:15', 'MINLEN:40'])
        else:
           self.__CallCommand('trimmomatic', ['java', '-jar', self.__trimmomatic, 'SE', '-threads', self.__threads, '-trimlog', self.trimmomatic + "/" + 'trimLog.txt',
                              self.input, self.trimmomatic + "/" + self.name + '_paired.fastq.gz',
                              'LEADING:3', 'TRAILING:3', 'SLIDINGWINDOW:4:15', 'MINLEN:40'])
        if self.paired:
           self.__CallCommand('rm', ['rm', self.trimmomatic + "/" + self.name + "_unpaired_1.fastq.gz",
                              self.trimmomatic + "/" + self.name + "_unpaired_2.fastq.gz"])
           self.input  = self.trimmomatic + "/" + self.name + "_paired_1.fastq.gz"
           self.input2 = self.trimmomatic + "/" + self.name + "_paired_2.fastq.gz"
        else:
           self.input = self.trimmomatic + "/" + self.name + "_paired.fastq.gz"
    
    """ Aligners """ 
    def runBWA(self, bwa):
        """ Align reads against the reference using bwa."""
        self.__ranBWA = True
        self.__ifVerbose("Running BWA.")
        self.__logFH.write("########## Running BWA. ##########\n")
        bwaOut = self.outdir + "/bwa"
        self.__CallCommand('mkdir', ['mkdir', '-p', bwaOut])
        self.__ifVerbose("   Building BWA index.")
        self.__bwaIndex(bwaOut + "/index")
        self.__alnSam = bwaOut + "/bwa.sam"
        self.__bwaLongReads(bwaOut)
        self.__ifVerbose("") 
        self.__processAlignment()
          
    def __bwaIndex(self, out):
        """ Make an index of the given reference genome. """ 
        self.__CallCommand('mkdir', ['mkdir', '-p', out])
        self.__CallCommand('cp', ['cp', self.reference, out + "/ref.fa"])
        self.reference = out + "/ref.fa"
        self.__CallCommand('bwa index', [self.__bwa, 'index', self.reference])
        self.__CallCommand('CreateSequenceDictionary', ['java', '-jar', self.__picard, 
                           'CreateSequenceDictionary', 'R='+self.reference,'O='+ out + "/ref.dict"])
        self.__CallCommand('samtools faidx', [self.__samtools, 'faidx', self.reference ])

    def __bwaLongReads(self, out):
        """ Make use of bwa mem """
        if self.paired:
            self.__ifVerbose("   Running BWA mem on paired end reads.")
            self.__CallCommand(['bwa mem', self.__alnSam], [self.__bwa, 'mem','-t',self.__threads,'-R', 
                               "@RG\tID:" + self.name + "\tSM:" + self.name + "\tPL:ILLUMINA", 
                                self.reference, self.input, self.input2])
        else:
            self.__ifVerbose("   Running BWA mem on single end reads.")
            self.__CallCommand(['bwa mem', self.__alnSam], [self.__bwa, 'mem','-t', self.__threads, '-R', 
                               "@RG\tID:" + self.name + "\tSM:" + self.name + "\tPL:ILLUMINA", 
                                self.reference, self.input])       

    def __processAlignment(self):
        """ Filter alignment using GATK and Picard-Tools """
        self.__ifVerbose("Filtering alignment with GATK and Picard-Tools.")
        self.__logFH.write("########## Filtering alignment with GATK and Picard-Tools. ##########\n")
        GATKdir = self.outdir + "/GATK"
        self.__CallCommand('mkdir', ['mkdir', '-p', GATKdir])

        """ Convert SAM to BAM"""
        if (self.__ranBWA):
            self.__ifVerbose("   Running SamFormatConverter.")
            self.__CallCommand('SamFormatConverter', ['java', '-Xmx4g', '-jar', self.__picard, 'SamFormatConverter',  
                               'INPUT='+ self.__alnSam, 'VALIDATION_STRINGENCY=LENIENT', 
                               'OUTPUT='+ GATKdir +'/GATK.bam', ])
        else:
            self.__CallCommand('cp', ['cp', self.__alnSam, GATKdir +'/GATK.bam'])


        """ Run mapping Report and Mark duplicates using Picard-Tools"""
        self.__ifVerbose("   Running SortSam.")
        self.__CallCommand('SortSam', ['java', '-Xmx8g', '-Djava.io.tmpdir=' + self.tmp, '-jar', self.__picard, 'SortSam',  
                           'INPUT='+ GATKdir +'/GATK.bam', 'SORT_ORDER=coordinate', 'OUTPUT='+ GATKdir +'/GATK_s.bam', 
                           'VALIDATION_STRINGENCY=LENIENT', 'TMP_DIR=' + self.tmp])
        self.__ifVerbose("   Running MarkDuplicates.")
        self.__CallCommand('MarkDuplicates', ['java', '-Xmx8g', '-jar', self.__picard, 'MarkDuplicates',  
                           'INPUT='+ GATKdir +'/GATK_s.bam', 'OUTPUT='+ GATKdir +'/GATK_sdr.bam',
                           'METRICS_FILE='+ GATKdir +'/MarkDupes.metrics', 'ASSUME_SORTED=true', 
                           'REMOVE_DUPLICATES=false', 'VALIDATION_STRINGENCY=LENIENT'])
        self.__ifVerbose("   Running BuildBamIndex.")
        self.__CallCommand('BuildBamIndex', ['java', '-Xmx8g', '-jar', self.__picard, 'BuildBamIndex',  
                           'INPUT='+ GATKdir +'/GATK_sdr.bam', 'VALIDATION_STRINGENCY=LENIENT'])

        """ Re-alignment around InDels using GATK """
        self.__ifVerbose("   Running RealignerTargetCreator.")
        self.__CallCommand('RealignerTargetCreator', ['java', '-Xmx32g', '-jar', self.__gatk, '-T', 
                           'RealignerTargetCreator', '-I', GATKdir +'/GATK_sdr.bam', '-R', self.reference, 
                           '-o', GATKdir +'/GATK.intervals', '-nt', '12'])
        self.__ifVerbose("   Running IndelRealigner.")
        self.__CallCommand('IndelRealigner', ['java', '-Xmx4g', '-jar', self.__gatk, '-T', 'IndelRealigner', '-l', 
                           'INFO', '-I', GATKdir +'/GATK_sdr.bam', '-R', self.reference, '-targetIntervals', 
                           GATKdir +'/GATK.intervals', '-o', GATKdir +'/GATK_sdrc.bam'])
        self.__ifVerbose("   Running BaseRecalibrator.")
        self.__CallCommand('BaseRecalibrator', ['java', '-Xmx16g', '-jar', self.__gatk, '-T', 'BaseRecalibrator', 
                           '-I', GATKdir +'/GATK_sdrc.bam', '-R', self.reference, '--knownSites', 
                           self.snplist,'--maximum_cycle_value', '1600', '-o', GATKdir +'/GATK_Resilist.grp','-nct', '8'])
        self.__ifVerbose("   Running PrintReads.")
        self.__CallCommand('PrintReads', ['java', '-Xmx4g', '-jar', self.__gatk, '-T', 'PrintReads', 
                           '-I', GATKdir +'/GATK_sdrc.bam', '-R', self.reference, '-BQSR', 
                           GATKdir +'/GATK_Resilist.grp', '-o', GATKdir +'/GATK_sdrcr.bam','-nct', '8'])
        self.__ifVerbose("   Running SortSam.")
        self.__CallCommand('SortSam', ['java', '-Xmx8g', '-Djava.io.tmpdir=' + self.tmp, '-jar', self.__picard,'SortSam',  
                           'INPUT='+ GATKdir +'/GATK_sdrcr.bam', 'SORT_ORDER=coordinate', 'TMP_DIR=' + self.tmp, 
                           'OUTPUT='+ GATKdir +'/GATK_sdrcs.bam', 'VALIDATION_STRINGENCY=LENIENT'])
        self.__ifVerbose("   Running BuildBamIndex.")
        self.__CallCommand('BuildBamIndex', ['java', '-Xmx8g', '-jar', self.__picard, 'BuildBamIndex', 
                           'INPUT='+ GATKdir +'/GATK_sdrcs.bam', 'VALIDATION_STRINGENCY=LENIENT'])

        """ Filter out unmapped reads """
        self.__finalBam = self.fOut + '/'+ self.name + '_sdrcsm.bam'
        self.__ifVerbose("   Running samtools view.")
        self.__CallCommand('samtools view', [self.__samtools, 'view', '-bhF', '4', '-o', self.__finalBam, 
                           GATKdir +'/GATK_sdrcs.bam'])
        self.__ifVerbose("   Running BuildBamIndex.")
        self.__CallCommand('BuildBamIndex', ['java', '-Xmx8g', '-jar', self.__picard, 'BuildBamIndex', 'INPUT='+ self.__finalBam, 
                           'VALIDATION_STRINGENCY=LENIENT'])
        self.__ifVerbose("")
        self.__CallCommand('rm', ['rm', '-r', self.tmp])
    
    """ Callers """

    def runGATK(self):
        if os.path.isfile(self.__finalBam):
            self.__ifVerbose("Calling SNPs/InDels with GATK.")
            self.__logFH.write("########## Calling SNPs/InDels with HaplotypeCaller. ##########\n")
            GATKdir = self.outdir + "/GATK"
            samDir = self.outdir + "/SamTools"
            self.__CallCommand('mkdir', ['mkdir', '-p', samDir])

            """ Call SNPs/InDels with VarDict """
            self.__ifVerbose("   Running HaplotypeCaller.")
            self.__CallCommand('HaplotypeCaller', ['java', '-Xmx4g', '-jar', self.__gatk, '-T', 'HaplotypeCaller',
                               '-R', self.reference, '-I', self.__finalBam, '-o',  GATKdir +'/gatk.g.vcf',
                               '--emitRefConfidence', 'GVCF',  '-stand_call_conf', '20.0', '-nct', '4'])
            self.__CallCommand('GenotypeGVCFs', ['java', '-Xmx4g', '-jar', self.__gatk, '-T', 'GenotypeGVCFs',
                               '-R', self.reference,'--variant', GATKdir +'/gatk.g.vcf', '-o', GATKdir +'/gatk.vcf',])
            self.__CallCommand('ReadBackedPhasing', ['java', '-Xmx4g', '-jar', self.__gatk, '-T', 'ReadBackedPhasing',
                               '-R', self.reference, '-I', self.__finalBam, '-V',  GATKdir +'/gatk.vcf', '-o', GATKdir +'/gatk_phased.vcf',
                               '--phaseQualityThresh', '20.0',  '-enableMergeToMNP', '-maxDistMNP', '2'])
            self.__CallCommand(['merge vcfs', self.fOut + "/" + self.name +'_gatk.vcf'],
                               [self.__merger, GATKdir +'/gatk.vcf', GATKdir +'/gatk_phased.vcf'])  
            self.__CallCommand(['vcf-annotate filter', self.fOut + "/" + self.name +'_gatk_marked.vcf'], 
                               [self.__vcfannotate, '--filter', 'SnpCluster=3,10/Qual=20/MinDP=10/MinMQ=20', self.fOut + "/" + self.name +'_gatk.vcf'])
            self.__CallCommand(['vcftools remove-filtered-all', self.fOut + "/" + self.name +'_gatk_filtered.vcf'], 
                                   [self.__vcftools, '--vcf', self.fOut + "/" + self.name +'_gatk_marked.vcf',
                                   '--stdout', '--exclude-bed', self.__excluded, '--remove-filtered-all', '--recode', '--recode-INFO-all'])
            self.__CallCommand(['samtools depth', samDir + '/coverage.txt'],
                                [self.__samtools,'depth', self.__finalBam])
            self.__CallCommand(['bedtools coverage', samDir + '/bed_1_coverage.txt' ],
                                [self.__bedtools, 'coverage', '-abam', self.__finalBam, '-b', self.__bedlist_one])
            self.__CallCommand(['bedtools coverage', samDir + '/bed_2_coverage.txt' ],
                                [self.__bedtools, 'coverage', '-abam', self.__finalBam, '-b', self.__bedlist_two])
            self.__CallCommand(['sort', samDir + '/bed_1_sorted_coverage.txt' ],['sort', '-nk', '6', samDir + '/bed_1_coverage.txt'])
            self.__CallCommand(['sort', samDir + '/bed_2_sorted_coverage.txt' ],['sort', '-nk', '6', samDir + '/bed_2_coverage.txt'])

            """ Set final VCF file. """
            
            if not self.__finalVCF: 
                self.__finalVCF = self.fOut + "/" + self.name +'_gatk_filtered.vcf'
        else:
            # print error
            pass  
       
    def annotateVCF(self):
        """ Annotate the final VCF file """
        cwd = os.getcwd()
        if self.__finalVCF:
           self.__ifVerbose("Annotating final VCF.")
           self.__CallCommand(['SnpEff', self.fOut + "/" + self.name +'_annotation.txt'],
                                ['java', '-Xmx4g', '-jar', self.__annotator, 'NC_000962', self.__finalVCF])
           self.__annotation = self.fOut + "/" + self.name +'_annotation.txt'
           self.__ifVerbose("parsing final Annotation.")
           self.__CallCommand(['parse annotation', self.fOut + "/" + self.name +'_Final_annotation.txt'],
                               ['python', self.__parser, self.__annotation, self.name, self.mutationloci])
        else:
            self.__ifVerbose("Use SamTools, GATK, or Freebayes to annotate the final VCF.")
        self.__CallCommand('rm', ['rm',  cwd + "/snpEff_genes.txt"])
        self.__CallCommand('rm', ['rm',  cwd + "/snpEff_summary.html"])

    def runLineage(self):
        """ Run lineage Analysis """
        self.__ifVerbose("Running Lineage Analysis")
        self.__final_annotation = self.fOut + "/" + self.name +'_Final_annotation.txt'
        self.__CallCommand(['lineage parsing', self.fOut + "/" + self.name +'_Lineage.txt'],
                              ['python', self.__lineage_parser, self.__lineages, self.__final_annotation, self.__lineage, self.name])
        count1 = 0
        count2 = 0
        count3 = 0
        if os.path.isfile(self.__lineage):
           fh1 = open(self.__lineage,'r')   
           for line in fh1:
               lined = line.rstrip("\r\n")
               i = datetime.now()
               if "No Informative SNPs" in lined:
                   self.__logFH2.write(i.strftime('%Y/%m/%d %H:%M:%S') + "\t" + "Input:" + "\t" + self.name + "\t" + "no clear lineage classification\n")
                   self.__unclear = "positive"
               elif "no precise lineage" in lined or "mixed lineage(s)" in lined:
                   self.__logFH2.write(i.strftime('%Y/%m/%d %H:%M:%S') + "\t" + "Input:" + "\t" + self.name + "\t" + "no clear lineage classification\n")
                   self.__mixed = "positive"
               elif "no concordance" in lined:
                   self.__logFH2.write(i.strftime('%Y/%m/%d %H:%M:%S') + "\t" + "Input:" + "\t" + self.name + "\t" + "no clear lineage classification\n")
                   self.__mixed = "positive"
               elif len(lined) < 3:
                  self.__logFH2.write(i.strftime('%Y/%m/%d %H:%M:%S') + "\t" + "Input:" + "\t" + self.name + "\t" + "no clear lineage classification\n")
                  self.__unclear = "positive" 
           fh1.close()
        if os.path.isfile(self.fOut + "/" + self.name +'_Final_annotation.txt'):
           fh2 = open(self.fOut + "/" + self.name +'_Final_annotation.txt','r')
           for lines in fh2:
              lined = lines.rstrip("\r\n").split("\t")
              if lined[16] == "rrs":
                 count1 += 1
              elif lined[16] == "rrl":
                 count2 += 1
           if count1 > 5 or count2 > 5 :
              self.__logFH2.write(i.strftime('%Y/%m/%d %H:%M:%S') + "\t" + "Input:" + "\t" + self.name + "\t" + "mixed species suspected\n")
              self.__mixed = "positive"
           fh2.close()
        if os.path.isfile(self.fOut + "/" + self.name +'_Final_annotation.txt'):
           if os.path.isfile(self.fOut + "/" + self.name +'_GATK.vcf'):
              self.__CallCommand('vcf_parser', ['python', self.__vcf_parser, self.fOut + "/" + self.name +'_Final_annotation.txt',
                                                self.fOut + "/" + self.name +'_GATK.vcf', self.__qlog])
           elif os.path.isfile(self.fOut + "/" + self.name +'_SamTools.vcf'):
               self.__CallCommand('vcf_parser', ['python', self.__vcf_parser, self.fOut + "/" + self.name +'_Final_annotation.txt',
                                                self.fOut + "/" + self.name +'_SamTools.vcf', self.__qlog])

    def runCoverage(self):
        """ Run Genome Coverage Statistics """
        cov = ""
        notes = []
        dele = True
        wid  = ""

        self.__ifVerbose("Running Genome Coverage Statistics")
        samDir = self.outdir + "/SamTools"
        i = datetime.now()
        self.__CallCommand(['coverage estimator', self.fOut + "/" + self.name + '_Coverage.txt'],
                            ['python', self.__coverage_estimator, samDir + '/coverage.txt', self.name])
        self.__CallCommand(['genome region coverage estimator', samDir + '/genome_region_coverage_1.txt'],
                            ['python', self.__resis_parser, samDir + '/bed_1_sorted_coverage.txt', samDir + '/coverage.txt', self.name])
        self.__CallCommand(['genome region coverage estimator', samDir + '/genome_region_coverage_2.txt'],
                            ['python', self.__resis_parser, samDir + '/bed_2_sorted_coverage.txt', samDir + '/coverage.txt', self.name])
        self.__CallCommand(['cat' , samDir + '/genome_region_coverage.txt'],['cat', samDir + '/genome_region_coverage_1.txt', samDir + '/genome_region_coverage_2.txt'])
        self.__CallCommand(['sort', self.fOut + "/" + self.name + '_genome_region_coverage.txt' ],['sort', '-nk', '3', samDir + '/genome_region_coverage.txt'])
        self.__CallCommand('sed',['sed', '-i', '1d', self.fOut + "/" + self.name + '_genome_region_coverage.txt'])
        
        if os.path.isfile(self.fOut + "/" + self.name + '_Coverage.txt'):
          fh2 = open(self.fOut + "/" + self.name + '_Coverage.txt','r')
          for line in fh2:
             if line.startswith("Sample"):
                continue
             cov_str = line.split("\t")
             cov = cov_str[1]
             wid = cov_str[2]
          if cov != '' and int(cov) < 10:
            self.__low = "positive"
            self.__logFH2.write(i.strftime('%Y/%m/%d %H:%M:%S') + "\t" + "Input:" + "\t" + self.name + "\t" + "low genome coverage depth\n")
          if wid != '' and float(wid) < 94.99:
            self.__low = "positive"
            self.__logFH2.write(i.strftime('%Y/%m/%d %H:%M:%S') + "\t" + "Input:" + "\t" + self.name + "\t" + "low genome coverage width\n")
          fh2.close()                   
         
    def cleanUp(self):
        """ Clean up the temporary files, and move them to a proper folder. """
        i = datetime.now()
        self.__CallCommand('rm', ['rm', '-r', self.outdir])
        self.__CallCommand('rm', ['rm',  self.fOut + "/" + self.name +'_annotation.txt'])
        self.__CallCommand('rm', ['rm',  self.fOut + '/'+ self.name + '_sdrcsm.bai'])
        self.__CallCommand('rm', ['rm',  self.__finalBam])
        self.__CallCommand('chmod', ['chmod', '777',  self.fOut])
        if self.__mixed == "positive":
           self.__CallCommand('mv', ['mv', self.fOut, self.flog])
        if self.__unclear == "positive":
           self.__CallCommand('mv', ['mv', self.fOut, self.flog])
        if self.__low == "positive":
           self.__CallCommand('mv', ['mv', self.fOut, self.flog])
        if os.path.isfile(self.fOut + "/" + self.name + '.log'):
           self.__logFH.close()
           self.__logged = False
           fh4 = open(self.__log,'r')
           for line in fh4:
               lines = line.rstrip("\r\n")
               if "Exception" in lines:
                  self.__exception = "positive"
           fh4.close()        
        if self.__exception == "positive":
           self.__logFH2.write(i.strftime('%Y/%m/%d %H:%M:%S') + "\t" + "Input:" + "\t" + self.name + "\t" + "Exception in analysis\n")
           self.__CallCommand('mv', ['mv', self.fOut, self.flog])
    def __ifVerbose(self, msg):
        """ If verbose print a given message. """
        if self.verbose: print msg
