
  library(ggplot2)
  library(msa)
  data <- read.csv("Data/GISDDprimer.csv")  
  
  
  library(Biostrings)
  
  D1 <- "ACCACAAAGAAAGACTTAGGGATTGGCCATGTGGCTGTTGAAAACCACCACCATGCCGCAATGCTGGACGTGGACTTACATCCAGCTTCAGCCTGGACCCTCTATGCAGTGGCCACAACAATTATCACTCCCATGATGAGGCACACAATCGAAAACACAACGGCAAACATTTCCCTGACAGCCATTGCAAACCAGGCAGCTATATTGATGGGACTTGACAAAGGATGGCCAATATCAAAGATGGACATAGGAGTTCCACTTCTCGCCTTGGGGTGCTATTCCCAGGTGAATCCACTGACGCTGACAGCGGCGGTATTGATGCTAGTGGCTCATTACGCTATAATTGGACCTGGACTGCAAGCAAAAGCTACTAGAGAAGCTCAAAAAAGGACAGCGGCCGGAATAATGAAAAATCCAACCGTTGATGGAATCGTTGCAATAGATTTGGACCCTGTGGTTTATGATGCGAAATTTGAGAAACAACTAGGCCAAATAATGTTGCTGATACTATGCACATCACAGATCCTCTTGATGCGGACTACATGGGCCTTGTGCGAATCCATCACGTTGGCCACTGGACCTCTGACCACGCTCTGGGAGGGATCTCCAGGAAAATTTTGGAACACCACGATAGCGGTTTCCATGGCAAACATTTTCAGAGGAAGTTATCTAGCAGGAGCAGGTCTGGCCTTCTCATTAATGAAATCTCTAGGAGGAGGTAGGAGAGGCACGGGAGCCCAAGGGGAAACACTGGGAGAGAAATGGAAAAGACGACTGAACCAACTGAGCAAGTCAGAATTTAACACCTATAAAAGGAGTGGGATTATGGAAGTGGACAGATCCGAAGCCAAAGAGGGACTGAAAAGAGGAGAAACAACCAAACATGCAGTGTCGAGAGGAACCGCTAAATTGAGGTGGTTTGTGGAGAGGAACCTTGTGAAACCAGAAGGGAAAGTCATAGACCTCGGTTGTGGAAGAGGTGGCTGGTCATATTATTGCGCTGGGCTGAAGAAAGTCACAGAAGTGAAGGGATATACAAAAGGAGGACCTGGACATGAAGAACCAATCCCAATGGCGACCTATGGATGGAACCTAGTAAAGCTGCATTCCGGGAAAGACGTATTCTTTATACCACCTGAGAAATGTGACACCCTTTTGTGTGATATTGGTGAGTCCTCTCCAAACCCAACTATAGAGGAAGGAAGAACGCTACGCGTCCTAAAGATGGTGGAACCATGGCTCAGAGGAAACCAATTTTGCATAAAAATTCTGAATCCCTACATGCCAAGTGTGGTGGAAACTCTGGAGCAAATGCAAAGAAAACATGGAGGGATGCTAGTGCGAAATCCACTTTCAAGAAATTCCACTCATGAAATGTATTGGGTTTCATGTGGAACAGGAAACATTGTGTCAGCAGTAAACATGACATCCAGAATGTTGCTAAATCGATTCACAATGGCTCACAGGAAACCAACATATGAAAGAGACGTGGACCTAGGCGCCGGAACAAGACACGTGGCAGTGGAACCAGAGGTAGCCAACCTAGATATCATTGGCCAGAGGATAGAGAACATAAAACATGAACACAAGTCAACATGGCATTATGATGAGGACAATCCATACAAAACATGGGCCTATCATGGATCATATGAGGTCAAGCCATCAGGATCAGCCTCATCCATGGTCAATGGCGTGGTGAAACTGCTCACCAAACCATGGGATGTCATCCCCATGGTCACACAA"  
  seq_up <- "TCAATATGCTGAAACGCGCGAGAAACCG"  
  alignment <- pairwiseAlignment(D1, seq_up, type="global-local") 
  print(alignment)
  
  ggplot(GISDDprimer)+
    geom_point(aes(Position_start, Primer))
  
  
  upseq <- "aatttt"
  downseq <- "ggggccc"
  Forward_primer <- c(
    upseq,
    nchar(upseq),
    find_Gc(upseq),
    Tm_GC(upseq,ambiguous=TRUE,variant="Primer3Plus",Na=50,mismatch=TRUE)[["Tm"]],
    Tm_NN(upseq,Na=50)[["Tm"]]
  )
  Reverse_primer <- c(
    downseq,
    nchar(downseq),
    find_Gc(downseq),
    Tm_GC(downseq,ambiguous=TRUE,variant="Primer3Plus",Na=50,mismatch=TRUE)[["Tm"]],
    Tm_NN(downseq,Na=50)[["Tm"]]
  )
  
  

  
  library(rprimer)
  library(Biostrings)

  
  
  library(openPrimeR)
  
  
  library(TmCalculator)
  ntseq <- "TCAATATGCTGAAACGCGCGAGAAACCG"
  Tm_GC(ntseq,ambiguous=TRUE,variant="Primer3Plus",Na=50,mismatch=TRUE)
  Tm_NN(ntseq,Na=50)
  Tm_Wallace(ntseq,ambiguous = TRUE)
  GC(ntseq,ambiguous = TRUE)
  
###############################################################################
#loading templates
  # Specify a FASTA file containing the templates:
  fasta.file <- system.file("extdata", "IMGT_data", "templates", 
                            "Homo_sapiens_IGH_functional_exon.fasta", package = "openPrimeR")
  # Load the template sequences from 'fasta.file'
  seq.df<- read_templates(fasta.file)
  hdr.structure <- c("ACCESSION", "GROUP", "SPECIES", "FUNCTION")
  seq.df <- read_templates(fasta.file, hdr.structure, delim = "|", id.column = "GROUP")
###############################################################################
#Loading and writing settings
  list.files(system.file("extdata", "settings", package = "openPrimeR"), pattern = "*\\.xml")
  settings.xml <- system.file(
    "extdata", 
    "settings", 
    "A_Taq_PCR_design.xml", 
    ?chjpackage = "openPrimeR")
  settings <- read_settings(settings.xml)
###############################################################################
#Individual binding regions
  l.fasta.file <- system.file("extdata", "IMGT_data", "templates", 
                              "Homo_sapiens_IGH_functional_leader.fasta", package = "openPrimeR")
  template.df <- assign_binding_regions(
    seq.df, fw = l.fasta.file, rev = NULL)
  
  template.df <- assign_binding_regions(seq.df, fw = c(1,50), rev = c(1,40))
###############################################################################
#desing
  design.settings <- settings
  constraints(design.settings) <- constraints(design.settings)[!grepl(
    "gc_clamp", names(constraints(design.settings)))]
#desing
  optimal.primers <- design_primers(
    template.df[1:2,], 
    mode.directionality = "fw",
                                    
    settings = design.settings)
###############################################################################
#analyzing primers
  primer.location <- system.file("extdata", "IMGT_data", "primers", "IGHV", 
                                 "Ippolito2012.fasta", package = "openPrimeR")
  # Load the primers
  primer.df <- read_primers(primer.location, fw.id = "_fw")
  ###############################################################################
  
   
  seq.df <- read_templates("Data/1L1_CDS_duplicated_20210124.fasta")
  primer.df <- read_primers("Data/GISDDprimer.fasta")
  
  
  constraint.df <- check_constraints(
    primer.df, 
    template.df, 
    settings, 
    active.constraints = names(constraints(settings)))
  
  
  
  
  constraint.df.site <- check_restriction_sites(
    primer.df,
    template.df,
    adapter.action = c("warn", "rm"),
    selected = NULL,
    only.confident.calls = TRUE,
    updateProgress = NULL
  )
  
  library(rprimer)
  library(Biostrings)
  myAlignment <- readDNAMultipleAlignment("Data/1L1_CDS_duplicated_20210124.fasta", format = "fasta")
  
  myConsensusProfile <- consensusProfile(myAlignment, ambiguityThreshold = 0.05)
  
  plotData(myConsensusProfile)
  
  myOligos <- designOligos(myConsensusProfile)
  
  plotData(myOligos)
  
  myAssays <- designAssays(myOligos)  
  
  
  library(Biostrings)
  primer <- "TCAATATGCTGAAACGCGCGAGAAACCG"
  primer.str <- DNAString(primer)
  sequences <- readDNAStringSet("Data/db_demo.fas","fasta")
  #sequences <- DNAString("AAAAAAAAAATCAATATGCTGAAACGCGCGAGAAACCGCGCTTTGGTGTGCTAGGTGTGAC")
  max.mismatch = 3
  min.mismatch = 0
  matches <- vmatchPattern(primer,
                           sequences,
                           max.mismatch = max.mismatch, 
                           min.mismatch = min.mismatch,fixed = T)
  if (length(matches) == 0) {
    rc_primer <- reverseComplement(primer)
    matches <- vmatchPattern( rc_primer, 
                              sequences, 
                              max.mismatch = max.mismatch, 
                              min.mismatch = min.mismatch,
                              fixed = T)
  }
  return(matches)
  
  