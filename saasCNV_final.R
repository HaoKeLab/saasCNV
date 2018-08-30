library(shiny)
library(shinyjs)
library(saasCNV)
setwd("c:/Users/Jason Ding/Desktop/temp")

vector_choices <- 1:22 #specifies the range of numbers the dropdown menu (that chooses the chromosome to view)can display 

ui <- fluidPage(
  useShinyjs(),
  h1("saas-CNV"), #header
  fluidRow(
    column(4,
           wellPanel(
             h2("Parameters"),
             #wiget that allows user to input file
             fileInput(inputId = "vcf.file", label = "Please input VCF file", multiple = FALSE, accept = NULL, width = NULL,
                       buttonLabel = "Browse...", placeholder = "No file selected"),
             #widget that allows user to specify the minimum number of probes a segment span
             textInput(inputId = "min.snps", label = "the minimum number of probes a segment span (for cluster plot)", value = "10"),
             hr(), #horizontal line (for formatting purposes)
             #displays the widget to choose which chromosome to view after you click 'Run analysis', this will be disabled until the graphs load.
             conditionalPanel(
               condition = "input.clicks >= 1", 
               h5("Please select the chromosome number to view:"),
               div(style="display:inline-block;", actionButton("prevBin", label = "<", width='40px')),
               div(style="display:inline-block;", 
                   selectizeInput("chr.num", "", choices = as.list(vector_choices), selected = 1, width='60px')),
               div(style="display:inline-block;", actionButton("nextBin", label = ">", width='40px'))
             )
           ),
           #runs the analysis on the data
           actionButton(inputId = "clicks", label = "Run analysis"), 
           #creates the 3 download buttons once the 'Run analysis' button is clicked, this will be disabled until the graphs load. 
           conditionalPanel(
             condition = "input.clicks >= 1",
             tags$hr(),
             downloadButton(outputId = "downloadDSP", label = "Download Diagnostic Scatter Plot"),
             downloadButton(outputId = "downloadGWP", label = "Download Genome Wide Plot"),
             tags$hr(),
             downloadButton(outputId = "download.TAR", label = "Download Result File")
           )
    ),
    column(8, 
           h3(textOutput("txt3")),
           plotOutput("diag.seg")
    )),
  fluidRow(  
    column(6,       
           h3(textOutput("txt2")),
           plotOutput("cluster")
    ),
    column(6,
           h3(textOutput("txt1")),
           plotOutput("genome"))
  )
)

server <- function(input, output, session) {
  #disables the download buttons and the widget that chooses which chromosome analysis to display
  disable("downloadGWP")
  disable("downloadDSP")
  disable("download.TAR")
  disable("chr.num")
  
  #code for the widget that chooses which chromosome analysis to display
  observeEvent(input$prevBin, {
    current <- which(vector_choices == input$chr.num)
    if(current > 1){
      updateSelectInput(session, "chr.num", selected = vector_choices[current - 1])
    }
    
  })
  observeEvent(input$nextBin, {
    current <- which(vector_choices == input$chr.num)
    if(current < length(vector_choices)){
      updateSelectInput(session, "chr.num", selected = vector_choices[current + 1])
    }
  })
  
  #runs the whole saasCNV program once the 'Run analysis' button is clicked
  observeEvent(input$clicks, {
    
    #this changes the 'Run analysis' button to say update after one click
    updateActionButton(session, "clicks", label = "Update", icon = icon("refresh")) 
    
    #reads in the vcf file chosen from the widget that inports the vcf file
    vcf_table <- read.delim(input$vcf.file$datapath, as.is=TRUE)
    
#------------------------------------------------------------------------------------------------------------------------    
   #all of these are parameters for the NGS.CNV function. The commented out ones use the default setting
    output.dir <- file.path(getwd(), "test_saasCNV")
    sample.id <- "WES_0116"
    #myNGS.CNV(vcf=vcf_table, output.dir=output.dir, sample.id=sample.id,
    #min.chr.probe=100,
    #min.snps=10,
    #joint.segmentation.pvalue.cutoff=1e-4,
    #max.chpts=30,
    #do.merge=TRUE, use.null.data=TRUE, num.perm=1000, maxL=2000,
    #merge.pvalue.cutoff=0.05,
    #do.cnvcall.on.merge=TRUE,
    #cnvcall.pvalue.cutoff=0.05,
    #do.plot=TRUE, cex=0.3, ref.num.probe=1000,
    #do.gene.anno=TRUE,
    #gene.anno.file="refGene_hg19.txt.gz",
    #seed=123456789,
    #verbose=TRUE) 
#------------------------------------------------------------------------------------------------------------------    
    #Start of NGS.CNV
    
    #making progression bar
    withProgress(message = "", value = 0, {
      vcf = vcf_table
      do.GC.adjust = FALSE
      gc.file = system.file("extdata", "GC_1kb_hg19.txt.gz", package = "saasCNV")
      min.chr.probe = 100
      min.snps = as.numeric(input$min.snps)
      joint.segmentation.pvalue.cutoff = 1e-04
      max.chpts = 30
      do.merge = TRUE
      use.null.data = TRUE
      num.perm = 1000
      maxL = NULL
      merge.pvalue.cutoff = 0.05
      do.cnvcall.on.merge = TRUE
      cnvcall.pvalue.cutoff = 0.05
      do.plot = TRUE
      cex = 0.3
      ref.num.probe = NULL
      do.gene.anno = FALSE
      gene.anno.file = NULL 
      seed = NULL
      verbose = TRUE 
      
      incProgress(1/11, detail = "Doing data input ...") #updates progression bar
      if (!file.exists(output.dir)) 
        dir.create(output.dir)
      seq.data <- cnv.data(vcf = vcf, min.chr.probe = min.chr.probe, 
                           verbose = verbose)
      incProgress(1/11, detail = "Doing GC adjustment ...") #updates progression bar
      if (do.GC.adjust) {
        gc <- read.delim(file = gc.file, as.is = TRUE)
        seq.data <- GC.adjust(data = seq.data, gc = gc, maxNumDataPoints = 10000)
      }
      data.dir <- paste0(output.dir, "/mid_data")
      if (!file.exists(data.dir)) 
        dir.create(data.dir)
      save(list = "seq.data", file = file.path(data.dir, "seq.data.RData"))
      
      incProgress(1/11, detail = "Doing joint segmentation ...") #updates progression bar
      seq.segs <- joint.segmentation(data = seq.data, min.snps = min.snps, 
                                     global.pval.cutoff = joint.segmentation.pvalue.cutoff, 
                                     max.chpts = max.chpts, verbose = verbose)
      write.table(seq.segs, file = file.path(data.dir, "seq.segs.txt"), 
                  quote = FALSE, row.names = FALSE, sep = "\t")
      incProgress(1/11, detail = "Doing segment plots ...") #updates progression bar
      if (do.plot) {
        plot.dir <- paste0(output.dir, "/mid_seg_plot")
        if (!file.exists(plot.dir)) 
          dir.create(plot.dir)
        chrs <- sub("^chr", "", unique(seq.segs$chr))
        for (chr in chrs) {
          incProgress(1/242, detail = paste("Doing diagnostic segment plot for chromosome ", chr)) #updates progression bar
          png(filename = file.path(plot.dir, paste0("seg_chr", 
                                                    chr, ".png")), width = 240 * 5, height = 240 * 
                4)
          diagnosis.seg.plot.chr(data = seq.data, segs = seq.segs, 
                                 sample.id = sample.id, chr = chr, cex = cex)
          dev.off()
        }
      }
      incProgress(1/11, detail = "Doing diagnosis plot ...") #updates progression bar
      if (do.plot) {
        diagnosis.plot.dir <- file.path(output.dir, "mid_diagnosis_plot")
        if (!file.exists(diagnosis.plot.dir)) 
          dir.create(diagnosis.plot.dir)
        chrs <- sub("^chr", "", unique(seq.segs$chr))
        png(filename = file.path(diagnosis.plot.dir, "diag_QQ_plot.png"), 
            width = 240 * 3, height = 240 * 2)
        diagnosis.QQ.plot(data = seq.data, chrs = chrs)
        dev.off()
      }   
      incProgress(1/11, detail = "Merging segments ...") #updates progression bar
      if (do.merge) {
        seq.segs.merge <- merging.segments(data = seq.data, segs.stat = seq.segs, 
                                           use.null.data = use.null.data, N = num.perm, maxL = maxL, 
                                           merge.pvalue.cutoff = merge.pvalue.cutoff, seed = seed, 
                                           verbose = verbose)
        write.table(seq.segs.merge, file = file.path(data.dir, 
                                                     "seq.segs.merge.txt"), quote = FALSE, row.names = FALSE, 
                    sep = "\t")
        if (do.plot) {
          plot.dir <- paste0(output.dir, "/mid_seg_merge_plot")
          if (!file.exists(plot.dir)) 
            dir.create(plot.dir)
          chrs <- sub("^chr", "", unique(seq.segs.merge$chr))
          for (chr in chrs) {
            incProgress(1/242, detail = paste("Doing diagnosis segment plot after merge for chr", chr)) #updates progression bar
            png(filename = file.path(plot.dir, paste0("merge_seg_chr", 
                                                      chr, ".png")), width = 240 * 5, height = 240 * 
                  4)
            diagnosis.seg.plot.chr(data = seq.data, segs = seq.segs.merge, 
                                   sample.id = sample.id, chr = chr, cex = cex)
            dev.off()
          }
        }
      }
      incProgress(1/11, detail = "Doing CNV calling ...") #updates progression bar
      cat("CNV calling ...\n")
      if (do.merge == TRUE & do.cnvcall.on.merge == TRUE) {
        seq.cnv <- cnv.call(data = seq.data, sample.id = sample.id, 
                            segs.stat = seq.segs.merge, maxL = maxL, N = num.perm, 
                            pvalue.cutoff = cnvcall.pvalue.cutoff, seed = seed)
      }
      else {
        seq.cnv <- cnv.call(data = seq.data, sample.id = sample.id, 
                            segs.stat = seq.segs, maxL = maxL, N = num.perm, 
                            pvalue.cutoff = cnvcall.pvalue.cutoff, seed = seed)
      }
      res.dir <- paste0(output.dir, "/mid_res")
      incProgress(2/33, detail = "Writing tables ...") #updates progression bar
      if (!file.exists(res.dir)) 
        dir.create(res.dir)
      write.table(seq.cnv, file = file.path(res.dir, "seq.cnv.txt"), 
                  quote = FALSE, row.names = FALSE, sep = "\t")
      incProgress(2/33, detail = "Doing gene annotations ...") #updates progression bar
      if (do.gene.anno) {
        cat("gene annotation ...\n")
        gene <- read.delim(file = gene.anno.file, as.is = TRUE, 
                           comment.char = "")
        seq.cnv.anno <- reannotate.CNV.res(res = seq.cnv, gene = gene, 
                                           only.CNV = TRUE)
        write.table(seq.cnv.anno, file = file.path(res.dir, "seq.cnv.anno.txt"), 
                    quote = FALSE, row.names = FALSE, sep = "\t")
      }
      incProgress(2/33, detail = "Finishing up ...") #updates progression bar
      png(filename = file.path(res.dir, "cnv_gw_plot.png"), width = 240 * 
            6, height = 240 * 4)
      genome.wide.plot(data = seq.data, segs = seq.cnv, sample.id = sample.id, 
                       chrs = sub("^chr", "", unique(seq.cnv$chr)), cex = cex)
      dev.off()
      pdf(file = file.path(res.dir, "cnv_cluster_plot.pdf"), width = 8, 
          height = 8)
      diagnosis.cluster.plot(segs = seq.cnv, chrs = sub("^chr", 
                                                        "", unique(seq.cnv$chr)), min.snps = min.snps, max.cex = 3, 
                             ref.num.probe = ref.num.probe)
      dev.off() #end of NGS.CNV
    
    #displays graphs and their headers  
    output$txt1 <- renderText("Genome Wide Plot")
    gwp <- function(){genome.wide.plot(seq.data, seq.cnv, "MY SAMPLE" , 1:22, cex = 0.3)}
    output$genome <- renderPlot(gwp())
    output$txt2 <- renderText("Diagnosis Scatter Plot")
    dsp <- function(){diagnosis.cluster.plot(seq.cnv, 1:22, input$min.snps, max.cex = 3, ref.num.probe = NULL)}
    output$cluster <- renderPlot(dsp())
    output$txt3 <- renderText(paste("Diagnosis Segment Plot: Chr", input$chr.num, sep = ""))
    output$diag.seg <- renderPlot(diagnosis.seg.plot.chr(data = seq.data, segs = seq.segs, 
                                                         sample.id = sample.id, chr = input$chr.num, cex = cex))
    #links the download button to the picture
    output$downloadGWP <- downloadHandler(
      filename =  "mygwp.png", 
      # content is a function with argument file. content writes the plot to the device
      content = function(file) {
        png(file) # open the png device
        gwp()
        dev.off()  # turn the device off
      }, 
      contentType = "image/png"
    )
    
    #links the download button to the picture
    output$downloadDSP <- downloadHandler(
      filename <- function() {paste("myDSP", "png", sep=".")},
      content = function(file) {
        png(file)
        dsp()
        dev.off()
      },
      contentType = "image/png"
    )
    
    #links the download button to the tar file
    output$download.TAR <- downloadHandler(
      filename = function() {
        "saas-CNV.tar"
      },
      content = function(file) {
        tar(file, "test_saasCNV")
      }
    )
    }) #end of progression bar
    
   #enables the download buttons and the widget.
    enable("downloadGWP")
    enable("downloadDSP")
    enable("download.TAR")
    enable("chr.num")
  })
}
shinyApp(ui = ui, server = server)