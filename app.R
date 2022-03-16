library(shiny)
library(ggplot2)
library(readr)

# Data pre-processing ----
## RNA-seq counts
df <- read.table(file = "data/all_colocalisation_shiny_df_hum.txt", header=T, stringsAsFactors = F, sep="\t")

#here adding symbols to ENSG. Just58 genes see ensg_mapping.R


assayNames <- sort(unique(df$assay))
names(assayNames) <- assayNames
traitNames <- sort(unique(df$trait))
names(traitNames) <- traitNames
geneNames <- sort(unique(df$gene))
names(geneNames) <- geneNames

colnames(df)[colnames(df)=="gene"|colnames(df)=="gene.peak"] <- "gene or peak"

df$color = as.character(cut(df$R2,breaks=c(0,0.2,0.4,0.6,0.8,1),
                            labels=c('blue4','skyblue','darkgreen','orange','red'),
                            include.lowest=TRUE))
df$color <- ifelse(df$snp == df$gwas_lead,'purple',df$color)

df$shape = ifelse(df$snp == df$gwas_lead, 23,
                  ifelse(df$is_top == 1, 22, 21))

df$size = ifelse(df$snp == df$gwas_lead, 3,
                 ifelse(df$is_top == 1, 3, 2))

# Define UI ----
ui <- fluidPage(

  titlePanel("Colocalization between Treg QTLs and immune GWAS"),

  helpText("Visualise the p-values for QTL and GWAS studies across
            all SNPs in a locus"),

  navbarPage(

    "Treg coloc",

    tabPanel("Analysis", {

      sidebarLayout(position = "left",
                    sidebarPanel(

                      selectInput("assayNames",
                                  p("Select an assay"),
                                  choices = assayNames,
                                  selected = assayNames[1]),
                      selectInput("traitNames",
                                        p("Select a trait"),
                                        choices = traitNames,
                                        selected = traitNames[1]),
                      selectInput("geneNames",
                                  p("Select a gene or peak"),
                                  choices = geneNames,
                                  selected = geneNames[1])

                    ),

                    mainPanel(
                      h4("Colocalisation plot",
                         plotOutput("colocPlot"),
                         align="center"),
                      h4("Coloc associated values",
                         tableOutput("colocTable"),
                         align="center")
                    )
      )

    },

    tabPanel("About", {mainPanel(
      p("The work presented here is described in ",
        a("this paper.",
          href = "https://www.biorxiv.org/content/10.1101/654632v1"),
      "If you plan on using any plots or information please cite this work."),
      p("All LD scores are calculated based on the EUR population from 1000G phase 3."),
      p("The GWAS summary stats used are: allergic diseases (",
        a("ALL",
          href = "https://pubmed.ncbi.nlm.nih.gov/29083406/"),
        "), ankylosing spondylitis (",
        a("AS",
          href = "https://pubmed.ncbi.nlm.nih.gov/23749187/"),
        "), asthma (",
        a("AST",
          href = "https://pubmed.ncbi.nlm.nih.gov/29273806/"),
        "), celiac disease (",
        a("CEL",
          href = "https://pubmed.ncbi.nlm.nih.gov/22057235/"),
        "), multiple sclerosis (",
        a("MS",
          href = "https://pubmed.ncbi.nlm.nih.gov/24076602/"),
        "), primary biliary cirrhosis (",
        a("PBC",
          href = "https://pubmed.ncbi.nlm.nih.gov/26394269/"),
        "), psoriasis (",
        a("PS",
          href = "https://pubmed.ncbi.nlm.nih.gov/23143594/"),
        "), rheumatoid arthritis (",
        a("RA",
          href = "https://pubmed.ncbi.nlm.nih.gov/24390342/"),
        "), systemic lupus erythematosus (",
        a("SLE",
          href = "https://pubmed.ncbi.nlm.nih.gov/26502338/"),
        "), type 1 diabetes (",
        a("T1D",
          href = "https://pubmed.ncbi.nlm.nih.gov/25751624/"),
        "), vitiligo (",
        a("VIT",
          href = "https://pubmed.ncbi.nlm.nih.gov/27723757/"),
        "), inflammatory bowel disease, Crohn’s disease and ulcerative colitis (",
        a("IBD, CD and UC",
          href = "https://pubmed.ncbi.nlm.nih.gov/28067908/"),
        ")."
                )
        )
    })
    )
  )
  )


# Define server logic ----
server <- function(session, input, output) {

  # Update based on the assay change event
  observeEvent(
    input$assayNames,
    updateSelectInput(session, "traitNames", "Select a trait",
                      choices = unique(df$trait[df$assay==input$assayNames])))

  observeEvent(
    input$traitNames,
    updateSelectInput(session, "geneNames", "Select a gene or peak",
                      choices = unique(df$gene[df$assay==input$assayNames & df$trait==input$traitNames])))

  output$colocPlot <- renderPlot({
    p1 <-
    ggplot(df[df$assay==input$assayNames & df$trait==input$traitNames & df$gene==input$geneNames,],
           aes(x=-log10(pvalue), y=-log10(nom_pvalue))) +
      geom_point(aes(size=size,shape=shape, fill=color),alpha=0.8)+
      scale_fill_identity() +
      scale_size_identity() +
      scale_shape_identity() +
      ylab(bquote(-log10(P[GWAS]))) +
      xlab(bquote(-log10(P[QTL]))) +
      theme_classic(base_size=14)

    legend_box = data.frame(x = 0.8, y = seq(0.4, 0.2, -0.05))

    p1 = cowplot::ggdraw(p1) +
        geom_rect(data = legend_box,
                  aes(xmin = x, xmax = x + 0.05, ymin = y, ymax = y + 0.05),
                  color = "black",
                  fill = rev(c("blue4", "skyblue", "darkgreen", "orange", "red"))) +
      cowplot::draw_label("0.8", x = legend_box$x[1] + 0.05, y = legend_box$y[1], hjust = -0.3, size = 12) +
      cowplot::draw_label("0.6", x = legend_box$x[2] + 0.05, y = legend_box$y[2], hjust = -0.3, size = 12) +
      cowplot::draw_label("0.4", x = legend_box$x[3] + 0.05, y = legend_box$y[3], hjust = -0.3, size = 12) +
      cowplot::draw_label("0.2", x = legend_box$x[4] + 0.05, y = legend_box$y[4], hjust = -0.3, size = 12) +
      cowplot::draw_label(parse(text = "r^2"), x = legend_box$x[1] + 0.05, y = legend_box$y[1], vjust = -2, size = 14)


    p2 <-
      ggplot(df[df$assay==input$assayNames & df$trait==input$traitNames & df$gene==input$geneNames,],
             aes(y=-log10(pvalue), x=start_snp)) +
      geom_point(aes(size=size,shape=shape, fill=color),alpha=0.8)+
      scale_fill_identity() +
      scale_size_identity() +
      scale_shape_identity() +
      ylab(bquote(-log10(P[GWAS]))) +
      scale_x_continuous(labels=function(x){sprintf('%.3f',x/1e6)})+
      xlab(paste0('chr',unique(df[df$assay==input$assayNames & df$trait==input$traitNames & df$gene==input$geneNames,]$chromosome_snp),' (Mb)'))+
      theme_classic(base_size=14)

    p3 <-
      ggplot(df[df$assay==input$assayNames & df$trait==input$traitNames & df$gene==input$geneNames,],
             aes(y=-log10(nom_pvalue), x=start_snp)) +
      geom_point(aes(size=size,shape=shape, fill=color),alpha=0.8)+
      scale_fill_identity() +
      scale_size_identity() +
      scale_shape_identity() +
      ylab(bquote(-log10(P[QTL]))) +
      scale_x_continuous(labels=function(x){sprintf('%.3f',x/1e6)})+
      xlab(paste0('chr',unique(df[df$assay==input$assayNames & df$trait==input$traitNames & df$gene==input$geneNames,]$chromosome_snp),' (Mb)'))+
      theme_classic(base_size=14)

    p2 = p2 + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
    p4 = cowplot::plot_grid(p2, p3, align = "v", nrow = 2, rel_heights=c(0.8,1))
    p5 = cowplot::plot_grid(p1, p4)
    return(p5)

  })

  output$colocTable <- renderTable(unique(df[df$assay==input$assayNames & df$trait==input$traitNames & df$gene==input$geneNames,c(1,2,3,20:24)])
                                   , digits = 4)
}


# Run the app ----
shinyApp(ui = ui, server = server)
