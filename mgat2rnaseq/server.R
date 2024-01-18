#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#library(shiny)

# Define server logic required to draw a histogram
function(input, output, session) {

    output$boxplots <- renderPlot({

       plot_boxes(gene = input$gene)

    })
    
    output$table = renderDataTable({
      ann_df
    
    })
    
    output$heatmap = renderPlot({
      generate_heatmap(process_dds_vst(generate_path_gene_list(input$pathway)))
    })
    
    output$gsea_chow <- renderImage({
      list(
        src = "GSEA_chow.png",  # Use the file name of the image
        contentType = "image/png",
        width = 1000
      )
    }, deleteFile = FALSE)
    
    
    output$gsea_lf <- renderImage({
      list(
        src = "GSEA_lf.png",  # Use the file name of the image
        contentType = "image/png",
        width = 1000
      )
    }, deleteFile = FALSE)
    
    output$gsea_hf <- renderImage({
      list(
        src = "GSEA_hf.png",  # Use the file name of the image
        contentType = "image/png",
        width = 1000
      )
    }, deleteFile = FALSE)
    
    

}
