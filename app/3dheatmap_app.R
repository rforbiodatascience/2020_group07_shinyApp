# https://softeng.oicr.on.ca/hardeep_nahal/2017/04/25/RBioc/

library(shiny)
library(d3heatmap)
library(RColorBrewer)
library(shinythemes)
library(pheatmap)
library(tidyverse)

# LOAD DATA
# -----------------------------
proteome_data_wide <- read.csv("/home/rebekabato/DTU_Studies/4-Semester/R_for_Bioscience/Exam_Project/02_proteome_data_wide_aug.csv")


# PREPARE THE DATA FOR HEATMAP
# -----------------------------
# sort by patient_ID
# proteome_data_wide <- proteome_data_wide[order(proteome_data_wide$patient_ID)]

# set row names attribute to patient_ID
row.names(proteome_data_wide) <- proteome_data_wide$patient_ID

# and deselect patient_ID, by now it is duplicated
proteome_data <- proteome_data_wide %>% 
  select(-patient_ID) %>% 
  data.matrix()

# save refseq accession numbers as array - > possible choices for selectInput
refSeq_acc_numbers = colnames(proteome_data_wide)[0:-1]


# SHINY APP
# ---------------------------

# USER INTERFACE ------------
# we use two wrapper functions to implement d3heatmap in shiny ->
#     d3heatmapOutput - > create a UI element whenever the app runs
#     renderD3heatmap - > render the actual heatmap (used later in the server seccion)
ui <- fluidPage(
  
  # Add title
  titlePanel("Expression of different protein isoforms in breast cancer data."),
  
  # Give an interesting theme to your app
  theme = shinytheme("united"),
  
  # Sidebar layout, elements: side panel + main panel
  sidebarLayout(
    
    # Sidebar panel for inputs
    sidebarPanel(
      
      # Input: check if clustering is needed or not
      checkboxInput(inputId = "cluster", label = "Apply clustering")
    ),
    
    # Main panel for displaying outputs
    mainPanel(
      
      # Output: d3heatmap
      d3heatmapOutput(outputId = "heatmap", width = "100%", height = "700px")
    )
  )
)


# SERVER -------------------
server <- function(input, output) {

  # d3heatmap of the breast cancer data
  # with clustering option
  
  output$heatmap <- renderD3heatmap({d3heatmap(proteome_data,
                                               colors = "RdYlBu",
                                               dendrogram = if(input$cluster) "both" else "none")})
}

# CALL SHINY APP -----------
shinyApp(ui = ui, server = server)