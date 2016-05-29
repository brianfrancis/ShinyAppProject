
shinyUI(
        navbarPage("qPCR Relative Quantification Calculator",
                tabPanel("Instructions",
                      h4('Motivation'),
                      p('The number of copies of an allele in a transgenic animal 
                        can be determined via qPCR assays.  
                        The technique requires comparing the cycle threshold 
                        (Ct values) returned by the assay to that of a sample 
                        where the copies are known (a positive template control or PTC) 
                        to get the relative quantification.'),
                      p('The application assists in doing these calculations and 
                        in grouping the samples in a meaningful way (i.e., by genetic line).  
                        Instrument software frequently has minimal labeling of the attributes.  
                        In the case of the dataset considered here, that includes only the 
                        unique identification for samples which is a concatenation of the 
                        order # under which the samples are being processed and the 
                        sample number for each animal.'),
                      p('However animals of different genetic lines may be tested using 
                        the same assays on the same PCR plate.  Therefore the genetic 
                        line of the samples and which PTC values should be used to 
                        calculate the relative quantification are additional parameters 
                        needed to perform the calculations.'),
                      p('For the configurations below you can enter anything you want
                        for the lines and associated controls.  The calculations will not 
                        be correction but it should be sufficient to evaluate the project.'),
                      h4('Steps to produce relative quantification:'),
                      tags$ol(tags$li("On the configurations page, select each order # and enter a line.
                                        Click update to save the line with the order."), 
                              tags$li("Next on the configurations page, select each line-assay 
                                      to display the availble PTCs to associated.  Check off 
                                      one or more of the PTCs."), 
                              tags$li("To view the final relative quantifications, 
                                      select Output and click Update Calcuations.
                                      This will display the table of results"),
                              tags$li("To see a scatter plot of these results, 
                                      select the Plot sub-tab and update the 
                                      line-target drop down to view a plot for each.")
                      )
                ),
                tabPanel("Configurations",                                             
                         numericInput('maxcycles', 'Maximum Cycles', 40, 
                             min = 35, max = 50, step = 1),
                        fluidRow(column(6,wellPanel(  
                                        h3("Enter Animal Line for each Order"),
                                        uiOutput('selectorder'),
                                        textInput("linename", "Enter Line:", ""),
                                        actionButton("addLineButton", "Update"))),
                                 column(6,tableOutput('orderlist'))
                        ),                          
                        fluidRow(column(6,wellPanel(
                                        h3("Select Controls for each Line Assay"),
                                        uiOutput('selectlineassay'),
                                        uiOutput('checkcontrols'))
                                ),
                                column(6, tableOutput('linetargetlist')
                                )
                        )
                ),            
                tabPanel("PTC Summary",
                        dataTableOutput('ptcSummary')
                        ),
                tabPanel("Output",
                         column(4, offset=8,actionButton("refreshButton", "Update Calculations")),
                         tabsetPanel(type = "tabs", 
                                     tabPanel("Summary", dataTableOutput('output')), 
                                     tabPanel("Plot", uiOutput('selectlineassaygraph'),               
                                              plotOutput('graph'))
                                     )
                        
                         )
                )
)

