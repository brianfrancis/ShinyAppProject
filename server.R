
###Libraries

library(plyr)
library(ggplot2)

###Read and process data
        ##Notes:
        ## data are multiplexed so target and endogenous control are in the same wells
        ## endogenous controls are rAPOB or mAPOB
        ## sample controls must have "PTC" in their Sample.Name
        ## unknown samples whould be entered as [order]-[sample identifier]
        

##read raw data
rawdata <- read.csv("Data.csv")

##subset data to columns we care about
processeddata <- rawdata[,c('Well','Well.Position', 'Sample.Name', 'Target.Name',
                                'CT')]

##change CT from factor to character
processeddata$CT <- as.character(processeddata$CT)

##mark the assay control wells (endogenous controls)
processeddata$isControlAssay <- grepl('APOB', toupper(processeddata$Target.Name))

##mark the assay control wells
processeddata$isSampleControl <- grepl('PTC', toupper(processeddata$Sample.Name))

##parse Sample.Name field in order nbr and sample identifier for unknown samples
x <- as.character(processeddata[processeddata$isSampleControl==FALSE, 
                                c('Sample.Name')])
x <- strsplit(x,"-")
x <- data.frame(t(sapply(x, `[`)))
processeddata$OrderNbr[processeddata$isSampleControl==FALSE] <- as.character(x[,1])
processeddata$SampleIdent[processeddata$isSampleControl==FALSE] <- x[,2]

## get a unique list of the orders in on the plate
orderlist <- as.data.frame(unique(processeddata$OrderNbr[!is.na(processeddata$OrderNbr)]))
names(orderlist) <- c("OrderNbr")
orderlist$Line <- NA

## get a unique list of the orders and gene targeting asssays on the plate
ordertargetlist <- as.data.frame(unique(processeddata[processeddata$isSampleControl==FALSE  &
                                                              processeddata$isControlAssay==FALSE,
                                                     c("OrderNbr", "Target.Name")]))
names(ordertargetlist) <- c("OrderNbr", "Target")

## get a unique list of the PTCs on the plate
controllist <- as.data.frame(unique(processeddata[processeddata$isSampleControl==TRUE,
                                                  c("Sample.Name", "Target.Name")]))
controllist
names(controllist) <- c("control", "target")

##function to replace any undetermined values with the user maximum number of cycles
##also convert to numeric
replaceUND <- function(maxcycles, data) {
        x <- data
        x[x$CT=='Undetermined',]$CT <- maxcycles
        x$CT <- as.numeric(x$CT)
        x
}

##function to calculate delta CT (endo results minus gene target result)
##subtract the control assay from the target assay to get the delta CT
##return a data frame with one row per well
getdeltaCT <- function(data) {
        target <- data[data$isControlAssay==FALSE,]
        target <- subset(target, select = -c(isControlAssay))
        target <- rename(target, c("CT"="targetCT"))
        
        control <- data[data$isControlAssay==TRUE,c('Well', 'CT')]
        
        control <- rename(control, c("CT"="controlCT"))
        
        x <- merge(target,control)
        x$deltaCT <- x$controlCT - x$targetCT
        x <- rename(x, c("Target.Name"="Target"))
        x
}

##function to get the median PTC for the replicates for an individual PTC
getMedianPTC <- function(data){
        ptcdata <- data[data$isSampleControl==TRUE,]
        x <- aggregate(deltaCT ~ Sample.Name + Target, data = ptcdata, mean)
        names(x) <- c("Control", "Target", "Median Delta CT")
        x
        
}

##function to plot the delta delta by sample
plotDeltaDelta <- function(data, linetarget){
        x <- data[data$linetarget==linetarget,]
        g <- ggplot(x, aes(x = paste(OrderNbr, SampleIdent, sep='-'), y = DeltaDeltaCT)) 
        g <- g + geom_point() 
        g <- g + labs(title = paste("Relative Quantification of qPCR: ", linetarget),x = "Sample", y="Delta Delta CT")
        g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
        g
}


## function to trim leading and trailing space on a string
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

shinyServer(
        function(input, output) {
                
                ##  Reactive values ##########
                
                ##list of reactive values
                values <- reactiveValues(df_orderlist = orderlist
                                         ##store line targets (populates drop down)
                                         ,df_linetargetlist = NULL
                                         ##store line targets and their associated controls
                                         ,df_linetargetcontrolist = NULL                                         
                                         ,df_finaloutput = NULL
                                         )
                
                ## all data has any UND values replaced basesd on user defined value
                alldata <- reactive({
                                        x <- replaceUND(input$maxcycles, processeddata)
                                        x <- getdeltaCT(x)
                                        x
                                })      
                
                ## get a list median PTCs
                ptcSummary <- reactive({
                                        getMedianPTC(alldata())                        
                                })
                
                ##parse out the target selected when user selects a line-target
                selectedtarget <- reactive({x <- input$lineassay
                                            trim(unlist(strsplit(x,"[|]"))[2])
                })
                
                ##  Observe events ##########
        
                ## action when user add a line
                        ##append line to order list
                        ##recalculate unique list of line targets
                observeEvent(input$addLineButton, {
                        x <- values$df_orderlist
                        x$Line[x$Order==input$ordernbr] <- input$linename
                        values$df_orderlist <- x
                        
                        y <- ordertargetlist
                        m <- merge(x,y)
                        u <- unique(m[,2:3])
                       ## u <- u[order(u),]
                        u$linetarget <- paste(u[,1],'|',u[,2])
                       u <- cbind(u,list(NA))
                       names(u)[4] <- "controllist"
                       u$Controls <- NA
                
                       ##store target list for later lookup
                       values$df_linetargetlist <- u
                       
                       ##merge target list with existing target control list (if it's been used)
                       v <- data.frame(values$df_linetargetcontrollist)
                       if (nrow(v)==0) {
                               print ('No rows')
                               values$df_linetargetcontrollist <- u
                       } else {
                                print('some rows')
                               values$df_linetargetcontrollist <- 
                                       merge(v, u[,1:3], all.y=TRUE)
                       }
                })                  
                
                ##when the selected line target changes, update the list of controls that can be selected
                observeEvent(input$lineassay, {
                        output$checkcontrols <- renderUI({
                                checkboxGroupInput('controls', 'Select controls for line-assay:',
                                                   as.character(controllist[
                                                           controllist$target==selectedtarget(),1]))
                        })
                })                
                
                ##action when user selects a control for a line-target
                        ##store control with line-target list
                        ##recalculate mean ptc for the line-target (average of the
                        ## median of the replicate ptcs calculated earlier)
                observeEvent(input$controls, {
                                x <- values$df_linetargetcontrollist
                                
                                ##store the selected controls with each line-target
                                x[[x$linetarget==input$lineassay,4]] <- list(input$controls)
                                x[x$linetarget==input$lineassay,"Controls"] <- paste(input$controls, sep='', collapse='; ')
                                
                                ##if multiple ptcs selected take the mean of the replicate medians
                                ##store with the selected line-target
                                y <- ptcSummary()
                                y <- y[y$Target==selectedtarget(),]
                                ptc <- y[y$Control %in% input$controls, 3]
                                x[x$linetarget==input$lineassay,"PTCdeltaCT"] <- mean(ptc)
                                
                                values$df_linetargetcontrollist <- x                                                        
                                }
                             )
                
                ## action when user refreshes output
                        ##recalculate delta delta
                observeEvent(input$refreshButton, {
                        m <- values$df_linetargetcontrollist[!is.na(values$df_linetargetcontrollist$Line),]
                        m <- merge(m, ordertargetlist)
                        m <- merge(m, values$df_orderlist)
                        m <- m[,c("OrderNbr", "Line", "Target", "PTCdeltaCT", "linetarget")]
                        m <- merge(alldata(), m)
                        m$DeltaDeltaCT <- m$deltaCT - m$PTCdeltaCT
                        m <- m[,c("Well", "Well.Position", "Line", "Target","OrderNbr", "SampleIdent", "DeltaDeltaCT", "linetarget")]
                        m <- m[order(m$Well),]
                        values$df_finaloutput <- m                        
                        output$output <- renderDataTable({values$df_finaloutput})
                        
                })
                
                ##when the selected line-target changes for the plot, update the plot
                observeEvent(input$lineassaygraph, {
                        output$graph <- renderPlot({plotDeltaDelta(values$df_finaloutput, input$lineassaygraph)})
                })
                
                ## Return output for display in UI
                
                ##return table to display order list                
                output$orderlist <- renderTable({values$df_orderlist})
                
                ##create drop down tot select and order nbr to enter the appropriate line
                output$selectorder <- renderUI({
                        selectInput('ordernbr', 'Select Order', as.character(orderlist[,1]))
                                    })
                        
                ##return table to display the list of lines and targets
                output$linetargetlist <- renderTable({values$df_linetargetcontrollist[,-c(3,4)]})
               
                
                ##populate line-target drop down to allow associating PTCs
                output$selectlineassay <- renderUI({
                        selectInput('lineassay', 'Select Line-Assay', values$df_linetargetlist$linetarget)
                })
                
                                
                ##return PTC summary for display
                output$ptcSummary <- renderDataTable({ptcSummary()})
                
                ##populate drop down for selecting line-target for the plot
                output$selectlineassaygraph <- renderUI({
                        selectInput('lineassaygraph', 'Select Line-Assay', values$df_linetargetlist$linetarget)
                })
                
                
        }
)


