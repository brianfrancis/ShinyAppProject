##problem with linking delta ptc back to line target list (going to every assay??)

library(plyr)
library(ggplot2)

rawdata <- read.csv("Data.csv")

processeddata <- rawdata[,c('Well','Well.Position', 'Sample.Name', 'Target.Name',
                                'CT')]

##change CT from factor to character
processeddata$CT <- as.character(processeddata$CT)

##mark the assay control wells
processeddata$isControlAssay <- grepl('APOB', toupper(processeddata$Target.Name))

##mark the assay control wells
processeddata$isSampleControl <- grepl('PTC', toupper(processeddata$Sample.Name))

x <- as.character(processeddata[processeddata$isSampleControl==FALSE, 
                                c('Sample.Name')])

x <- strsplit(x,"-")
x <- data.frame(t(sapply(x, `[`)))

processeddata$OrderNbr[processeddata$isSampleControl==FALSE] <- as.character(x[,1])
processeddata$SampleIdent[processeddata$isSampleControl==FALSE] <- x[,2]

orderlist <- as.data.frame(unique(processeddata$OrderNbr[!is.na(processeddata$OrderNbr)]))
names(orderlist) <- c("OrderNbr")
orderlist$Line <- NA

ordertargetlist <- as.data.frame(unique(processeddata[processeddata$isSampleControl==FALSE  &
                                                              processeddata$isControlAssay==FALSE,
                                                     c("OrderNbr", "Target.Name")]))
names(ordertargetlist) <- c("OrderNbr", "Target")

controllist <- as.data.frame(unique(processeddata[processeddata$isSampleControl==TRUE,
                                                  c("Sample.Name", "Target.Name")]))
controllist
names(controllist) <- c("control", "target")

##replace any undetermined values with the user maximum number of cycles
##also convert to numeric
replaceUND <- function(maxcycles, data) {
        x <- data
        x[x$CT=='Undetermined',]$CT <- maxcycles
        x$CT <- as.numeric(x$CT)
        x
}

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

getMedianPTC <- function(data){
        ptcdata <- data[data$isSampleControl==TRUE,]
        x <- aggregate(deltaCT ~ Sample.Name + Target, data = ptcdata, mean)
        names(x) <- c("Control", "Target", "Median Delta CT")
        x
        
}

##plot the delta delta by sample
plotDeltaDelta <- function(data, linetarget){
        x <- data[data$linetarget==linetarget,]
        g <- ggplot(x, aes(x = paste(OrderNbr, SampleIdent, sep='-'), y = DeltaDeltaCT)) 
        g <- g + geom_point() 
        g <- g + labs(title = paste("Relative Quantification of qPCR: ", linetarget),x = "Sample", y="Delta Delta CT")
        g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
        g
}


##trim leading and trailing space on a string
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

shinyServer(
        function(input, output) {
                values <- reactiveValues(df_orderlist = orderlist
                                         ,df_linetargetlist = NULL
                                         ,df_controllist = controllist
                                         ,df_output = NULL
                                         )
                
                alldata <- reactive({
                                        x <- replaceUND(input$maxcycles, processeddata)
                                        x <- getdeltaCT(x)
                                        x
                                })      
                
                ptcSummary <- reactive({
                                        getMedianPTC(alldata())                        
                                })
                                
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
                       ##u$deltaCT <- NA
                        values$df_linetargetlist <- u
                })  
                
                selectedtarget <- reactive({x <- input$lineassay
                                            trim(unlist(strsplit(x,"[|]"))[2])
                                            })
                
                observeEvent(input$controls, {
                                x <- values$df_linetargetlist
                                
                                ##store the selected controls with each line-target
                                x[[x$linetarget==input$lineassay,4]] <- list(input$controls)
                                x[x$linetarget==input$lineassay,"Controls"] <- paste(input$controls, sep='', collapse='; ')
                                
                                ##if multiple ptcs selected take the mean of the replicate medians
                                ##store with the selected line-target
                                y <- ptcSummary()
                                y <- y[y$Target==selectedtarget(),]
                                ptc <- y[y$Control %in% input$controls, 3]
                                x[x$linetarget==input$lineassay,"PTCdeltaCT"] <- mean(ptc)
                                
                                values$df_linetargetlist <- x                                                        
                                }
                             )
                
                ##update the control 
                observeEvent(alldata(), {
                        ##recalcualte the median controls and store with each line target (for loop)
                        
                })
                
                observeEvent(input$refreshButton, {
                                           
                        m <- merge(values$df_linetargetlist[!is.na(values$df_linetargetlist$Line),], ordertargetlist)
                        m <- m[,c("OrderNbr", "Line", "Target", "PTCdeltaCT", "linetarget")]
                        m <- merge(alldata(), m)
                        m$DeltaDeltaCT <- m$deltaCT - m$PTCdeltaCT
                        m <- m[,c("Well", "Well.Position", "Line", "Target","OrderNbr", "SampleIdent", "DeltaDeltaCT", "linetarget")]
                        values$df_output <- m
                        output$output <- renderDataTable({values$df_output})
                        
                })
                
                                
                output$orderlist <- renderTable({values$df_orderlist})
                output$selectorder <- renderUI({
                        selectInput('ordernbr', 'Select Order', as.character(values$df_orderlist[,1]))
                                    })
                        
        
                output$linetargetlist <- renderTable({values$df_linetargetlist[,-c(3,4)]})
               
                
                
                output$selectlineassay <- renderUI({
                        selectInput('lineassay', 'Select Line-Assay', values$df_linetargetlist$linetarget)
                })
                
                
                ##change this to be line assay and then allow them to pick control
                observeEvent(input$lineassay, {
                                output$checkcontrols <- renderUI({
                                        checkboxGroupInput('controls', 'Select controls for line-assay:',
                                                           as.character(values$df_controllist[
                                                                   values$df_controllist$target==selectedtarget(),1]))
                                })
                })
                
                output$alldata <- renderDataTable({alldata()})
                output$ptcSummary <- renderDataTable({ptcSummary()})
                
                output$selectlineassaygraph <- renderUI({
                        selectInput('lineassaygraph', 'Select Line-Assay', values$df_linetargetlist$linetarget)
                })
                
                
                observeEvent(input$lineassaygraph, {
                        output$graph <- renderPlot({plotDeltaDelta(values$df_output, input$lineassaygraph)})
                })
                
                
        }
)


