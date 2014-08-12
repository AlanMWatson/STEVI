


loadImages <- function(fullMonty = FALSE, count=FALSE, outline=FALSE, colocalization=FALSE, overlay=FALSE, ...){
          
          require("EBImage")
          
          imagedir <<- list.dirs(full.names=FALSE, recursive=FALSE)
          #imagedir <<- list.files()
          imagefiles <<- list.files(pattern=".tif", recursive=TRUE)
          dirnumber <<- length(imagedir)
          
          #         if(!file.exists("tmp")){
          #                    dir.create("tmp")
          #          }
          
          if(count == TRUE | fullMonty == TRUE){
                    
                    countnucleus(...)    
          }
          
          if(outline == TRUE | fullMonty == TRUE){
                    
                    outlinecells(...)
          }
          
          if(colocalization == TRUE | fullMonty == TRUE){
                    
                    colocalization(...)
          }
          
          if(overlay == TRUE | fullMonty == TRUE){
                    
                    overlay(...)
          }
          
          statQuant(...)
          
          exportIndCells(...)
          
          organizeForPlot(...)
          
}         


countnucleus <- function(threshold = 0.2, brush = 5, hipass = FALSE, writeimage=FALSE, writereport=FALSE, ...){
          
          ## Find directory that has images of nucleus to use for cell count and load color images
          ## into the 4th dimention of a 4 dimention image array.  imagearray[pixel(x), pixel(y), color(RGB), image]
          ## This function is called by load images by default when count==TRUE
          
          DAPI <- grepl(pattern="dapi", imagedir, ignore.case=TRUE)
          DAPI <- which(DAPI)
          nuc <- grepl(pattern="nuc", imagedir, ignore.case=TRUE)
          nuc <- which(nuc)
          
          
          if(length(DAPI) == 0 & length(nuc) == 0){
                    stop("No folders exist containing DAPI or Nucleus images:  Make a folder named DAPI or Nucleus and add images for cell counting")
                    
          }
          
          else if(length(DAPI) >= 1 & length(nuc) >= 1){
                    stop("More than 1 folder exists for DAPI or Nucleus images:  Only 1 folder may exist labeled DAPI or Nucleus")
          }
          
          
          else if(length(DAPI) == 1){
                    DAPI <- grepl(pattern="dapi", imagefiles, ignore.case=TRUE)
                    location <- as.list(which(DAPI))
                    n <- length(location)
                    nucimages <<- combine(readImage(imagefiles[location[[1]]], colormode=Color))
                    for(i in 2:n){
                              nucimages <<- combine(nucimages, readImage(imagefiles[location[[i]]], colormode=Color))
                    }
                    
          }
          
          else if(length(nuc) == 1){
                    nuc <- grepl(pattern="nuc", imagefiles, ignore.case=TRUE)
                    location <- as.list(which(nuc))
                    n <- length(location)
                    nucimages <<- combine(readImage(imagefiles[location[[1]]], colormode=Color))
                    for(i in 2:n){
                              nucimages <<- combine(nucimages, readImage(imagefiles[location[[i]]], colormode=Color))
                    }
          }
          
          ## All images of DAPI/Nucleus are collapsed into 1 dimentional B/W images = the average of the intensities
          ## of each RGB pixel.
          
          nucimagesbw <<- (nucimages[,,1,1] + nucimages[,,2,1] + nucimages[,,3,1])/3
          for(i in 2:n){
                    nucimagesbw <<- combine(nucimagesbw, (nucimages[,,1,i] + nucimages[,,2,i] + nucimages[,,3,i])/3)
          }
          colorMode(nucimagesbw) <<- Grayscale
          
          
          ##Maximize contrast in all wholecellimagebw
          for(i in 1:n){
                    r <- range(nucimagesbw[,,,i])
                    nucimagesbw[,,,i] <<- nucimagesbw[,,,i] * (1/r[2])
                    display(nucimagesbw[,,,i])
          }
          
          
          ## Counting of nuclei using threshold only to estimate cell number **Generally a high estimate**
          
          nucimagesthreshold <<- nucimagesbw
          for(i in 1:n){
                    nucimagesthreshold[,,,i] <- nucimagesthreshold[,,,i] > threshold
          }
          
          
          ## Option for Hi pass filter prior to counts specified by hipass argument (default=FALSE)
          if(hipass == TRUE){
                    for(i in 1:n){
                              
                              fhi <- matrix(1, nc=3, nr=3)
                              fhi[2,2] <- -8
                              nucimagesthreshold[,,,i] <<- filter2(nucimagesthreshold[,,,1], fhi)
                    }
          }
          
          
          nuccounthi <<- nucimagesthreshold
          for(i in 1:n){
                    nuccounthi[,,,i] <<- bwlabel(nuccounthi[,,,i])
          }
          
          for(i in 1:n){
                    display(nuccounthi[,,,i]/max(nuccounthi))
                    if(writeimage == TRUE) writeImage(nuccounthi[,,,i]/max(nuccounthi), paste0("count_", basename(imagefiles[location[[i]]]),"_hi.tif"), "tiff")
          }
          
          
          ## Counting of nuclei using threshold/erode/dilate to estimate to include close cells or eliminate debris
          ## **Note** this may be a more acurate number to use, but images should first be examined as a guide.
          
          
          nuccountlow <<- nucimagesthreshold
          kern <- makeBrush(brush, shape="disc")
          
                    
          for(i in 1:n){
                    nuccountlow[,,,i] <<- dilate(erode(nuccountlow[,,,i], kern), kern)
          }
          
          for(i in 1:n){
                    nuccountlow[,,,i] <<- bwlabel(nuccountlow[,,,i])
          }
          
          for(i in 1:n){
                    display(nuccountlow[,,,i]/max(nuccountlow))
                    if(writeimage == TRUE) {writeImage(nuccountlow[,,,i]/max(nuccountlow), paste0("count_", basename(imagefiles[location[[i]]]), "_low.tif"), "tiff")}
          }
          
          
          ## Print the number of nuclei determined for each file.
          report <<- data.frame(file.name=0, low.count=0, high.count=0, low.image.output="<NA>", high.image.output="<NA>", stringsAsFactors = FALSE)
          for(i in 1:n){
                    report[i,1] <<- imagefiles[location[[i]]]
                    report[i,2] <<- max(nuccountlow[,,,i])
                    report[i,3] <<- max(nuccounthi[,,,i])
                    if(writeimage == TRUE) {report[i,4] <<- paste0("count_", basename(imagefiles[location[[i]]]),"_low.tif")}
                    if(writeimage == TRUE) {report[i,5] <<- paste0("count_", basename(imagefiles[location[[i]]]),"_hi.tif")}
                                                                  
                                                                  
          }
          if(writereport == TRUE) write.csv(report, file="report.csv")
          print(report)
}


outlinecells <- function(wholecell="STAT", wholethresh = .25, writeimage=FALSE, writereport=FALSE, ...){
          
          ## Find directory='wholecell' specified when calling the function that has images and load color images
          ## into the 4th dimention of a 4 dimention image array.  imagearray[pixel(x), pixel(y), color(RGB), image]
          
          #if(missing(wholecell)){
          #          stop(cat("A character string for 'wholecell=??' must be specified in function 'outlinecells()' cooresponding to the directory containing images of cell bodies"))
          #}
          
          whole <- grepl(pattern=wholecell, imagedir, ignore.case=TRUE)
          whole <- which(whole)
          
          if(length(whole) == 0){
                    stop(cat("No folders exist named", wholecell, ": Make a folder named", wholecell, "and add images coorespnding to those for nucleus counting"))
                    
          }
          
          else if(length(whole) > 1){
                    stop(cat("More than 1 folder exists named", wholecell, ": Only 1 folder may exist labeled", wholecell))
          }
          
          
          else if(length(whole) == 1){
                    whole <- grepl(pattern=wholecell, imagefiles, ignore.case=TRUE)
                    location <- as.list(which(whole))
                    n <- length(location)
                    wholecellimages <<- combine(readImage(imagefiles[location[[1]]], colormode=Color))
                    for(i in 2:n){
                              wholecellimages <<- combine(wholecellimages, readImage(imagefiles[location[[i]]], colormode=Color))
                    }
                    
          }
          
          ## All images of whole cell staining protein (STAT) are collapsed 
          ## into 1 dimentional B/W images = the average of the intensities of each RGB pixel.
          
          wholecellimagebw <<- (wholecellimages[,,1,1] + wholecellimages[,,2,1] + wholecellimages[,,3,1])/3
          for(i in 2:n){
                    wholecellimagebw <<- combine(wholecellimagebw, (wholecellimages[,,1,i] + wholecellimages[,,2,i] + wholecellimages[,,3,i])/3)
          }
          colorMode(wholecellimagebw) <<- Grayscale
          
          wholecellimagebwcontrast <<- wholecellimagebw
          
          ##Maximize contrast in all wholecellimagebw
          for(i in 1:n){
                    r <- range(wholecellimagebwcontrast[,,,i])
                    wholecellimagebwcontrast[,,,i] <<- wholecellimagebwcontrast[,,,i] * (1/r[2])
          }
          
          ## Setting the threshold that is used to outline cells currently set at > "wholethresh" quantile
          wholecellthresh <<- wholecellimagebwcontrast
          for(i in 1:n){
                    wholecellthresh[,,,i] <<- wholecellimagebwcontrast[,,,i] > quantile(wholecellimagebwcontrast[,,,i], wholethresh)
          }
          
          ctmask <<- wholecellthresh
          for(i in 1:n){
                    ctmask[,,,i] <<- opening(wholecellthresh[,,,i], makeBrush(5, shape='disc'))
          }
          
          if(writeimage == TRUE & writereport == TRUE){
                    report <- data.frame(ghost.cell.image="<NA>", outline.cell.nucleus.image="<NA>", stringsAsFactors = FALSE)
          }
          
          cmask <<- ctmask
          for(i in 1:n){
                    cmask[,,,i] <<- propagate(wholecellimagebwcontrast[,,,i], seeds=nuccountlow[,,,i], mask=ctmask[,,,i])
                    display(cmask[,,,i]/max(cmask[,,,i]))
                    if(writeimage == TRUE) {writeImage(cmask[,,,i]/max(cmask[,,,i]), paste0("ghost_", basename(imagefiles[location[[i]]]), "_low.tif"), "tiff")}
                    if(writereport == TRUE) {report[i,1] <- paste0("ghost_", basename(imagefiles[location[[i]]]), "_low.tif")}
          }
          
          
          ##Produce a composite of cellbody and nucleus
          
          composite <<- wholecellimages
          for(i in 1:n){
                    composite[,,3,i] <<- nucimagesbw[,,,i]              ##Makes Nucleus blue
                    composite[,,1,i] <<- wholecellimagebwcontrast[,,,i]         ##Makes wholecell red
                    composite[,,2,i] <<- 0                              ##Makes green channel 0
          }
          
          outlinenuccel <<- composite
          for(i in 1:n){
                    outlinenuccel[,,,i] <<- paintObjects(nuccountlow[,,,i], composite[,,,i], col='yellow')
                    outlinenuccel[,,,i] <<- paintObjects(cmask[,,,i], outlinenuccel[,,,i], col='purple')
                    display(outlinenuccel[,,,i])
                    if(writeimage == TRUE) {writeImage(outlinenuccel[,,,i], paste0("outline_nuc_", wholecell, "_", i, "_low.tif"), "tiff")}
                    if(writereport == TRUE) {report[i,2] <- paste0("outline_nuc_", wholecell, "_", i, "_low.tif")}
                    
          }
          
          if(writeimage == TRUE & writereport == TRUE){
                    counts <- read.csv(file="report.csv", row.names=1)
                    counts <- cbind(counts, report)
                    write.csv(counts, file="report.csv")
          } 
                  
}

colocalization <- function(selection="YFV", wholecell="STAT", colocBrush = 25, colocThresh = 0.75, writeimage=FALSE, writereport=FALSE, ...){
          
          # Find directory='selection' specified when calling the function that has images and load color images
          ## into the 4th dimention of a 4 dimention image array.  imagearray[pixel(x), pixel(y), color(RGB), image]
          
          #if(missing(selection)){
          #          stop(cat("A character string for 'wholecell=??' must be specified in function 'outlinecells()' cooresponding to the directory containing images of cell bodies"))
          #}
          
          select <- grepl(pattern=selection, imagedir, ignore.case=TRUE)
          select <- which(select)
          
          if(length(select) == 0){
                    stop(cat("No folders exist named", selection, ": Make a folder named", selection, "and add images coorespnding to those for nucleus counting"))
                    
          }
          
          else if(length(select) > 1){
                    stop(cat("More than 1 folder exists named", selection, ": Only 1 folder may exist labeled", selection))
          }
          
          
          else if(length(select) == 1){
                    select <- grepl(pattern=selection, imagefiles, ignore.case=TRUE)
                    location <- as.list(which(select))
                    n <- length(location)
                    selectionimages <<- combine(readImage(imagefiles[location[[1]]], colormode=Color))
                    for(i in 2:n){
                              selectionimages <<- combine(selectionimages, readImage(imagefiles[location[[i]]], colormode=Color))
                    }
                    
          }
          
          ## All images of "selection" are collapsed into 1 dimentional B/W images = the average of the intensities
          ## of each RGB pixel.
          
          selectionimagesbw <<- (selectionimages[,,1,1] + selectionimages[,,2,1] + selectionimages[,,3,1])/3
          for(i in 2:n){
                    selectionimagesbw <<- combine(selectionimagesbw, (selectionimages[,,1,i] + selectionimages[,,2,i] + selectionimages[,,3,i])/3)
          }
          colorMode(selectionimagesbw) <<- Grayscale
          
          ##Maximize contrast in all wholecellimagebw
          
          selectionimagesbwcontrast <<- selectionimagesbw
          for(i in 1:n){
                    r <- range(selectionimagesbwcontrast[,,,i])
                    selectionimagesbwcontrast[,,,i] <<- selectionimagesbwcontrast[,,,i] * (1/r[2])
                    display(selectionimagesbwcontrast[,,,i])
          }
          
          ## Threshold then Erode away "selection" background and leave positive cells.
          selectionimagespos <<- selectionimagesbwcontrast
          kern <- makeBrush(colocBrush, shape="disc")
          for(i in 1:n){
                    selectionimagespos[,,,i] <<- erode(selectionimagespos[,,,i] > quantile(selectionimagespos[,,,i], colocThresh), kern)
                    display(selectionimagespos[,,,i])
          }
          
          ## Ask whether positive areas overlap with identified whole cells and save the identifying cell number in a vector "cell"
          ## col#='image', row#='cell identification', 0='no overlap', 1='overlap'.  The sum of a column is the number of positive cells
          cell <<- matrix(data=NA, nrow=max(cmask), ncol=n)
          sums <<- data.frame(selection="<NA>", stringsAsFactors=FALSE)
          colnames(sums) <<- paste0("number.of.", selection, " pos.cells")
          for(i in 1:n){
                    m <- max(cmask[,,,i])
                    for(c in 1:m){
                              cmasktmp <- cmask[,,,i] == c
                              cmasktmp <- selectionimagespos[,,,i][cmasktmp]
                              if(sum(cmasktmp >= 1)){
                                        cell[c,i] <<- 1
                              }
                              else{
                                        cell[c,i] <<- 0
                              }
                    }
                    sums[i,1] <<- sum(cell[,i], na.rm=TRUE)
          }
          
          print(sums)
          
          if(writereport == TRUE){
                    report <- read.csv("report.csv", row.names=1)
                    report <- cbind(report, sums)
                    write.csv(report, file="report.csv")
          }
          
          write.csv(cell, file="cell.csv")

          
}


overlay <- function(writeimage=FALSE, writereport=FALSE, ...){
          ##Produce a composite of cellbody and nucleus
          
          n <- dim(nucimagesbw)[4]
          composite <<- wholecellimages
          if(writereport == TRUE){names <- data.frame(overlay.image.name="<NA>", stringsAsFactors=FALSE)}
          
          for(i in 1:n){
                    composite[,,3,i] <<- nucimagesbw[,,,i]                      ##Makes Nucleus blue
                    composite[,,1,i] <<- wholecellimagebwcontrast[,,,i]         ##Makes wholecell red
                    composite[,,2,i] <<- selectionimagesbwcontrast[,,,i]        ##Makes selection green
                    display(composite[,,,i])
                    if(writeimage == TRUE) {writeImage(composite[,,,i], paste0("overlay_", i, ".tif"), "tiff")}
                    if(writereport == TRUE){names[i,] <- paste0("overlay_", i, ".tif")}
          }
          
          if(writereport == TRUE){
                    report <- read.csv("report.csv", row.names=1)
                    report <- cbind(report, names)
                    write.csv(report, file="report.csv")
          }
          
}


statQuant <- function(hipassQuant = FALSE, ...){
          
          ## Take as average of whole cell marker intensity infected and non-infected cells and report those values.
          
          imageNum <- ncol(cell)
          imageCol <- imageNum*4
          imageCol2 <- imageNum*2
          imageCol3 <- imageNum*3
          rows <- nrow(cell)
          expressionInfected <<- matrix(ncol=imageCol, nrow=rows)
          expressionMeans <<- matrix(ncol=imageCol2, nrow=rows)
          expressionSums <<- matrix(ncol=imageCol2, nrow=rows)
          expressionArea <<- matrix(ncol=imageCol3, nrow=rows)
          colNames <<- character()
          colMean <<- character()
          colSum <<- character()
          colArea <<- character()
          for(a in 1:imageNum){
                    
                    l <- length(colNames)
                    m <- length(colMean)
                    s <- length(colSum)
                    z <- length(colArea)
                    
                    colNames[l+1] <<- paste0("image.", a, ".uninfect.stat.sum")
                    colNames[l+2] <<- paste0("image.", a, ".uninfect.stat.mean")
                    colNames[l+3] <<- paste0("image.", a, ".infect.stat.sum")
                    colNames[l+4] <<- paste0("image.", a, ".infect.stat.mean")
                    
                    colMean[m+1] <<- paste0("image.", a, ".stat.mean")
                    colMean[m+2] <<- paste0("image.", a, ".YFV.mean")
                    
                    colSum[s+1] <<- paste0("image.", a, ".stat.sum")                    
                    colSum[s+2] <<- paste0("image.", a, ".YFV.sum")
                    
                    colArea[z+1] <<- paste0("image.", a, ".cell.pixel.area")
                    colArea[z+2] <<- paste0("image.", a, ".stat.per.area")
                    colArea[z+3] <<- paste0("image.", a, ".yfv.per.area")
                    
                    
          }
          
          colnames(expressionInfected) <<- colNames
          colnames(expressionMeans) <<- colMean
          colnames(expressionSums) <<- colSum
          colnames(expressionArea) <<- colArea
          
          expressionInfected <<- as.data.frame(expressionInfected)
          expressionMeans <<- as.data.frame(expressionMeans)
          expressionSums <<- as.data.frame(expressionSums)
          expressionArea <<- as.data.frame(expressionArea)
          
          ## Option for Hi pass filter prior to intensity acquisition from selection (default=FALSE)
          if(hipassQuant == TRUE){
                    n <- dim(selectionimagesbw)[4]
                    for(i in 1:n){
                              
                              fhi <- matrix(1, nc=3, nr=3)
                              fhi[2,2] <- -7
                              selectionimagesbw[,,,i] <<- filter2(selectionimagesbw[,,,i], fhi)
                    }
          }
                    
          for(u in 1:imageNum){
                    
                    cellNum <- sum(!is.na(cell[,u]))
                    
                    for(i in 1:cellNum){
                                        cellLocation <- cmask[,,,u] == i
                                        
                                        cellArea <- sum(cellLocation)
                                        statSum <- sum(wholecellimagebw[,,,u][cellLocation])
                                        statMean <- mean(wholecellimagebw[,,,u][cellLocation])
                                        YFVSum <- sum(selectionimagesbw[,,,u][cellLocation])
                                        YFVMean <- mean(selectionimagesbw[,,,u][cellLocation])
                                        
                                        if(cell[i,u] == 0){
                                                  expressionInfected[i, (u*4)-3] <<- statSum
                                                  expressionInfected[i, (u*4)-2] <<- statMean
                                        }
                                        else if(cell[i,u] == 1){
                                                  
                                                  expressionInfected[i, (u*4)-1] <<- statSum
                                                  expressionInfected[i, (u*4)] <<- statMean
                                                  
                                        }
                                        
                                        expressionMeans[i, (u*2)-1] <<- statMean
                                        expressionMeans[i, (u*2)] <<- YFVMean
                                        
                                        expressionSums[i, (u*2)-1] <<- statSum
                                        expressionSums[i, (u*2)] <<- YFVSum
                                        
                                        expressionArea[i, (u*3)-2] <<- cellArea
                                        expressionArea[i, (u*3)-1] <<- statSum / cellArea
                                        expressionArea[i, (u*3)] <<- YFVSum / cellArea

                    }
                    
                              
          }
          
          write.csv(expressionInfected, "expressionInfected.csv")
          write.csv(expressionMeans, "expressionMeans.csv")
          write.csv(expressionSums, "expressionSums.csv")
          write.csv(expressionArea, "expressionArea.csv")
          
          
          expressionInfectedSummary <- matrix(ncol=imageCol, nrow=3)
          row.names(expressionInfectedSummary) <- c("number.of.cells", "average", "standard.deviation")
          colnames(expressionInfectedSummary) <- colNames
          
          for(i in 1:imageCol){
                    
                    number <- sum(!is.na(expressionInfected[,i]))
                    average <- mean(expressionInfected[,i], na.rm=TRUE)
                    stddev <- sd(expressionInfected[,i], na.rm=TRUE)
                    
                    expressionInfectedSummary[1,i] <- number
                    expressionInfectedSummary[2,i] <- average
                    expressionInfectedSummary[3,i] <- stddev
                    
          }
          
          write.csv(expressionInfectedSummary, "expressionInfectedSummary.csv")
          
}

exportIndCells <- function(...){
          
          colNam <<- colnames(expressionInfected)
          meansValues <<- colNam[grepl(pattern="mean", colNam, ignore.case=TRUE)]
          l <<- length(meansValues)
          r <<- nrow(expressionInfected)
          organizedInfected <<- matrix(ncol=l, nrow=r)
          colnames(organizedInfected) <<- colnames(expressionInfected[meansValues])
          for(i in 1:l){
                    
                    wanted <<- !is.na(expressionInfected[,meansValues[i]])
                    wanted <<- expressionInfected[wanted, meansValues[i]]
                    wantedLength <- length(wanted)
                    for(a in 1:wantedLength){
                              
                              organizedInfected[a,i] <<- wanted[a]
                    }
                    write.csv(organizedInfected, "organizedInfected.csv")
          }
}


######################################################################################################
## Script to load data files produced from running previous functions                               ##
## This can be used after running an analysis and wishing to revist plotting of acquired data       ##
## or wishing to finish an incomplete analysis that requires only the data dumped to disk           ##
######################################################################################################

loadDataFiles <- function(...){
          
          if(file.exists("cell.csv") == TRUE){
                    cell <<- read.csv("cell.csv", row.names=1)
          }
          
          if(file.exists("expressionInfected.csv") == TRUE){
                    expressionInfected <<- read.csv("expressionInfected.csv", row.names=1)
          }
          if(file.exists("expressionInfectedSummary.csv") == TRUE){
                    expressionInfectedSummary <<- read.csv("expressionInfectedSummary.csv", row.names=1)
          }
          if(file.exists("expressionMeans.csv") == TRUE){
                    expressionMeans <<- read.csv("expressionMeans.csv", row.names=1)
          }
          if(file.exists("expressionSums.csv") == TRUE){
                    expressionSums <<- read.csv("expressionSums.csv", row.names=1)
          }
          if(file.exists("forPlot.csv") == TRUE){
                    forPlot <<- read.csv("forPlot.csv", row.names=1)
          }
          if(file.exists("forPlot2.csv") == TRUE){
                    forPlot2 <<- read.csv("forPlot2.csv", row.names=1)
          }
          if(file.exists("organizedInfected.csv") == TRUE){
                    organizedInfected <<- read.csv("organizedInfected.csv", row.names=1)
          }
          if(file.exists("report.csv") == TRUE){
                    report <<- read.csv("report.csv", row.names=1)
          }
}


#######################################################################################################################
## Places all values for expressionMeans into single columns for easy plotting.  Each cell is identified by the      ##
## image that it was obtained from, whether it was consider to be infected and the pixel intensities of stat and YFV ##
#######################################################################################################################

organizeForPlot <- function(...){
          
          forPlot <<- data.frame(imagenum=0, infected=0, stat=0, yfv=0, stringsAsFactors=FALSE)
          
          imageNum <- ncol(cell)
          
          n <- 1
          for(i in 1:imageNum){
                    
                    cellNum <- sum(!is.na(cell[,i]))
                    
                    for(p in 1:cellNum){
                              forPlot[n, 1] <<- i
                              forPlot[n, 2] <<- cell[p, i]
                              forPlot[n, 3] <<- expressionMeans[p,(i*2)-1]
                              forPlot[n, 4] <<- expressionMeans[p, (i*2)]
                              n <- n+1
                    }
          }
          
          write.csv(forPlot, "forPlot.csv")
}