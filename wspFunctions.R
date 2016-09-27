require(XML)

########################
#### findGateInWspGates(gateList, wspGates)
########################
findGateInWspGates = function (gateList, wspGates){
  wspCol = 2
  for(wspCol in c(1:ncol(wspGates))){
    if(compareGates(gate1=gateList, gate2=wspGates[,wspCol]))
      print(wspG)
  }
}

########################
######### compareGates(gate1, gate2)
########################

compareGates = function(gate1, gate2){
  if(gate1$gateType != gate2$gateType) {return (FALSE)}
  if(length(gate1$dimNames) != length(gate2$gateValues)) {return(FALSE)}
  idx=1
  for (idx in c(1:length(gate1$dimNames))){
    for (valueName in names(gate2$gateValues[gate1$dimNames[idx]][[1]]) ){
      class(wspGates[,2])
      gate1$gateDimensions[[idx]]
      
    }
  }
}

########################
#### extractGates(xml)
########################

extractGates = function(x) {
  populationName = xmlAttrs(x)[["name"]]
  gateType=names(x[['Gate']])
  switch(gateType ,
         EllipsoidGate={
           print(paste("found EllipsoidGate"))
           list(populationName=populationName, gateType=gateType,gateValues=NA)
         },
         RectangleGate={
           dimIdxs=grep('dimension',names(x[['Gate']][[1]]))
           gateNames = rep(NA,length(dimIdxs))
           gateValues = rep(list(),length(dimIdxs))
           gnIdx=1
           dimIdx=1
           for(dimIdx in dimIdxs){
             gateNames[gnIdx] = as.character(
               xmlAttrs(
                 x[['Gate']][[gateType]][[dimIdx]][['fcs-dimension']]
               )
             )
             gateAttrsNames = names(xmlAttrs(x[['Gate']][[gateType]][[dimIdx]]))
             tmpValues = rep(NA,length(gateAttrsNames))
             atIdx=1
             for(attr in gateAttrsNames){
               tmpValues[atIdx]=as.numeric(as.character(xmlAttrs(x[['Gate']][[gateType]][[dimIdx]])[[attr]]))
               atIdx=atIdx+1
             }
             names(tmpValues)=gateAttrsNames
             gateValues[gnIdx] = list(tmpValues)
             gnIdx=gnIdx+1
           }
           names(gateValues)=gateNames
           
           return(list(populationName=populationName, gateType=gateType,gateValues=gateValues))
           
         },
         PolygonGate={
           # get the positions of the dimension names from the first gate object. We assume that there can only be one.
           dimIdxs=grep('dimension',names(x[['Gate']][[1]]))
           gateNames = rep(NA,length(dimIdxs))
           gateValues = rep(list(),length(dimIdxs))
           gnIdx=1
           dimIdx=1
           # get the fcs dimension names
           for(dimIdx in dimIdxs){
             gateNames[gnIdx] = as.character(
               xmlAttrs(
                 x[['Gate']][[gateType]][[dimIdx]][['fcs-dimension']]
               )
             )
             gnIdx=gnIdx+1
           }
           # A Polygon gate has many gating:vertices entries with two values for each of the dimensions previously defined.
           # get the indices in the x object (again we assume that there is only one gate)
           vertices=grep('vertex',names(x[['Gate']][[1]]))
           vIdx=3
           gateVertices = matrix(ncol=2,nrow=length(vertices))
           count=1
           for(vIdx in vertices){
             gateVertices[count,1]=as.numeric(as.character(xmlAttrs(x[['Gate']][[gateType]][[vIdx]][[1]])[[1]]))
             gateVertices[count,2]=as.numeric(as.character(xmlAttrs(x[['Gate']][[gateType]][[vIdx]][[2]])[[1]]))
             count=count+1
           }
           colnames(gateVertices)=gateNames
           return(list(populationName=populationName, gateType=gateType,gateValues=gateVertices))
           },
         {
           print(paste("not found gate type: ",gateType))
           return(NA)
         }
  )
  #print(gateType)
  #assign("globalIdx", globalIdx+1, envir = .GlobalEnv)
  return(NA)
  }


########################
########## readWsp(fileName)
########################
readWsp = function(fileName, sampleName="001 160518G29a1 SCM WLSM.fcs"){
  data=xmlTreeParse(fileName,useInternalNodes = T)
  populations =getNodeSet(data,paste0("//SampleNode[@name='",sampleName,"']//Population"))
  #populations = getNodeSet(data,"//Population")
  #extractGates(x=populations[[1]])
  wspGates = xmlSApply(populations,FUN =extractGates,simplify=T )
  colnames(wspGates)=unlist(wspGates['populationName',])
  return(wspGates)
}


########################
########## xmlEntry2internalGatelist(gating_xml)
########################

xmlEntry2internalGatelist = function(gating_xml){
  xmlGateType=names(gating_xml)
  
  switch(xmlGateType ,
         RectangleGate={
           dimIdxs=grep('dimension',names(gating_xml[[xmlGateType]]))
           gateNames = rep(NA,length(dimIdxs))
           gateValues = rep(list(),length(dimIdxs))
           gnIdx=1
           dimIdx=1
           for(dimIdx in dimIdxs){
             gating_xml[[xmlGateType]][[dimIdx]][['fcs-dimension']][['name']]
             gateNames[gnIdx] = gating_xml[[xmlGateType]][[dimIdx]][['fcs-dimension']][['name']]
             gateAttrsNames = names(gating_xml[[xmlGateType]][[dimIdx]]$.attrs)
             tmpValues = rep(NA,length(gateAttrsNames))
             atIdx=1
             for(attr in gateAttrsNames){
               tmpValues[atIdx]=as.numeric(gating_xml[[xmlGateType]][[dimIdx]]$.attrs[[attr]])
               atIdx=atIdx+1
             }
             names(tmpValues)=gateAttrsNames
             gateValues[gnIdx] = list(tmpValues)
             gnIdx=gnIdx+1
           }
           names(gateValues)=gateNames
           xmlGate= list(gateType=xmlGateType, 
                         gateValues=gateValues)
           return(xmlGate)
         },
         PolygonGate={
           dimIdxs=grep('dimension',names(gating_xml[[xmlGateType]]))
           gateNames = rep(NA,length(dimIdxs))
           gateValues = rep(list(),length(dimIdxs))
           gnIdx=1
           dimIdx=1
           for(dimIdx in dimIdxs){
             gating_xml[[xmlGateType]][[dimIdx]][['fcs-dimension']][['name']]
             gateNames[gnIdx] = gating_xml[[xmlGateType]][[dimIdx]][['fcs-dimension']][['name']]
             gnIdx=gnIdx+1
           }
           
           vertices=grep('vertex',names(gating_xml[[xmlGateType]]))
           vIdx=3
           gateVertices = matrix(ncol=2,nrow=length(vertices))
           count=1
           for(vIdx in vertices){
             gateVertices[count,1]=as.numeric(as.character(xmlAttrs(x[['Gate']][[gateType]][[vIdx]][[1]])[[1]]))
             gateVertices[count,2]=as.numeric(as.character(xmlAttrs(x[['Gate']][[gateType]][[vIdx]][[2]])[[1]]))
             count=count+1
           }
           
           colnames(gateVertices)=gateNames
           xmlGate = list(gateType=gateType,gateValues=gateVertices)
           return(xmlGate)
         },
         {
           print(paste("not found gate type: ",gateType))
           
         })
  return(NA)
}

########################
########## compare2Gates(xmlGate, wspGate)
########################

compare2Gates = function(xmlGate, wspGate){
  #print("XmlGate")
  #print(xmlGate)
  #print("WspGate")
  #print(wspGate)
  if(xmlGate$gateType != wspGate$gateType){return (FALSE)}
  #if(xmlGate$gateType != "RectangleGate"){return (FALSE)}
  
  #xmlGate$gateValues[1]
  if(xmlGate$gateType =="RectangleGate"){
    gateVIdx=1
    for(gateVIdx in c(1:length(xmlGate$gateValues))){
      wspValues=wspGate$gateValues[ names(xmlGate$gateValues[gateVIdx])][[1]]
      if(is.null(wspValues)) {return(FALSE)}
      xmlValues=xmlGate$gateValues[gateVIdx][[1]]
      #print("wspValue")
      #print(wspValues)
      #print("xmlValues")
      #print(xmlValues)
      #print(names(xmlValues[1]))
      for(thIdx in c(1:length(xmlValues))){
        if(is.na(wspValues[names(xmlValues[thIdx])])){return(FALSE)}
        if(wspValues[names(xmlValues[thIdx])] != xmlValues[thIdx]){
          return(FALSE)
        }
      }
    }
    return(TRUE)
  }
  
  if(xmlGate$gateType == "PolygonGate"){
    if(!all(colnames(wspGate$gateValues) == colnames(xmlGate$gateValues))){return(FALSE)}
    if (all(wspGate$gateValues[,1] == xmlGate$gateValues[,1]) & all(wspGate$gateValues[,2] == xmlGate$gateValues[,2])){
      return(TRUE)
    }
    return(FALSE)
  }
  return(FALSE)
}

########################
########## shiftFunction(x,n)
########################
shiftFunction <- function(x,n){
  c(x[(n+1):length(x)],x[1:n])
}

########################
########## findName(xmlGate, wspData)
########################

findName = function (xmlGate, wspData){
  wspIdx=1
  for (wspIdx in c(1:ncol(wspData))){
    wspGate = wspData[,wspIdx]
    if(compare2Gates(xmlGate=xmlGate, wspGate=wspGate)) {
      #print("Found")
      #print(wspGate$populationName)
      return(wspGate$populationName)
    }
  }
  return(NA)
}

PlotNode = function (fsom, 
                     id, 
                     markers = fsom$map$colsUsed, 
                     colorPalette = grDevices::colorRampPalette(c("#00007F", "blue", 
                                                                  "#007FFF", "cyan", 
                                                                  "#7FFF7F", "yellow", 
                                                                  "#FF7F00", "red", "#7F0000")), 
                     main = paste0("Cluster ", id),
                     plotLables=TRUE,
                     magnifier=1, ...) {
  coords=matrix(c(100, 0, 100, 100, 0, 100), nrow = 3)
  graphics::plot(1, type = "n", xlab = "", ylab = "", 
                 xlim = c(-10, 10), ylim = c(-3, 3), 
                 asp = 1, bty = "n", xaxt = "n", 
                 yaxt = "n", main = main)
  labels <- fsom$prettyColnames[markers]
  n <- length(markers)
  if (is.function(colorPalette)) {
    colorPalette <- colorPalette(n)
  }
  data <- rbind(apply(fsom$map$medianValues[, markers, drop = FALSE], 
                      2, min), fsom$map$medianValues[id, markers, drop = FALSE], 
                apply(fsom$map$medianValues[, markers, drop = FALSE], 
                      2, max))
  print(data)
  graphics::stars(data*magnifier, col.segments = colorPalette, locations = coords, 
                  draw.segments = TRUE, add = TRUE, inches = FALSE, len=magnifier, ...)
  graphics::symbols(coords[, 1], coords[, 2], circles = c(magnifier, magnifier, magnifier), 
                    inches = FALSE, bg = "transparent", bty = "n", 
                    add = TRUE)
  if(plotLables){
    angle <- 2 * pi/n
    angles <- seq(angle/2, 2 * pi, by = angle)
    left <- (angles > (pi/2) & angles < (3 * pi/2))
    x <- c(2, -2)[left + 1]
    y_tmp <- c(seq(-2, 2, by = 4/(sum(!left) + 1))[-c(1, sum(!left) + 
                                                        2)], seq(2, -2, by = -4/(sum(left) + 1))[-c(1, sum(left) + 
                                                                                                      2)])
    y <- shiftFunction(y_tmp, max((cummax(y_tmp) < 0) * seq_along(y_tmp)))
    for (i in seq_along(labels)) {
      graphics::text(x = x[i]*magnifier, y = y[i]*magnifier, labels = labels[i], 
                     adj = c(as.numeric(left)[i], 0.5), cex = 0.5)
      graphics::lines(x = magnifier * c(x[i] + c(-0.2, 0.2)[left[i] + 1],c(1.5, -1.5)[left[i] + 1], cos(angles[i])), 
                      y = magnifier * c(y[i], y[i], sin(angles[i])), 
                      col = colorPalette[i], lwd = 2)
    }
  }
}
