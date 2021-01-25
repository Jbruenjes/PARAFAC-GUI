#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(plotly)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    
    #increase upload size limit
    options(shiny.maxRequestSize=30*1024^2)
    
    #load in packages
    library(plotly)
    library(tidyverse)
    library(staRdom)
    library(DT)

    #define functions
    #define scatter removal function
    remove_scatter<-function(eem_list)
    {
        
        ex<-eem_list$ex
        em<-eem_list$em
        x<-eem_list$x
        #taken from EEMR package and adjusted for input below AND above scatter region with two independent values
        
        ray1 <- eem_list[["scatter"]][["ray1"]]
        ray2 <- eem_list[["scatter"]][["ray2"]]
        raman1 <- eem_list[["scatter"]][["raman1"]]
        raman2 <- eem_list[["scatter"]][["raman2"]]
        
        ind1 <- mapply(function(x) em <= x, 1 * ex + ray1[1])
        ind2 <- mapply(function(x) em <= x, 1 * ex + ray1[2])
        ind3 <- ifelse(ind1 + ind2 == 1, NA, 1)
        x <- x * ind3
        
        ind1 <- mapply(function(x) em <= x, 2 * ex + ray2[1])
        ind2 <- mapply(function(x) em <= x, 2 * ex + ray2[2])
        ind3 <- ifelse(ind1 + ind2 == 1, NA, 1)
        x <- x * ind3
        #raman
        ex_raman <- -(ex / (0.00036 * ex - 1))
        
        ind1 <- mapply(function(x) em <= x, 1 * ex_raman + raman1[1])
        ind2 <- mapply(function(x) em <= x, 1 * ex_raman + raman1[2])
        ind3 <- ifelse(ind1 + ind2 == 1, NA, 1)
        x <- x * ind3
        
        ind1 <- mapply(function(x) em <= x, 2 * ex_raman + raman2[1])
        ind2 <- mapply(function(x) em <= x, 2 * ex_raman + raman2[2])
        ind3 <- ifelse(ind1 + ind2 == 1, NA, 1)
        x <- x * ind3
        
        
        eem_list$x <-x
        
        attr<-attributes(eem_list)
        attr$is_scatter_corrected <- TRUE
        attributes(eem_list)<-attr
        
        eem_list
        
    }
    
    #define eem cutting function
    cuteem<-function(x, range.em.ind.min, range.em.ind.max,range.ex.ind.min,range.ex.ind.max){
        
        x$x<-x$x[c(range.em.ind.min:range.em.ind.max),
                 c(range.ex.ind.min:range.ex.ind.max)]
        x$ex<-x$ex[range.ex.ind.min:range.ex.ind.max]
        x$em<-x$em[range.em.ind.min:range.em.ind.max]
        
        x
    }
    #show sse of individual samples
    
    residuals_sse<-function (pfmodel, eem_list, select = NULL, cores = parallel::detectCores(logical = FALSE)/2) 
    {
        #originallly from staRdom package.
        #removed code to bind values for compounds and samples to reduce calculation time. original function from staRdom
        pfmodel <- norm2A(pfmodel)
        if (!is.null(select)) {
            eem_list <- eem_extract(eem_list, sample = select, keep = TRUE, 
                                    verbose = FALSE)
        }
        if (!all(eem_names(eem_list) %in% rownames(pfmodel$A)) | 
            length(eem_list) == 0) {
            pfmodel <- A_missing(eem_list, pfmodel, cores = cores)
        }
        what <- which(rownames(pfmodel$A) %in% (eem_list %>% eem_names()))
        pfmodel$A <- as.data.frame(pfmodel$A)[what, ]
        res_data <- lapply(pfmodel$A %>% rownames(), function(sample) {
            comps <- lapply(pfmodel$A %>% colnames(), function(component) {
                pfmodel$B[, component] %*% t(pfmodel$C[, component]) * 
                    pfmodel$A[sample, component]
            })
            names(comps) <- pfmodel$C %>% colnames()
            fit <- comps %>% Reduce("+", .)
            eem <- eem_list[[which(eem_list %>% eem_names == sample)]]
            samp <- eem$x[eem$em %in% rownames(pfmodel$B), eem$ex %in% 
                              rownames(pfmodel$C)]
            res <- samp - fit
            
            colnames(samp) <- rownames(pfmodel$C)
            
            rownames(res) <- rownames(pfmodel$B)
            colnames(res) <- rownames(pfmodel$C)
            res <- res %>% data.frame() %>% mutate(type = "residual", 
                                                   em = rownames(pfmodel$B)) %>% gather(ex, value, -em, 
                                                                                        -type) %>% mutate(ex = substr(ex, 2, 4))
            bind_rows(list(res)) %>% mutate(Sample = sample)
        }) %>% bind_rows()
    }
    
    ##########
    spectralvarianceplot<- function(eem_list){
        range.ex<-eem_list[[1]]$ex
        range.em<-eem_list[[1]]$em
        
        
        data<-lapply(eem_list,"[[",3) %>% lapply(as.vector) %>% as.data.frame()
        
        csum<-colSums(data, na.rm=T)
        
        mat_norm<-vector("list", length(eem_list))
        
        for (i in 1:length(mat_norm)){
            mat_norm[[i]]<- data[,i]/ csum[i]
        }
        
        mat_norm<- as.data.frame(mat_norm) 
        samplenames<- sapply(eem_list, function(x){x$sample})
        
        
        sample_ind_vec<-apply(mat_norm, FUN=which.max, MARGIN =1) %>% as.integer()
        sample_name_vec<-sapply(sample_ind_vec, function(x){ x<-samplenames[x]})
        
        
        spectralvariance_matrix<-matrix(data=apply(mat_norm,MARGIN = 1, sd, na.rm = T),nrow=length(range.em),
                                        ncol=length(range.ex))
        
        
        spectralvariance_names<-matrix(data= paste0("Sample ",sample_ind_vec, " ", sample_name_vec),
                                       nrow=length(range.em),
                                       ncol=length(range.ex))
        
        
        
        #contour plot the stdev of all eems
        plot_ly(
            x =  range.ex,
            y =  range.em,
            z = spectralvariance_matrix,
            text=spectralvariance_names,
            type = "contour",contours=list(coloring ="heatmap"),ncontours=50,
            colorscale = "Viridis",line=list(color="black", width=0.4)) %>% 
            layout ( title = "Spectral Variance of all EEMs",
                     xaxis = list ( title = "Ex [nm]"), 
                     yaxis = list( title = "Em [nm]"))
        
    }
    
    ####
    read_absorbance_shiny <- function(abs_data, order = TRUE, recursive = TRUE, dec = NULL, 
                                      sep = NULL, verbose = FALSE, cores = parallel::detectCores(logical = FALSE)){
        #From staRdom, adjusted to take a list of uploaded files as input argument.
        
        
        cl <- makeCluster(min(cores, length(abs_data)), type = "PSOCK")
        clusterExport(cl, c("dec", "sep", "verbose"), 
                      envir = environment())
        clusterEvalQ(cl, require(dplyr))
        clusterEvalQ(cl, require(stringr))
        abs_data <- parLapply(cl, abs_data, function(tab) {
            tryCatch({
                rawdata <- readLines(tab)
                data <- rawdata %>% sapply(str_remove, pattern = "([^0-9]*$)")
                first_number <- min(which((substr(data, 1, 1) %>% 
                                               grepl("[0-9]", .))))
                last_number <- max(which((substr(data, 1, 1) %>% 
                                              grepl("[0-9]", .))))
                if (is.null(sep) | is.null(dec)) {
                    nsepdec <- data[first_number] %>% str_extract_all("[^-0-9eE]") %>% 
                        unlist()
                    example_number <- data[first_number] %>% str_extract("([-]?[0-9]+[.,]?[0-9]+[eE]?[-0-9]+)$")
                    if (is.null(dec) & length(nsepdec) > 1) 
                        dec <- example_number %>% str_replace("([-0-9eE]+)([.,]?)([-0-9eE]*)", 
                                                              "\\2")
                    if (is.null(sep)) 
                        sep <- gsub(pattern = dec, replacement = "", 
                                    x = data[first_number], fixed = TRUE) %>% 
                        str_extract(paste0("[^-0-9eE", dec, 
                                           "]"))
                    if (verbose) 
                        warning("processing", tab, ": using", 
                                sep, "as field separator and", dec, 
                                "as decimal separator.", fill = TRUE)
                }
                data <- str_split(data, sep)
                table <- data[(first_number):last_number] %>% unlist() %>% 
                    matrix(ncol = length(data[[first_number]]), byrow = TRUE) %>% 
                    data.frame(stringsAsFactors = FALSE) %>% mutate_all(gsub, 
                                                                        pattern = ifelse(dec != "", dec, "."), 
                                                                        replacement = ".", fixed = TRUE)
                table <- table %>% mutate_all(as.numeric)
                attr(table, "location") <- rep(tab, ncol(table) - 
                                                   1)
                if (ncol(table) == 2) {
                    samples <- tab %>% basename() %>% str_replace_all(regex(".txt$|.csv$", 
                                                                            ignore_case = TRUE), "")
                }
                else {
                    samples <- rawdata[[1]] %>% str_split(sep) %>% 
                        unlist() %>% matrix(ncol = length(.), byrow = TRUE) %>% 
                        data.frame(stringsAsFactors = FALSE) %>% .[-1]
                }
                table <- table %>% setNames(c("wavelength", 
                                              samples))
            }, error = function(err) {
                stop("Error while reading ", tab, ": ", 
                     err)
            })
        })
        stopCluster(cl)
        locations <- lapply(abs_data, function(tab) {
            attr(tab, "location")
        }) %>% unlist()
        if (length(abs_data) == 1) 
            abs_data <- abs_data[[1]] %>% as.data.frame()
        else abs_data <- abs_data %>% list_join(by = "wavelength")
        if (order) 
            abs_data <- abs_data %>% arrange(wavelength)
        attr(abs_data, "location") <- locations
        abs_data
    }
    ###
    ## from staRDOM package, adjusted to hide print
    plot_sh<-function (fits) 
    {
        sel <- 0
        table <- lapply(fits, function(fit) {
            sel <<- sel + 1
            c <- fit %>% lapply(eempf_comp_mat)
            tab <- lapply(c, function(c1) {
                nc1 <- length(c1)
                nc2 <- 0
                lapply(c1, function(c2) {
                    nc2 <<- nc2 + 1
                    c2 <- c2 %>% mutate(comps = nc1, comp = paste0("Comp.", 
                                                                   nc2))
                }) %>% bind_rows()
            }) %>% bind_rows() %>% mutate(selection = sel) %>% group_by(comps, 
                                                                        comp) %>% mutate(max_pos = which.max(value), max_em = em[max_pos], 
                                                                                         max_ex = ex[max_pos]) %>% mutate(exn = ifelse(em == 
                                                                                                                                           max_em, ex, NA), emn = ifelse(ex == max_ex, em, NA)) %>% 
                filter(!is.na(emn) | !is.na(exn)) %>% mutate(ex = exn, 
                                                             em = emn) %>% select(-exn, -emn, -max_pos, -max_em, 
                                                                                  -max_ex) %>% ungroup() %>% mutate_at(vars(em, ex, 
                                                                                                                            value), as.numeric)
        }) %>% bind_rows()
        pl1 <- table %>% mutate(selection = factor(selection, ordered = FALSE)) %>% 
            ggplot() + geom_line(data = . %>% filter(!is.na(ex)), 
                                 aes(x = ex, y = value, colour = selection, group = selection), 
                                 linetype = 2) + geom_line(data = . %>% filter(!is.na(em)), 
                                                           aes(x = em, y = value, colour = selection, group = selection), 
                                                           linetype = 1) + labs(x = "Wavelength (nm)", y = "Loading") + 
            theme(legend.position = "none", axis.text.x = element_text(angle = 90, 
                                                                       hjust = 1)) + facet_grid(. ~ comp)
        pl1
    }
    ##
    eem_plot_contour_lines <- function (eem_list, s_ind, event.data, includelines, source){
        # isolate sample ex/em and value
        ex<-eem_list[[s_ind]]$ex
        em<-eem_list[[s_ind]]$em
        x<-eem_list[[s_ind]]$x
        
        if (includelines){
            
            #Distinguish between eem + click on eem to show ex/em lines
            if (is.null(event.data)) {
                
                #Distinguish if eem is plotted as contour or surface plot
                xinput <- min(ex)
                yinput <- min(em)
                
                em_x<- x[,match(xinput,ex)]
                
                ex_x<- x[match(yinput,em),]
                
                
                cont<-plot_ly(x=ex, y=em, z=x, type = "contour",contours=list(coloring ="heatmap"),ncontours=25,
                              colorscale = "Viridis",line=list(color="black", width=0.4),source = source) %>% 
                    layout(title = paste0( s_ind, "/", length(eem_list)," ", eem_list[[s_ind]]$sample[1]),
                           yaxis = list(showticklabels=T, showticks=F))
                
                
                
                emplot<-plot_ly( x=em_x , y= em, type="scatter", mode="lines", showlegend =F, line = list(color = "black"), visible= F,source = source) %>% 
                    layout(yaxis = list(title = "Em [nm]",range = c(min(em),max(em))), xaxis = list(side = "top", autorange = "reversed"))
                
                explot<-plot_ly( x= ex, y=ex_x , type="scatter", mode="lines", showlegend =F, line = list(color = "black"), visible= F,source = source)%>% 
                    layout(xaxis = list(title = "Ex [nm]",range = c(min (ex),max(ex))))
                
                dummy<-plot_ly(x=1, y=1, type="scatter", mode="markers", visible= F,source = source) %>% 
                    layout(yaxis = list(showticklabels=F, showgrid=F, showticks=F, zeroline= F),
                           xaxis = list(showticklabels=F, showgrid=F, showticks=F, zeroline= F))
                
                subplot(emplot, cont, dummy, explot, shareX = T, shareY = T, nrows = 2, heights = c(0.85,0.15), widths = c(0.15, 0.85), titleX = T, titleY = T)
                
                
                
                
            } else {
                #isolates x and y coordinate of click event
                xinput<-event.data$x
                yinput<-event.data$y
                #plot eem as contour + chosen ex/em lines
                em_x<- x[,match(xinput,ex)]
                
                ex_x<- x[match(yinput,em),]
                
                range.em.ind.min<-match(input$Emrange[1],em)
                range.em.ind.max<-match(input$Emrange[2],em)
                range.ex.ind.min<-match(input$Exrange[1],ex)
                range.ex.ind.max<-match(input$Exrange[2],ex)
                
                #cut original eem_list to setting applied on re() for comparison in the plot  
                #eem_list_red<- eem_list_main$eem_list
                
                # eem_cut<-cuteem(eem_list_red[[s_ind]],range.em.ind.min, range.em.ind.max,range.ex.ind.min,range.ex.ind.max)
                
                
                # em_old<-eem_cut$em #isolate emission of sample from unprocessed eem_list
                # ex_old<-eem_cut$ex #isolate excitation of sample from unprocessed eem_list
                # x_old<-eem_cut$x
                #isolate sample index
                # ex_x_old <- x_old[match(yinput,em_old),]
                #  em_x_old <- x_old[,match(xinput,ex_old)]
                
                #this plot is a bit unorganized and not easy to read, but it renders faster than creating 4 independent plots and combine them later with subplot
                
                subplot(
                    (plot_ly( x=em_x , y= em, type="scatter", mode="lines", showlegend =F, line = list(color = "black"), visible= T,source = source)%>% 
                         #add_trace(x=em_x_old, y=em_old, type = "scatter",mode="lines",
                         #         line = list(color = "grey", dash = 'dot'),showlegend = FALSE) %>% 
                         layout(yaxis = list(title = "Em [nm]",range = c(min(em),max(em))), xaxis = list(side = "top", autorange = "reversed", showticklabels = F))),
                    
                    (plot_ly(x=ex, y=em, z=x, type = "contour",contours=list(coloring ="heatmap"),ncontours=25,
                             colorscale = "Viridis",line=list(color="black", width=0.4),source = source) %>% 
                         layout(title = paste0( s_ind, "/", length(eem_list)," ", eem_list[[s_ind]]$sample[1]),
                                yaxis = list(showticklabels=T, showticks=F))%>% 
                         add_segments(x = xinput, xend = xinput, y = min(em), yend = max(em), inherit = FALSE,showlegend = FALSE,line = list(color = "black"))%>% 
                         add_segments(x = min(ex), xend = max(ex), y = yinput, yend = yinput, inherit = FALSE,showlegend = FALSE,line = list(color = "black"))),
                    
                    (plot_ly(x=1, y=1, type="scatter", mode="markers", visible= F,source = source) %>% 
                         layout(yaxis = list(showticklabels=F, showgrid=F, showticks=F, zeroline= F),
                                xaxis = list(showticklabels=F, showgrid=F, showticks=F, zeroline= F))),
                    
                    
                    (plot_ly( x= ex, y=ex_x , type="scatter", mode="lines", showlegend =F, line = list(color = "black"), visible= T,source = source)%>% 
                         #add_trace(x=ex_old, y=ex_x_old, type = "scatter",mode="lines",
                         #          line = list(color = "grey", dash = 'dot'),showlegend = FALSE) %>% 
                         layout(xaxis = list(title = "Ex [nm]",range = c(min (ex),max(ex))), yaxis = list (  showticklabels = F))),
                    
                    shareX = T, shareY = T, nrows = 2, heights = c(0.85,0.15), widths = c(0.15, 0.85), titleX = T, titleY = T)
                
                
                
            }
            
        }
        
        
        else {
            
            plot_ly(  x =  ex,
                      y =  em,
                      z = x,
                      type = "contour",contours=list(coloring ="heatmap"),ncontours=25,
                      colorscale = "Viridis",line=list(color="black", width=0.4),source = source) %>% 
                layout(title = paste0( s_ind, "/", length(eem_list)," ", eem_list[[s_ind]]$sample[1]),
                       yaxis = list ( title = "Em [nm]"), xaxis = list ( title = "[Ex [nm]"))
            
        }
        
        
        
    }
    ###
    SSE_ABCmode <- function(pf, eem_list){
        
        res<-lapply(pf, function(x){ residuals_sse (x, eem_list)}) %>% 
            lapply (function (x) {subset(x,x$type == "residual")})  %>% 
            lapply ( function (x) { 
                samples<-unique(x$Sample)
                ex<-unique(x$ex)
                em<-unique(x$em)
                
                sampleSSE<-sapply(samples, function(y){
                    sub<-subset(x,x$Sample == y)
                    sub<-sub$value^2
                    sum(sub, na.rm = T)
                })
                
                exSSE<-sapply(ex, function(y){
                    sub<-subset(x,x$ex == y)
                    sub<-sub$value^2
                    sum(sub, na.rm = T)
                })
                emSSE<-sapply(em, function(y){
                    sub<-subset(x,x$em == y)
                    sub<-sub$value^2
                    sum(sub, na.rm = T)
                })
                list(sampleSSE, exSSE, emSSE)
                
                
            })
        
        mypalette <- c(
            "blue", "orange", "green", "grey", "red", "yellow", "purple", "pink"
        )
        names<-paste0("C ",sapply(parafac(), function (x) {
            ncol(x$A)
        }))
        
        p_sample<- plot_ly ()
        for (i in 1:length(pf)){
            p_sample <- add_trace(p_sample, x = c(1:length( res[[1]][[1]])) , y = res[[i]][[1]], type = "scatter", mode = "lines", name = names[i], text= names(res[[i]][[1]]), showlegend = F,
                                  line = list(color = mypalette[i]))%>%
                layout (xaxis = list (title = "Sample", range = c(1,length( res[[1]][[1]]) )), yaxis = list ( title = "SSE"))
        }
        
        p_ex<- plot_ly ()
        
        for (i in 1:length(pf)){
            p_ex <- add_trace(p_ex, x = as.integer(names(res[[i]][[2]])) , y = res[[i]][[2]], type = "scatter", mode = "lines", name = names[i], text= names(res[[i]][[2]]), showlegend = F,
                              line = list(color = mypalette[i]))%>%
                layout (xaxis = list (title = "Excitation [nm]", range = c ( min(eem_list[[1]]$ex), max(eem_list[[1]]$ex))), yaxis = list ( title = "SSE"))
        }
        p_em<- plot_ly ()
        
        for (i in 1:length(pf)){
            p_em <- add_trace(p_em, x = as.integer(names(res[[i]][[3]])) , y = res[[i]][[3]], type = "scatter", mode = "lines", name = names[i], text= names(res[[i]][[3]]), showlegend = T ,
                              line = list(color = mypalette[i]))%>%
                layout (xaxis = list (title = "Emission [nm]", range = c ( min(eem_list[[1]]$em), max(eem_list[[1]]$em))), yaxis = list ( title = "SSE"))
        }
        
        
        subplot(p_sample, p_ex, p_em, nrows = 3, titleX = T, titleY = T, margin = 0.02)
    }
    ###
    sse_rsq_plot<-function(pf){
        nmodels<- sapply(pf, function (x) {
            ncol(x$A)
        })
        
        if (!is.null(pf[[1]]$SSE)){
            sse<- sapply( pf, function (x){
                x[["SSE"]]
            })
            
            Rsq<- sapply( pf, function (x){
                x[["Rsq"]]
            })
        }
        
        p1<- plot_ly(x=nmodels, y=Rsq, type="scatter", mode="lines + markers", name="Rsq")%>% 
            layout(title = "Model Fit and Residual Error", yaxis=list(title="Model fit (Rsq)",showgrid = F),
                   xaxis = list(title="Components", dtick = 1))
        
        p2<-plot_ly(x=nmodels, y=sse, type="scatter", mode="lines + markers", name="SSE")%>% 
            layout(title = "Model Fit and Residual Error", yaxis=list(title="Sum of squared Error",showgrid = F),
                   xaxis = list(title="Components", dtick = 1))
        
        subplot(p2, p1, nrows = 2, shareX = T, titleY = T)
        
        
        
    }
    ###
    parafac_contourplot<-   function(pf,ncompsmax){
        
        ncomps<-ncol(pf[["A"]])
        
        em<-as.numeric(rownames(pf[["B"]]))
        ex<-as.numeric(rownames(pf[["C"]]))
        nex<-length(ex)
        nem<-length(em)
        
        y<-eempf_comp_mat(pf)
        #this next condition decides if the layout should be adapted in case of multiple pf models.
        if (ncol(pf$A) == ncomps){
            
            plots<-lapply(y, function(x) {
                
                
                
                p<-plot_ly(
                    x =  ex,
                    y =  em,
                    z = matrix(data=x$value, ncol= nex, nrow=nem),
                    type = "contour",contours=list(coloring ="heatmap"),ncontours=25,
                    colorscale = "Viridis",line=list(color="black", width=0.4), showscale=F
                )
                p
            })
        }
        else{
            plots<-lapply(y, function(x) {
                
                
                
                p<-plot_ly(
                    x =  ex,
                    y =  em,
                    z = matrix(data=x$value, ncol= nex, nrow=nem),
                    type = "contour",contours=list(coloring ="heatmap"),ncontours=25,
                    colorscale = "Viridis",line=list(color="black", width=0.4), showscale=F
                ) %>% layout(yaxis = list(showticklabels=F))
                p
            }) 
            
        }
        
        
        
        
        #if multiple models are calculated, an empty plots is added to make all plots fit in nicely in subplot with 1 column for 1 model
        if(ncomps < ncompsmax){
            if (ncol(pf$A) == ncomps){
                
                for (i in 1:(ncompsmax-ncomps)) {
                    plots[[ncomps+i]]<-plot_ly(type="scatter", mode="markers")%>% 
                        layout(yaxis = 
                                   list(showgrid = T, range= c(
                                       min(ex),max(ex)
                                   )),
                               xaxis=
                                   list(showgrid = T))
                    
                    #   x =  ex,
                    #    y =  em,
                    #    z = rep(0, (nex*nem)),
                    #    type = "contour",contours=list(coloring ="heatmap"),ncontours=25,
                    #    colorscale = "Viridis",line=list(color="black", width=0.4), showscale=F
                    #  )
                    
                }
            }
            else {
                for (i in 1:(ncompsmax-ncomps)) {
                    plots[[ncomps+i]]<-plot_ly(type="scatter", mode="markers")%>% 
                        layout(yaxis = 
                                   list(showticklabels=F,showgrid = T, range= c(
                                       min(ex),max(ex)
                                   )),
                               xaxis=
                                   list(showgrid = T))
                    
                }  
                
                
            }
        }
        
        subplot(plots, nrows=ncompsmax, shareX = T, shareY = T, margin = 0.005)
        
        
        
        
        
    }
    ###
    residuals_shiny<- function(pf, eem_list, output_plot = F, SSEonly = F){
        ex = as.integer(eem_list[[1]]$ex)
        nex = length(eem_list[[1]]$ex)
        em = as.integer(eem_list[[1]]$em)
        nem = length(eem_list[[1]]$em)
        
        model_mat<- eempf_comp_mat(pf) %>% lapply (function(x){x$value}) %>% as.data.frame()
        pf<-norm2A(pf)
        Amodes<- as.data.frame(t(pf[["A"]]))
        
        data<-model_mat
        model<- lapply(Amodes, function(Amodes){
            d<-lapply(1:length(Amodes), function (i){
                data[,i]<-Amodes[i]*model_mat[,i]
            }) %>% as.data.frame() %>% rowSums() 
            d<-matrix(data= d, nrow= nem, ncol =nex)
            return(d)
        })
        
        sample <- lapply(eem_list, function ( x){ x$x})
        
        res<- sample
        residual <- lapply(1:length(eem_list), function(i) {
            res[[i]] <- sample[[i]] - model[[i]] 
        }) 
        
        names(residual) <- names( model)
        names(sample) <- names( model)
        
        res<-list( "Sample" = sample,"Model" = model, "Residual" =residual)
        if (output_plot){
            samplenames <- eem_names(eem_list)
            nsamples<- length(eem_list)
            p_sample<- lapply(1:nsamples, function(i){
                
                samplename<- samplenames[i]
                
                plot_ly (x = ex, y = em, z = sample[[i]],
                         type = "contour",contours=list(coloring ="heatmap"),ncontours=25,
                         colorscale = "Viridis",line=list(color="black", width=0.4)) %>% 
                    layout(title = paste0( i, "/", nsamples," ", samplename),
                           yaxis = list(showticklabels=T, showticks=F, title = " Em [nm]"),
                           xaxis = list ( title = "Ex [nm]"))
            })
            p_model<- lapply(1:nsamples, function(i){
                
                samplename<- samplenames[i]
                
                plot_ly (x = ex, y = em, z = model[[i]],
                         type = "contour",contours=list(coloring ="heatmap"),ncontours=25,
                         colorscale = "Viridis",line=list(color="black", width=0.4)) %>% 
                    layout(title = paste0( i, "/", nsamples," ", samplename),
                           yaxis = list(showticklabels=T, showticks=F, title = " Em [nm]"),
                           xaxis = list ( title = "Ex [nm]"))
            })
            
            p_residual<- lapply(1:nsamples, function(i){
                
                samplename<- samplenames[i]
                
                plot_ly (x = ex, y = em, z = residual[[i]],
                         type = "contour",contours=list(coloring ="heatmap"),ncontours=40,
                         colorscale = "Viridis",line=list(color="black", width=0.4)) %>% 
                    layout(title = paste0( i, "/", nsamples," ", samplename),
                           yaxis = list(showticklabels=T, showticks=F, title = " Em [nm]"),
                           xaxis = list ( title = "Ex [nm]"))
            })
            
            res<-list( "Sample" = p_sample,"Model" = p_model, "Residual" =p_residual)
        }
        if (SSEonly){
            res<-sapply(residual, function (x){
                res<-x^2 %>%sum (na.rm = T)
            })
        }
        return (res)
        
    }
    ###
    contourplot_pf_single<- function( pf){
        model_mat<- eempf_comp_mat(pf) %>% lapply (function(x){x$value}) %>% as.data.frame()
        ncomps<- ncol(model_mat)
        ex<-as.integer(rownames(pf$C))
        nex <- length(ex)
        em<-as.integer(rownames(pf$B))
        nem <- length(em)
        
        if (ncomps>4){
            nrows <- 2
        }
        else{
            nrows <- 1
        }
        
        
        
        plot<- lapply(1:ncomps, function(i){
            
            plot_ly (x = ex, y = em,
                     z = matrix(data = model_mat[,i], ncol = nex, nrow = nem),
                     type= "contour",contours=list(coloring ="heatmap"),ncontours=25,
                     colorscale = "Viridis",line=list(color="black", width=0.4), showscale = FALSE ) %>%
                layout ( yaxis = list ( title = "Em [nm]"), xaxis = list ( title = "Ex [nm]"))
            
        }) %>% subplot (shareY = T, shareX = T, titleX = T, titleY = T, margin = 0.005, nrows = nrows)
        
    }
    ###
    
    lineplot_pf_single <- function (pf){
        ncomps<- ncol(pf$A)
        ex<-as.integer(rownames(pf$C))
        nex <- length(ex)
        em<-as.integer(rownames(pf$B))
        nem <- length(em)
        
        
        if (ncomps>4){
            nrows <- 2
        }
        else{
            nrows <- 1
        }
        
        plot<- lapply(1:ncomps, function(i){
            plot_ly (x = ex, y = pf$C[,i],type="scatter", mode="lines + markers", name = "Ex",line = list(color = "blue"), showlegend = F)%>%
                add_trace (p, x = em, y = pf$B[,i], name = "Em",line = list(color = "black"), showlegend = F) %>%
                layout (yaxis = list (showgrid = F), xaxis = list ( showgrid = T, title = "Wavelength [nm]"))
        }) %>% subplot (shareY = T, shareX = T, titleX = T, titleY = T, margin = 0.005, nrows = nrows)
        
    }
    
    

    #attempt to work with an eem_list as reactive value. does not work yet because ui observeing events are dependent on it
    eem_list_main <- reactiveValues( eem_list= list(), watchlist = list (), temp = list(), manualscatter = 0, scattervalue = data.frame (), setNA = list(),setNAcheck = F,
                                     residuals = list (),samplenames_res = c(), corplot = list (), levplot = list(), loadingsplot = list (), pfcontourplot = list (), pfspectrallines = list(),
                                     models = list (), modelnames = c(), modelncomps = list(), modelnameexport = c(), pfcontourplotsingle = list (), pflineplotsingle = list ()
                                     )
    
    removed_samples <- reactiveValues( samples= c(), ex1= c(), ex2= c(), em1= c(), em2= c(), undo= c(), eem_list_ex = list(), state= FALSE)
    
    
    
    ##########counter
    #counter for next/previous samples
    counter <- reactiveValues(countervalue = 0) # Defining & initializing the reactiveValues object
    
    observeEvent(input$add1, {
        counter$countervalue <- counter$countervalue + 1     # if the add button is clicked, increment the value by 1 and update it
        updateSelectInput(session, "sample", choices = 
                              setNames(
                                  c(1:length(eem_list_main$eem_list)),
                                  paste0("sample",c(1:length(eem_list_main$eem_list))," ",sapply( 
                                      eem_list_main$eem_list, function(x){x$sample})
                                  )), selected = 
                              setNames(
                                  c(1:length(eem_list_main$eem_list)),
                                  paste0("sample",c(1:length(eem_list_main$eem_list))," ",sapply( 
                                      eem_list_main$eem_list, function(x){x$sample})
                                  ))[as.integer(input$sample[1])+counter$countervalue]
        )
        
    })
    observeEvent(input$sub1, {
        counter$countervalue <- counter$countervalue - 1  # if the sub button is clicked, decrement the value by 1 and update it
        updateSelectInput(session, "sample", choices = 
                              setNames(
                                  c(1:length(eem_list_main$eem_list)),
                                  paste0("sample",c(1:length(eem_list_main$eem_list))," ",sapply( 
                                      eem_list_main$eem_list, function(x){x$sample})
                                  )), selected = 
                              setNames(
                                  c(1:length(eem_list_main$eem_list)),
                                  paste0("sample",c(1:length(eem_list_main$eem_list))," ",sapply( 
                                      eem_list_main$eem_list, function(x){x$sample})
                                  ))[as.integer(input$sample[1])+counter$countervalue]
        )
    })
    
    observeEvent(input$sample, { # switches counter to zero whenever dropdown menu is activated
        counter$countervalue <- 0
        
    })
    
    
    
    
    ####################################################################################################################################################
    uploadcheck<- reactiveValues(eems = F, spectralcorex = F,spectralcorem = F, absorbance = F, meta = F, stateoverview = NULL, state = 
                                     c("You can either upload raw datafiles or a previously corrected EEM List.", "If you choose raw datafiles and want to perform any of the Corrections,","a Metatable is needed containing sample, blank and optional absorbance names in first/second/third column.", "Progress:"))
    
    observeEvent(input$eemraw,{
        uploadcheck$state <- c(uploadcheck$state,paste0(length(eem_list_unprocessed()), " EEMs uploaded"))
        uploadcheck$stateoverview <- summary(eem_list_unprocessed())
    })
    observeEvent(input$spectralcorrectionex,{
        uploadcheck$state <- c(uploadcheck$state,"Spectral Correction File for Excitation uploaded")
        uploadcheck$stateoverview <- excor()
        updateCheckboxInput(session, "spectralcor", value = T)
    })
    observeEvent(input$spectralcorrectionem,{
        uploadcheck$state <- c(uploadcheck$state,"Spectral Correction File for Emission uploaded")
        uploadcheck$stateoverview <- emcor()
        updateCheckboxInput(session, "spectralcor", value = T)
    })
    observeEvent(input$absorbanceraw,{
        uploadcheck$state <- c(uploadcheck$state,paste0((ncol(abs())-1), " Absorbance Files uploaded"))
        uploadcheck$stateoverview <- abs()
        updateCheckboxInput(session, "ifecor", value = T)
    })
    observeEvent(input$meta,{
        uploadcheck$state <- c(uploadcheck$state,"Metatable uploaded")
        uploadcheck$stateoverview <- meta()
    })
    observeEvent(input$applycorrection,{
        uploadcheck$stateoverview <- summary(eem_list())
    })
    observeEvent(input$eeminput,{
        uploadcheck$state <- c(uploadcheck$state,"EEM List uploaded")
    })
    observeEvent(input$uploadfile,{
        uploadcheck$stateoverview <- summary(eem_list())
    })
    observeEvent(input$download_dreem_button,{
        req(eem_list())
        uploadcheck$stateoverview <- summary(eem_list())
        uploadcheck$state <- c(uploadcheck$state,"Dr.eem Dataset downloaded")
    })
    
    
    
    #options(shiny.maxRequestSize=30*1024^2)
    
    eem_list_unprocessed <- eventReactive(input$eemraw, {
        req(input$eemraw)
        showNotification("EEMs uploaded. Processing...")
        tryCatch({
            eem_list = list()
            
            for(nr in 1:length(input$eemraw[, 1])){
                eem_list[[nr]] <- eem_read(
                    file = input$eemraw[[nr, 'datapath']],
                    import_function = input$spectrophotometer
                )[[1]]
                eem_list[[nr]]$sample <- substr(input$eemraw$name[[nr]], 1, nchar(input$eemraw$name[[nr]])-4) 
            }
            class(eem_list)<- "eemlist"
            uploadcheck$eems <- T
            
            
            
            return(eem_list)
        },error=function(cond){uploadcheck$state <- c(uploadcheck$state,"EEM upload failed. Are the files in the right format?")})
    })
    
    
    excor <- eventReactive(input$spectralcorrectionex, {
        req(input$spectralcorrectionex)
        showNotification("File uploaded. Processing...")
        excor<-as.data.frame(read.csv (file = input$spectralcorrectionex[['datapath']]))
        excor[,1]<- round(excor[,1])
        
        uploadcheck$spectralcorex <- T
        
        message ("Spectral Correction File for Excitation uploaded")
        
        
        
        return(excor)
    })
    
    emcor <- eventReactive(input$spectralcorrectionem, {
        req(input$spectralcorrectionem)
        showNotification("File uploaded. Processing...")
        emcor<-as.data.frame(read.csv (file = input$spectralcorrectionem[['datapath']]))
        emcor[,1]<- round(emcor[,1])
        
        uploadcheck$spectralcorem <- T
        
        message ("Spectral Correction File for Emission uploaded")
        
        return(emcor)
    })
    #################
    abs <- eventReactive(input$absorbanceraw, {
        req(input$absorbanceraw)
        showNotification("Absorbance files uploaded. Processing...")
        abs_data <- c()
        abs_names <- c()
        for(i in 1:length(input$absorbanceraw[, 1])){
            
            abs_data <- c( abs_data, input$absorbanceraw[[i, 'datapath']])
            abs_names<- c(abs_names,substr(input$absorbanceraw$name[[i]], 1, nchar(input$absorbanceraw$name[[i]])-4)  )
        }
        
        abs_data<-read_absorbance_shiny(abs_data)
        
        
        colnames(abs_data)<-c("wavelength", abs_names)
        
        uploadcheck$absorbance <- T
        
        message ( paste0((ncol(abs_data)-1), " Absorbance Files uploaded"))
        
        
        return(abs_data)
    })
    
    ##########################
    meta <- eventReactive(input$meta, {
        req(input$meta)
        showNotification("File uploaded. Processing...")
        meta<-as.data.frame(read.csv (file =input$meta[['datapath']]))
        
        uploadcheck$meta <- T
        
        message ("Metadata File uploaded")
        
        
        return(meta)
    })
    
    
    
    
    eem_list<-eventReactive(c(input$applycorrection,input$uploadfile, input$download_dreem_button), {
        if (input$selectupload == "uploadraw"){
            
            
            meta<-meta()
            all_samples<-sapply(eem_list_unprocessed(), function(x){x$sample})
            meta.list <- split(meta, seq(nrow(meta)))
            
            eem_list<-eem_red2smallest(eem_list_unprocessed(), verbose = FALSE) %>% lapply( function(x) {
                
                x$ex <- round(x$ex)
                x$em <- round(x$em)
                x
            })
            class(eem_list) <- "eemlist"
            
            
            if (input$spectralcor){
                eem_list <- eem_spectral_cor(eem_list,excor(),emcor())
            }
            # blanks<-unique(meta[,2])
            
            blanks <- lapply(unique(meta[,2]), function (x){
                blank <-eem_list[[match(x, all_samples)]]
                return (blank)
            })
            class(blanks)<- "eemlist"
            
            
            if (input$blankcor){
                
                eem_list<-lapply(meta.list, function(x){
                    sample_ind<-match(x[[1]], all_samples)
                    blank_ind<-match(x[[2]], all_samples)
                    sample<-eem_list[[sample_ind]]
                    blank<-list(eem_list[[blank_ind]])
                    class(blank)<- class(eem_list)
                    res<-eem_remove_blank(sample, blank=blank)
                    res$sample <- x[[1]]
                    res
                })
                class(eem_list)<- "eemlist"
                eem_list <- eem_bind(blanks, eem_list)
            }
            
            if (input$ifecor){
                abs_data<-abs()
                #switch names to match eem_list
                oldnames<- colnames(abs_data)[-1]
                
                newnames<-c("wavelength",
                            sapply(oldnames, function(x){
                                abs_ind<- match(x, meta[,3])
                                meta[,1][abs_ind]
                            }))
                colnames(abs_data) <- newnames
                
                if (input$abs_blcor){
                    abs_data<- abs_blcor(abs_data,wlrange = c((max(abs_data[,1])-20),max(abs_data[,1])))
                }
                
                #set ranges for absorbance and eem data
                eem_list <- eem_list %>% eem_range(ex = c(min(abs_data[,1]),Inf), em = c(min(abs_data[,1]),Inf)) %>% eem_ife_correction(abs_data, cuvl = 1, unit= "absorbance")
                
            }
            
            
            
            
            if (input$ramannorm){
                all_samples<-sapply(eem_list, function(x){x$sample})
                
                eem_list<-lapply(meta.list,function(x){
                    sample_ind<-match(x[[1]], all_samples)
                    blank_ind<-match(x[[2]], all_samples)
                    sample<-eem_list[[sample_ind]]
                    blank<-list(eem_list[[blank_ind]])
                    class(blank)<- "eemlist"
                    class(sample)<- "eem"
                    res<-eem_raman_normalisation(sample, blank=blank)
                    res$sample <- x[[1]]
                    res
                    
                })
                class(eem_list)<- "eemlist"
                
                
            }
            
            eem_list_main$eem_list<- eem_list
            #this way watchlist has the appropriate data structure
            eem_list_main$watchlist <- eem_list
            eem_list_main$watchlist<- eem_list_main$watchlist[-(1:length(eem_list))]
            eem_list_main$state<- T
            
            updateSelectInput(session, "sample", choices = 
                                  setNames(
                                      c(1:length(eem_list_main$eem_list)),
                                      paste0("sample",c(1:length(eem_list_main$eem_list))," ",sapply( 
                                          eem_list_main$eem_list, function(x){x$sample})
                                      ))
            )
            
            
            
            return(eem_list)
        }
        else if (input$selectupload == "upload_eemlist"){
            if (!is.null(input$eeminput)){
                showNotification("EEMlist uploaded. Processing...")
                fileupload<-load(input$eeminput$datapath)
                fileupload<-eval(as.symbol(fileupload))
                file<-lapply(fileupload, function(x) {
                    
                    x$ex <- round(x$ex)
                    x$em <- round(x$em)
                    x
                })
                attributes(file)<-attributes(fileupload)
                eem_list_main$eem_list<- file
                #this way watchlist has the appropriate data structure
                eem_list_main$watchlist <- file
                eem_list_main$watchlist<- eem_list_main$watchlist[-(1:length(file))]
                eem_list_main$state<- T
                
                
                updateSelectInput(session, "sample", choices = 
                                      setNames(
                                          c(1:length(eem_list_main$eem_list)),
                                          paste0("sample",c(1:length(eem_list_main$eem_list))," ",sapply( 
                                              eem_list_main$eem_list, function(x){x$sample})
                                          ))
                )
                
                
                return(file)
            }
        }
            else if (input$selectupload == "download_dreem"){
                eem_list <- eem_load_dreem()

                
                eem_list_main$eem_list<- eem_list
                #this way watchlist has the appropriate data structure
                eem_list_main$watchlist <- eem_list
                eem_list_main$watchlist<- eem_list_main$watchlist[-(1:length(eem_list))]
                eem_list_main$state<- T

                return(eem_list)            
                }
            else{
                #Dummy list to prevent app from crashing if actionbutton is clicked without uploading file
                #should be changed to trycatch at some point
             
                return(list(list(sample =1, ex=1, em=1)))
            }
            
        
    })
    
    
    
    
    
    ########################################################################################################################################
    
    
    
    #######
    
    
    
    observe({
        updateSliderInput(session, "Exrange", 
                          min = min(sapply(eem_list(), function(x){min(x$ex)})),
                          max = max(sapply(eem_list(), function(x){max(x$ex)})),
                          step=(eem_list()[[1]]$ex[2]-eem_list()[[1]]$ex[1]),
                          value = c(min(eem_list()[[1]]$ex), max(eem_list()[[1]]$ex))
        )
        updateSliderInput(session, "Emrange", 
                          min = min(sapply(eem_list(), function(x){min(x$em)})),
                          max = max(sapply(eem_list(), function(x){max(x$em)})),
                          step=(eem_list()[[1]]$em[2]-eem_list()[[1]]$em[1]),
                          value = c(min(eem_list()[[1]]$em), max(eem_list()[[1]]$em))
        )
    })
    

    
    
    # observe if a sample should be removed from dataset
    observeEvent(input$remove, {
        req(re())
        s_ind<-as.integer(input$sample[1])
        
        removed_samples$samples<-c(removed_samples$samples,eem_list_main$eem_list[[s_ind]]$sample)
        
        
        eem_list_main$eem_list[[s_ind]]<-NULL
        updateSelectInput(session, "sample", choices = 
                              setNames(
                                  c(1:length(eem_list_main$eem_list)),
                                  paste0("sample",c(1:length(eem_list_main$eem_list))," ",sapply( 
                                      eem_list_main$eem_list, function(x){x$sample}))
                              ), selected = 
                              setNames(
                                  c(1:length(eem_list_main$eem_list)),
                                  paste0("sample",c(1:length(eem_list_main$eem_list))," ",sapply( 
                                      eem_list_main$eem_list, function(x){x$sample})
                                  ))[as.integer(input$sample[1])+counter$countervalue]
        )
        
        showNotification(
            paste0("Removed Sample: ",removed_samples$samples[length(removed_samples$samples)])
        )
        #message(c("Removed Sample: ",removed_samples$samples[length(removed_samples$samples)]))
        
        updateSelectInput(session,"watchlist_samples", choices=removed_samples$samples)
        
        
    })
    
    ###
    observeEvent(input$removemodel, {
        req(parafac())
        showModal(modalDialog(
            tagList(
                selectInput("deletemodelname", label = "Delete a Model", choices = eem_list_main$modelnames,
                            selected = eem_list_main$modelnames[[match(input$choose_pf_investigation, eem_list_main$modelnames)]] )
            ), 
            title="Delete a Model",
            footer = tagList(actionButton("confirmDelete", "Delete"),
                             modalButton("Cancel")
            )
        ))
    })
    
    
    
    
    #######
    # observe if a model should be removed from dataset
    observeEvent(input$confirmDelete, {
        
        modelindex <- match(input$deletemodelname, eem_list_main$modelnames)
        
        if (length(eem_list_main$modelnames) >1){
        showNotification(
            paste0("Deleted Model: ",input$deletemodelname)
        )
        
        
        eem_list_main$models[[modelindex]] <- NULL
        eem_list_main$modelnames[[modelindex]] <- NULL
        eem_list_main$modelncomps[[modelindex]] <- NULL
        eem_list_main$corplot[[modelindex]] <- NULL
        eem_list_main$levplot[[modelindex]] <- NULL
        eem_list_main$loadingsplot[[modelindex]] <- NULL
        eem_list_main$pfcontourplotsingle[[modelindex]] <- NULL
        eem_list_main$pflineplotsingle[[modelindex]] <- NULL

        
        updateSelectInput(session, "choose_pf_res", choices = eem_list_main$modelnames,
                          selected = eem_list_main$modelnames[[length(eem_list_main$modelnames)]]
                          
        )
        updateSelectInput(session, "choose_pf_investigation", choices = eem_list_main$modelnames,
                          selected = eem_list_main$modelnames[[length(eem_list_main$modelnames)]]
                          
        )
        updateSelectInput(session, "choose_pf_export", choices = eem_list_main$modelnames,
                          selected = eem_list_main$modelnames[[length(eem_list_main$modelnames)]]
                          
        )
        updateSelectInput(session, "choose_pf_sh", choices = eem_list_main$modelnames,
                          selected = eem_list_main$modelnames[[length(eem_list_main$modelnames)]]
                          
        )
        }
        else{
            showNotification(
                "You only have one Model!"
            )   
        }

        
        
        #message(c("Removed Sample: ",removed_samples$samples[length(removed_samples$samples)]))
        

        
    })
    ################Watchlist
    
    
    observeEvent(input$moveback_watchlist, {
        samplename<-input$watchlist_samples
        removed_samples_all<-removed_samples$samples
        ind<-match(samplename,removed_samples_all)
        if(!is.na(ind)){
            removed_samples_all<- removed_samples_all[-ind]
            if (is.null(removed_samples_all)){
                removed_samples$samples <- NULL 
            }
            else{
                removed_samples$samples <- removed_samples_all
                
            }
            
            
            removed_samples_list_names<- sapply(eem_list_main$watchlist, function(x) {x$sample})
            
            
            eem_list_main$eem_list[[(length(eem_list_main$eem_list)+1)]] <- eem_list_main$watchlist[[match(samplename, removed_samples_list_names)]]
            eem_list_main$watchlist[[match(samplename, removed_samples_list_names)]] <- NULL
            updateSelectInput(session,"watchlist_samples", choices=removed_samples$samples)
            updateSelectInput(session, "sample", choices = 
                                  setNames(
                                      c(1:length(eem_list_main$eem_list)),
                                      paste0("sample",c(1:length(eem_list_main$eem_list))," ",sapply( 
                                          eem_list_main$eem_list, function(x){x$sample}))
                                  ), selected = 
                                  setNames(
                                      c(1:length(eem_list_main$eem_list)),
                                      paste0("sample",c(1:length(eem_list_main$eem_list))," ",sapply( 
                                          eem_list_main$eem_list, function(x){x$sample})
                                      ))[as.integer(input$sample[1])+counter$countervalue]
            )
            
            
        }
        
        showNotification(
            paste0("Sample Moved to main Dataset: ",samplename)
        )
        #message(c("Sample Moved to main Dataset: ", samplename))
    })
    
    observeEvent(input$apply, {
        #samplename<-substr(names(input$sample), 1, 7)
        samplename <- re()[[as.integer(input$sample)]]$sample
        
        
        
        ex1<- as.integer(input$ex1)
        ex2<- as.integer(input$ex2)
        em1<- as.integer(input$em1)
        em2<- as.integer(input$em2)
        
        
        if (!is.na(sum(ex1,ex2,em1,em2))){
            ex<-ex1:ex2
            em<-em1:em2
            
            eem_list_main$setNA[[(1+ length(eem_list_main$setNA))]]<-list(samplename, ex, em)
            
            eem_list_main$setNAcheck <- T
            
        }
        
    })
    
    observeEvent(input$undo, {
        if (length(eem_list_main$setNA) == 1){
            eem_list_main$setNAcheck <- F
            
        }
        if (length(eem_list_main$setNA) > 0){
            eem_list_main$setNA[[length(eem_list_main$setNA)]]<-NULL
        }
        
        
    })
    
    
    
    ########
    
    
    re <- reactive({ 
        req(eem_list())
        req(eem_list_main$state)
        
        ############
        ray1 <- input$ray1
        ray2 <- input$ray2
        raman1 <- input$raman1
        raman2 <- input$raman2
        
        em<-eem_list()[[1]]$em
        ex<-eem_list()[[1]]$ex
        
        range.em.ind.min<-match(input$Emrange[1],em)
        range.em.ind.max<-match(input$Emrange[2],em)
        range.ex.ind.min<-match(input$Exrange[1],ex)
        range.ex.ind.max<-match(input$Exrange[2],ex)
        
        data.cor<- eem_list() %>% lapply( function(x) {
            x$scatter <- data.frame(ray1, ray2, raman1, raman2)
            x
        }) %>% lapply(remove_scatter) %>% lapply(cuteem,range.em.ind.min, range.em.ind.max,range.ex.ind.min,range.ex.ind.max)
        
        #set na
        
        if (eem_list_main$setNAcheck == T){
            for (i in 1:length(eem_list_main$setNA)){
                data.cor <- eem_setNA(data.cor, sample = eem_list_main$setNA[[i]][[1]], ex =eem_list_main$setNA[[i]][[2]], em =eem_list_main$setNA[[i]][[3]], interpolate = F)
            }
        }
        
        
        ###########
        
        
        if (input$interpolate==T){
            
            data.cor<- eem_interp(data=data.cor, type=input$interpolation.type[1])
            
        }
        
        attributes(data.cor)<-attributes(eem_list())
        
        
        if (!is.null(removed_samples$samples) & length(removed_samples$samples) > 0){
            
            for (i in 1:length(removed_samples$samples)){
                
                eem_list_main$watchlist[[i]]<-data.cor[[match(removed_samples$samples[i], sapply(data.cor, function(x){x$sample}))]]
                
                
            }
            data.cor <- eem_exclude(data.cor, list (ex = c(), em = c(), sample = removed_samples$samples ))
            eem_list_main$eem_list <- data.cor
            
            
        }
        
        return(data.cor)
    }) 
    
    
    output$uploadstate <-renderPrint({
        if (!is.null (  uploadcheck$state)){
            for (i in 1:length(uploadcheck$state)){
                print (  uploadcheck$state[i])
            }
            
        }
    })
    
    output$uploadsummary <- DT::renderDataTable({
        tryCatch({
            if(!is.null (  uploadcheck$stateoverview )){
                uploadcheck$stateoverview
            }
        },error=function(cond){message ("error")})
        
    })
    
    
    output$EEMplot <- renderPlotly({
        req(re())
        
        event.data <- event_data(event = "plotly_click", source = "imgLink") #isolate event (click)
        
        #isolate sample index and make sure only items on the less will be indiced
        s_ind<-as.integer(input$sample[1])+counter$countervalue
        
        if (s_ind > length(re())){
            s_ind <- length( re())
            
            counter$countervalue <- counter$countervalue -1
        }
        if (s_ind < 1){
            s_ind <- 1
            counter$countervalue <- counter$countervalue +1
        }
        
        eem_plot_contour_lines(re(), s_ind, event.data, input$includelines,source = "imgLink")
        
    })
    
    
    
    
    
    output$spectralvariance <- renderPlotly({
        
        if (input$checkbox[1]) {
            
            spectralvarianceplot(re())
        }
    })
    
    
    ##
    output$EEMplotwatchlist <- renderPlotly({
        req(is.character(removed_samples$samples))
        
        #isolate sample index
        s_ind<-
            match(input$watchlist_samples,removed_samples$samples)
        if (!is.na(s_ind)){
            event.data <- event_data(event = "plotly_click", source = "watchlist") #isolate event (click)
            
            eem_plot_contour_lines(re(), s_ind, event.data, input$includelines, source = "watchlist")
            
        }
        
        
        
        
    })
    
    
    
    
    
    
    #################
    # output$click <- renderPrint({
    #    d <- event_data("plotly_click", source="imgLink")
    #    if (is.null(d)){ }
    #    else {
    ##      updateTextInput(session,"ex1", value = 
    #                      eem_list()[[1]]$ex[which.min(abs(eem_list()[[1]]$ex - d$x[1]))])
    #updateTextInput(session,"ex2", value =  eem_list()[[1]]$ex[which.min(abs(eem_list()[[1]]$ex - d$x[1]))])
    #      updateTextInput(session,"em1", value =  eem_list()[[1]]$em[which.min(abs(eem_list()[[1]]$em - d$y[1]))])
    #updateTextInput(session,"em2", value =  eem_list()[[1]]$em[which.min(abs(eem_list()[[1]]$em - d$y[1]))])
    
    # }
    #print("Drag on plot to retrieve Ex and EM values")
    #print("Still a bit buggy, dragging and zooming sometimes switch")
    #     print( " ")
    #
    #  })
    #   observe({
    #      updateTextInput(session,"ex1", value =  eem_list()[[1]]$ex[which.min(abs(eem_list()[[1]]$ex - (event_data("plotly_click", source="imgLink"))$x[1]))])
    #      updateTextInput(session,"ex2", value =  eem_list()[[1]]$ex[which.min(abs(eem_list()[[1]]$ex - (event_data("plotly_click", source="imgLink"))$x[1]))])
    #      updateTextInput(session,"em1", value =  eem_list()[[1]]$em[which.min(abs(eem_list()[[1]]$em - (event_data("plotly_click", source="imgLink"))$y[1]))])
    #      updateTextInput(session,"em2", value =  eem_list()[[1]]$em[which.min(abs(eem_list()[[1]]$em - (event_data("plotly_click", source="imgLink"))$y[1]))])
    #    })
    
    ###############
    
    ###################page 2 event for parafac analysis
    
    parafac<-eventReactive(input$applyparafac,{
        
        showNotification("Calculating PARAFAC Model...")
        
        #isolate eem ex/em, values
        ctol <- as.numeric(input$ctol) # decrease tolerance in PARAFAC analysis
        nstart <- as.numeric(input$nstart) # increase number of random starts
        maxit <- as.numeric(input$maxit)# increase number of maximum interations
        pf <- tryCatch(
            {   
                eem_parafac(re() , comps = as.integer(input$nmodels), normalise = input$normalize, 
                            const = c(input$constraints, input$constraints, input$constraints), maxit = maxit, nstart = nstart, ctol = ctol,strictly_converging = input$str_conv,verbose=T)  
            },
            error=function(cond){
                message("Error. Increase nstarts/maxit or decrease tolerance. May have chosen too many Components")
            }
        )
        updateCheckboxInput(session, "reversenorm", value =input$normalize )
        updateCheckboxInput(session, "normalize_sh", value =input$normalize )
        updateCheckboxInput(session, "str_conv_sh", value =input$str_conv )
        updateTextInput(session, "ctol_sh", value =input$ctol )
        updateTextInput(session, "nstart_sh", value =input$nstart )
        updateTextInput(session, "maxit_sh", value =input$maxit )
        updateSelectInput(session, "constraints_sh", choices =c("Non-Negativity" ="nonneg",
                                                                "Unconstrained" = "uncons",
                                                                "Smoothed Non-Negativity" = "smonon"), selected = input$constraints)
        ################
        #Residuals
        ###
        
        
        
        samplenames_res<-setNames(
            c(1:length(re())),
            paste0("sample",c(1:length(re()))," ",sapply( 
                re(), function(x){x$sample})
            ))
        
        
        
        
        eem_list_main$models[[length(eem_list_main$models)+1]]<-pf
        if (nchar(input$modelname) >0){
            modelname <- input$modelname
        }
        else if (nchar(input$modelname) ==0 & input$normalize){
            modelname <-paste0("pf",(length(eem_list_main$modelnames)+1), "norm" )
        }
        else{
            modelname <-paste0("pf",(length(eem_list_main$modelnames)+1))  
        }
        
        eem_list_main$modelnames[[length(eem_list_main$modelnames)+1]]<-modelname
    
        eem_list_main$modelncomps[[length(eem_list_main$modelncomps)+1]] <-as.integer(input$nmodels)
        updateSelectInput(session, "sample_res", choices = samplenames_res
                          
        )
        updateSelectInput(session, "choose_pf_res", choices = eem_list_main$modelnames,
                          selected = modelname
                          
        )
        updateSelectInput(session, "choose_pf_investigation", choices = eem_list_main$modelnames,
                          selected = modelname
                          
        )
        updateSelectInput(session, "choose_pf_export", choices = eem_list_main$modelnames,
                          selected = modelname
                          
        )
        updateSelectInput(session, "choose_pf_sh", choices = eem_list_main$modelnames,
                          selected = modelname
                          
        )
        ###
        showNotification("Calculating Residuals...")
        eem_list_main$residuals[[length(eem_list_main$residuals)+1]] <- lapply(pf,residuals_shiny, eem_list =re(), output_plot = T)

        
        eem_list_main$samplenames_res <- samplenames_res
        
        #further calculate multiple plots for PARAFAC Model determination
        showNotification("Calculating Plots...")
        eem_list_main$corplot[[length(eem_list_main$corplot)+1]] <- lapply(pf,function(x){
            p<-eempf_corplot(x, normalisation = input$normalize)  %>% suppressWarnings() %>% ggplotly() %>% suppressWarnings()
            p
        }
        )
        eem_list_main$levplot[[length(eem_list_main$levplot)+1]] <- lapply(pf, function(x){
            p<-eempf_leverage_plot(eempf_leverage(x),qlabel=0.1) %>% ggplotly () %>% suppressWarnings()
            p
        })
        
        eem_list_main$loadingsplot[[length(eem_list_main$loadingsplot)+1]]<-lapply(pf, function(x){
            p<-eempf_comp_load_plot(x)
            p[[2]] %>% ggplotly () %>% suppressWarnings()
        })
        
        eem_list_main$pfcontourplotsingle[[length(eem_list_main$pfcontourplotsingle)+1]] <-lapply(pf, eempf_rescaleBC, newscale = "Fmax") %>% lapply(contourplot_pf_single)
        
        eem_list_main$pflineplotsingle[[length(eem_list_main$pflineplotsingle)+1]] <-lapply(pf, eempf_rescaleBC, newscale = "Fmax") %>% lapply(lineplot_pf_single)

        #####calculate contour and line plot of all models
        
        ncompsmax<-max(
            sapply(pf, function(x){
                ncol(x[["A"]])
            })
        )
        pf_newscale <- lapply(pf, eempf_rescaleBC, newscale = "Fmax")
        
        #determine maximum number of compounds

        
        eem_list_main$pfcontourplot <-lapply(pf, eempf_rescaleBC, newscale = "Fmax") %>% lapply(parafac_contourplot, ncompsmax = ncompsmax) %>% subplot()
        
        
        eem_list_main$pfspectrallines<-ggplotly(eempf_compare(pf_newscale)[[3]])
        
        assign("alldata_contour", eem_list_main$pfcontourplotsingle, envir = .GlobalEnv)
        assign("alldata_lev", eem_list_main$levplot, envir = .GlobalEnv)
        
        return(pf)
    })
    
  #  res <- reactive({
  #      req(parafac())
  #      lapply(parafac (),residuals_shiny, eem_list = re(), output_plot = T)
  #  })
    ##########################
    #create dynamiclly updated overview of available models
    
    observe({
        ####################################################
        
            req(nchar(input$choose_pf_res) != 1)
        updateRadioButtons(session, "choose_model", choices = 
                               eem_list_main$modelncomps[[match(input$choose_pf_investigation, eem_list_main$modelnames)]]#,
                           #selected = min(eem_list_main$modelncomps[[match(input$choose_pf_investigation, eem_list_main$modelnames)]])
        )
        
        updateRadioButtons(session, "choose_model2", choices = 
                               eem_list_main$modelncomps[[match(input$choose_pf_res, eem_list_main$modelnames)]]#,
                           #selected = min(eem_list_main$modelncomps[[match(input$choose_pf_investigation, eem_list_main$modelnames)]])                          
                           )
        updateRadioButtons(session, "choose_model3", choices = 
                               eem_list_main$modelncomps[[match(input$choose_pf_sh, eem_list_main$modelnames)]]#,
                           #selected = min(eem_list_main$modelncomps[[match(input$choose_pf_investigation, eem_list_main$modelnames)]])                          
        )
        updateRadioButtons(session, "choose_model4", choices = 
                               eem_list_main$modelncomps[[match(input$choose_pf_export, eem_list_main$modelnames)]]#,
                           #selected = min(eem_list_main$modelncomps[[match(input$choose_pf_investigation, eem_list_main$modelnames)]])                          
        )
    
        
    })


    
    #################################page 2
    #######################
    output$pfmodel <- renderPlotly({ 
        req(parafac())
        
        
        if (input$switchpfplot == "contour"){
            eem_list_main$pfcontourplot[[1]]
        }
        
        else {
            #this case means spectrallines are chosen. function from Stardom.
            
            eem_list_main$pfspectrallines[[1]]
        }
        
        
    })
    
    ###########################
    output$pfdiagnosis <- renderPlotly({
        
        p<-tryCatch({
            sse_rsq_plot(parafac())
            
            
        },
        error=function(cond){
        })
        p
        
    })
    
    output$SSEplot<- renderPlotly({ 
        
        SSE_ABCmode ( parafac(), re())
        
        
        
        
        
    })
    
    ################################
    
    ##########################
    #event for splithalf
    splithalf_re<-eventReactive(input$applysh,{
        req(re())
        showNotification("Calculating Splithalf...")
        
        ctol <- as.numeric(input$ctol_sh) # decrease tolerance in PARAFAC analysis
        nstart <- as.numeric(input$nstart_sh) # increase number of random starts
        maxit <- as.numeric(input$maxit_sh)# increase number of maximum interations
        sh <- splithalf(re() , comps = as.integer(input$choose_model3), normalise = input$normalize_sh, rand = input$random_sh,
                        const = c(input$constraints_sh, input$constraints_sh, input$constraints_sh), maxit = maxit, nstart = nstart, ctol = ctol,strictly_converging = input$str_conv_sh,verbose=T)
        return(sh)
    })
    
    output$shplot<-renderPlot({
        
        
        plot_sh(splithalf_re())
        
    })
    output$ssc_sh <-renderTable({
        
        splithalf_tcc(splithalf_re())
        
    })
    
    observe({
        req(splithalf_re())
        
        tcc_sh_table <- splithalf_tcc(splithalf_re())
        
        if (any(c(tcc_sh_table$tcc_ex,tcc_sh_table$tcc_em) <0.95)){
            text<-"Model not validated!"
        }
        else {
            text<-"Model validated!"
            
        }
        
        showModal(modalDialog(
            title = "Splifhalfvalidation",
            text
        ))
    })
    
    
    # model<-eventReactive(input$apply_model_invenstigation,{
    #    pf<-parafac()
    #    nmodels<-sapply(pf, function(x){ncol(x$A)})
    #    modelindex<-match(input$choose_model,nmodels)
    #    pf[[modelindex]]
    
    # })
    
    #  model2<-eventReactive(input$apply_model_invenstigation2,{
    #    pf<-parafac()
    #    nmodels<-sapply(pf, function(x){ncol(x$A)})
    #    modelindex<-match(input$choose_model2,nmodels)
    #    modelindex
    
    #  })
    
    
    ##########counter
    #counter for next/previous samples in residual plots
    counter_res <- reactiveValues(countervalue = 0) # Defining & initializing the reactiveValues object
    
    observeEvent(input$add1_res, {
        counter_res$countervalue <- counter_res$countervalue + 1     # if the add button is clicked, increment the value by 1 and update it
        updateSelectInput(session, "sample_res", choices = eem_list_main$samplenames_res, selected = 
                              eem_list_main$samplenames_res[as.integer(input$sample_res[1])+counter_res$countervalue]
        )
        
    })
    observeEvent(input$sub1_res, {
        counter_res$countervalue <- counter_res$countervalue - 1  # if the sub button is clicked, decrement the value by 1 and update it
        updateSelectInput(session, "sample_res", choices = eem_list_main$samplenames_res, selected = 
                              eem_list_main$samplenames_res[as.integer(input$sample_res[1])+counter_res$countervalue]
        )
    })
    
    observeEvent(input$sample_res, { # switches counter to zero whenever dropdown menu is activated
        counter_res$countervalue <- 0
        
    })
    
    
    output$corrplot <- renderPlotly({
        req(parafac())
        modelindex <- match(input$choose_pf_investigation, eem_list_main$modelnames)
        
        nmodels<-eem_list_main$modelncomps[[modelindex]]
        pfindex<-match(input$choose_model,nmodels)

        if (input$choose_plot == "corr"){
            eem_list_main$corplot[[modelindex]][[pfindex]]
            
        }
        else if (input$choose_plot == "lev"){
            eem_list_main$levplot[[modelindex]][[pfindex]]
            
        }
        else if (input$choose_plot == "loadings"){
            eem_list_main$loadingsplot[[modelindex]][[pfindex]]
            
        }
        else if (input$choose_plot == "contour"){
            eem_list_main$pfcontourplotsingle[[modelindex]][[pfindex]]
        }
        else if (input$choose_plot == "lines"){
            eem_list_main$pflineplotsingle[[modelindex]][[pfindex]]
        }
    })
    
    
    
    output$resplot <- renderPlotly({
        req(parafac())
        modelindex <- match(input$choose_pf_res, eem_list_main$modelnames)
        nmodels<-eem_list_main$modelncomps[[modelindex]]
        
        pfindex<-match(input$choose_model2,nmodels)
        sampleindex <- as.integer(input$sample_res[1])+counter_res$countervalue
        

        
        
        if (sampleindex > 0 & sampleindex< (length(eem_list_main$residuals[[modelindex]][[pfindex]][["Sample"]])+1)){
            
            if (input$choose_plot_res == "res_all"){
                
                
                
                subplot(eem_list_main$residuals[[modelindex]][[pfindex]][["Sample"]][[sampleindex]],
                        eem_list_main$residuals[[modelindex]][[pfindex]][["Model"]][[sampleindex]],
                        eem_list_main$residuals[[modelindex]][[pfindex]][["Residual"]][[sampleindex]],
                        shareY = T)
            }
            else {
                eem_list_main$residuals[[modelindex]][[pfindex]][["Residual"]][[sampleindex]]
            }
            
        }
        
        
        
    })



    output$downloadeem_list <- downloadHandler(
        filename = function() {
            "eem_list.RData"
        },
        content = function(file) {
            eem_list <- re()
            save(eem_list, file=file)
        }
    )
    
    
    output$downloadmodel <- downloadHandler(
        filename = function() {
            "model.csv"
        },
        content = function(file) {
            modelindex <- match(input$choose_pf_export, eem_list_main$modelnames)
            
            nmodels<-sapply(parafac(), function(x){ncol(x$A)})
            pfindex<-match(input$choose_model4,nmodels)
            model<-eem_list_main$models[[modelindex]][[pfindex]] 
            
            if (input$reversenorm){
                model<-norm2A(model)
            }
            
            
            export<-staRdom::eempf_export (model)
            write.csv(export, file)
        }
    )
    
    output$downloadopenfluor <- downloadHandler(
        filename = function() {
            "model_openfluor.txt"
        },
        content = function(file) {
            modelindex <- match(input$choose_pf_export, eem_list_main$modelnames)
            
            nmodels<-sapply(parafac(), function(x){ncol(x$A)})
            pfindex<-match(input$choose_model4,nmodels)
            model<-eem_list_main$models[[modelindex]][[pfindex]] %>% eempf_rescaleBC(newscale = "Fmax")
            
            if (input$reversenorm){
                model<-norm2A(model)
            }
            
            
            
            #export<-staRdom::eempf_openfluor(model,file.path(tempdir()) )
            staRdom::eempf_openfluor(model, file)
        }
    )
  #  
    
    
}
