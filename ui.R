#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(plotly)


shinyUI(
    navbarPage("EEM-PARAFAC",
               tabPanel("File Upload",
                        sidebarPanel(
                            radioButtons("selectupload", "Select file input", choices= c("Raw Data" ="uploadraw",
                                                                                         "Corrected EEM List" = "upload_eemlist",
                                                                                         "Use dr.eem dataset" = "download_dreem")),
                            conditionalPanel(
                                condition = "input.selectupload == 'uploadraw'",
                                radioButtons("selectrawmethod", "Upload Files", choices= c("Upload EEMs" = "uploadeemsraw",
                                                                                           "Spectral Correction Files" ="uploadspectralcor",
                                                                                           "Upload Absorbance Files" = "uploadifecor",
                                                                                           "Metadata" = "uploadmetaraw")),
                                conditionalPanel(
                                    condition = "input.selectrawmethod == 'uploadeemsraw'",
                                    selectInput("spectrophotometer", label = "Spectrophotometer", choices = 
                                                    c("cary", "aqualog", "shimadzu", "fluoromax4")),
                                    fileInput(
                                        inputId = "eemraw", 
                                        label = "Upload EEMs", 
                                        multiple = TRUE,
                                        accept = c("text/csv",
                                                   "text/comma-separated-values,text/plain",
                                                   ".csv")
                                    )
                                ),
                                conditionalPanel(
                                    condition = "input.selectrawmethod == 'uploadspectralcor'",
                                    fileInput(
                                        inputId = "spectralcorrectionex", 
                                        label = "Spectral Correction for Excitation", 
                                        multiple = F,
                                        accept = c("text/csv",
                                                   "text/comma-separated-values,text/plain",
                                                   ".csv")
                                    ),
                                    fileInput(
                                        inputId = "spectralcorrectionem", 
                                        label = "Spectral Correction for Emission", 
                                        multiple = F,
                                        accept = c("text/csv",
                                                   "text/comma-separated-values,text/plain",
                                                   ".csv")
                                    )
                                ),
                                conditionalPanel(
                                    condition = "input.selectrawmethod == 'uploadifecor'",
                                    
                                    checkboxInput("abs_blcor", label = "Perform baseline Correction", value = T),
                                    fileInput(
                                        inputId = "absorbanceraw", 
                                        label = "Absorbance Data (May take a minute to process)", 
                                        multiple = T,
                                        accept = c("text/csv",
                                                   "text/comma-separated-values,text/plain",
                                                   ".csv")
                                    )
                                    
                                ),
                                conditionalPanel(
                                    condition = "input.selectrawmethod == 'uploadmetaraw'",
                                    fileInput(
                                        inputId = "meta", 
                                        label = "Upload Metatable", 
                                        multiple = F,
                                        accept = c("text/csv",
                                                   "text/comma-separated-values,text/plain",
                                                   ".csv")
                                    )
                                ),
                                
                                
                                checkboxInput("blankcor", label = "Blank Correction", value = T),
                                checkboxInput("spectralcor", label = "Spectral Correction", value = F),
                                checkboxInput("ifecor", label = "IFE Correction", value = F),
                                checkboxInput("ramannorm", label = "Raman Normalisation", value = T),
                                
                                actionButton("applycorrection", label = "Apply Correction")
                                
                            ),
                            
                            
                            conditionalPanel(
                                condition = "input.selectupload == 'upload_eemlist'",
                                fileInput("eeminput", "Upload formatted eem_list", multiple=F,
                                          accept = ".RData"),
                                actionButton("uploadfile", "Load in Data")
                            ),
                            conditionalPanel(
                              condition = "input.selectupload == 'download_dreem'",
                        
                              actionButton("download_dreem_button", "Load in Data")
                            )
                            
                        ),
                        mainPanel(
                            verbatimTextOutput("uploadstate"),
                            dataTableOutput("uploadsummary")
                            
                        )
               ),
               tabPanel("Processing",
                        sidebarPanel(
                            selectInput("sample", "Sample:", choices=c(1:1)),
                            actionButton("sub1", "Previous"),
                            actionButton("add1", "Next"),
                            checkboxInput("includelines", label = "Show Ex and Em Lines by clicking (slower)"),
                            actionButton("remove","Remove this sample"),
                            
                            checkboxInput("checkbox", label = "Show spectralvariance", value = FALSE),
                            radioButtons("selectprocessing", "Select Method", choices= c("Change Range" ="processrange",
                                                                                         "Change Scatter Removal" = "processscatter",
                                                                                         "Cut Noise in specified Sample" = "processnoise",
                                                                                         "Interpolate Missing Numbers" = "processinterpolate")),
                            conditionalPanel(
                                condition = "input.selectprocessing == 'processrange'",
                                sliderInput("Exrange", label = h5("Excitation Range [nm] "), min = 1, 
                                            max = 2,step=1, value = c(1,2)),
                                sliderInput("Emrange", label = h5("Emission Range [nm] "), min = 1, 
                                            max = 2,step=1, value = c(1, 2))
                            ),
                            
                            conditionalPanel(
                                condition = "input.selectprocessing == 'processscatter'",
                                sliderInput("ray1", label = h5("Rayleigh 1 Scatter (below and above)"), min = -30, 
                                            max = 30,step=1, value = c(-10,10)),
                                sliderInput("raman1", label = h5("Raman 1 Scatter (below and above)"), min = -30, 
                                            max = 30,step=1, value = c(-10,10)),
                                sliderInput("ray2", label = h5("Rayleigh 2 Scatter (below and above)"), min = -30, 
                                            max = 30,step=1, value = c(-10,10)),
                                sliderInput("raman2", label = h5("Raman 2 Scatter (below and above)"), min = -30, 
                                            max = 30,step=1, value = c(-10,10))
                                
                            ),
                            
                            conditionalPanel(
                                condition = "input.selectprocessing == 'processnoise'",
                                textInput("ex1", "From Ex wavelength" ),
                                textInput("ex2", "To Ex wavelength" ),
                                textInput("em1", "From Em wavelength" ),
                                textInput("em2", "To Em wavelength" ),
                                actionButton("apply", "Apply to sample"),
                                actionButton("undo", "Undo")
                                
                                
                            ),
                            
                            conditionalPanel(
                                condition = "input.selectprocessing == 'processinterpolate'",
                                checkboxInput("interpolate", label = "Interpolate scatter. Careful, large computation", value = FALSE),
                                selectInput("interpolation.type", label = h5("Interpolation Method"),
                                            choices=c(0:4), selected= 1)
                            )
                            
                            
                            
                            
                            
                            
                            
                        ),
                        
                        mainPanel(
                            plotlyOutput("EEMplot",height = "700px"),
                            verbatimTextOutput("click"),
                            plotlyOutput("spectralvariance")
                        )
                        
               ),
               tabPanel("Watchlist",
                        sidebarPanel(
                            selectInput("watchlist_samples", "Removed Samples:", choices=c("None")),
                            checkboxInput("includelineswatchlist", label = "Show Ex and Em Lines by clicking (slower)"),
                            actionButton("moveback_watchlist", label= "Move selected sample to main Dataset")
                        ),
                        mainPanel(
                            
                            plotlyOutput("EEMplotwatchlist")
                            
                            
                            
                            
                        )
               ),
               tabPanel("PARAFAC",
                        sidebarPanel(
                            checkboxGroupInput("nmodels", "Number of Components:",
                                               c("2 Components" = 2,
                                                 "3 Components" = 3,
                                                 "4 Components" = 4,
                                                 "5 Components" = 5,
                                                 "6 Components" = 6,
                                                 "7 Components" = 7,
                                                 "8 Components" = 8,
                                                 "9 Components" = 9)),
                            textInput("nstart", "n of random starts", value= 30),
                            textInput("maxit", "maximum of interations", value= 2500),
                            textInput("ctol", "Tolerance", value= 10^-8),
                            checkboxInput("normalize", label = "Normalize Intensity", value = FALSE),
                            checkboxInput("str_conv", label = "Strictly Converging", value = FALSE),
                            selectInput("constraints", "Constraints:", choices=c("Non-Negativity" ="nonneg",
                                                                                 "Unconstrained" = "uncons",
                                                                                 "Smoothed Non-Negativity" = "smonon"
                            )),
                            textInput("modelname", "Modelname"),
                            actionButton("applyparafac", "Apply PARAFAC"),
                            radioButtons("switchpfplot", "PARAFAC plot", choices= c("Contour plot" ="contour",
                                                                                    "Spectral Lines" = "spectralllines"))
                            
                        ),
                        
                        mainPanel(
                            plotlyOutput("pfmodel", height = "800px"),
                            plotlyOutput("pfdiagnosis"),
                            plotlyOutput("SSEplot")
                            
                        )
                        
               ),
               tabPanel("Further Model Investigation",
                        sidebarPanel(
                            selectInput("choose_pf_investigation", "Choose Model", choices = c(1)),
                            radioButtons("choose_model", "Choose Component", choices= c("Calculate Model First!")),
                            radioButtons("choose_plot", "Choose plot", choices= 
                                             c("Contour plot" = "contour",
                                               "Lineplot" = "lines",
                                             "Component Correlation" = "corr",
                                               "Leverage" = "lev",
                                               "Component Loadings" = "loadings")),
                            actionButton("removemodel","Delete a Model"),
                            
                        ),
                        mainPanel(
                            plotlyOutput("corrplot",width = "800px", height = "400px")
                            
                        )
               ),
               tabPanel("Residuals",
                        sidebarPanel(
                            selectInput("choose_pf_res", "Choose Model", choices = c(1)),
                            radioButtons("choose_model2", "Choose Component", choices= c("Calculate Model First!")),
                            selectInput("sample_res", "Sample:", choices=c(1:1)),
                            actionButton("sub1_res", "Previous"),
                            actionButton("add1_res", "Next"),
                            radioButtons("choose_plot_res", "Choose plot", choices= 
                                             c(
                                                 "Residuals and Model" = "res_all",
                                                 "Residuals only" = "res_only"))
                        ),
                        mainPanel(
                            
                            plotlyOutput("resplot")
                            
                        )
               ),
               
               tabPanel("Splithalf Analysis",
                        sidebarPanel(
                            selectInput("choose_pf_sh", "Choose Model", choices = c(1)),
                            radioButtons("choose_model3", "Choose Component Number", choices= c("Calculate Model First!")),
                            textInput("nstart_sh", "Random Starts", value= 30),
                            textInput("maxit_sh", "Maximum  Iterations", value= 2500),
                            textInput("ctol_sh", "Tolerance", value= 10^-8),
                            checkboxInput("normalize_sh", label = "Normalize Intensity", value = FALSE),
                            checkboxInput("random_sh", label = "Assign Splits Randomly", value = FALSE),
                            checkboxInput("str_conv_sh", label = "Strictly Converging", value = FALSE),
                            
                            selectInput("constraints_sh", "Constraints:", choices=c("Non-Negativity" ="nonneg",
                                                                                    "Unconstrained" = "uncons",
                                                                                    "Smoothed and Non-Negativity" = "smonon")),
                            
                            actionButton("applysh", "Calculate Splithalf")
                        ),
                        mainPanel(
                            plotOutput("shplot"),
                            tableOutput("ssc_sh"),

                        )),
               tabPanel("Model Export",
                        sidebarPanel(
                            
                            selectInput("choose_pf_export", "Choose Model", choices = c(1)),
                            radioButtons("choose_model4", "Choose Component Number", choices= c("Calculate Model First!")),
                            #textInput("modelname", "Modelname", value = "Mymodel"),
                            checkboxInput("reversenorm", "Reverse normalization", value = F),
                            
                            downloadButton("downloadeem_list", 'Download EEMlist'),
                            downloadButton("downloadmodel", 'Download Model'),
                            downloadButton("downloadopenfluor", 'Download Model OpenFluor')
                            
                            
                            
                        ),
                        mainPanel(
                            plotOutput("test")
                            
                            
                        ))
    )
)

