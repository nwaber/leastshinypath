library(shiny)
library(sf)
library(terra)
library(leastcostpath)

create_conductance_matrix <- function(x) {
  cells <- which(!is.na(terra::values(x)))
  na_cells <- which(is.na(terra::values(x)))
  adj <- terra::adjacent(x = x, cells = cells, pairs = TRUE)
  adj <- adj[!adj[, 2] %in% na_cells, ]
  
  spatvals <- terra::values(x)[, 1]
  spatvals <- spatvals[adj[, 2]]
  
  ncells <- length(cells) + length(na_cells)
  cs_matrix <- Matrix::Matrix(data = 0, nrow = ncells, ncol = ncells)
  cs_matrix[adj] <- spatvals
  
  cs <- list(conductanceMatrix = cs_matrix, 
             costFunction = NA, 
             "maxSlope" = NA,
             exaggeration = FALSE,
             criticalSlope = NA,
             neighbours = NA, 
             nrow = terra::nrow(x), 
             ncol = terra::ncol(x), 
             "resolution" = terra::res(x), 
             "extent" = as.vector(terra::ext(x)),  
             crs = terra::crs(x, proj = TRUE))
  
  class(cs) <- "conductanceMatrix"
  
  return(cs)
}


ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(
        tabPanel("Cost Surface",
          div(
            h3("Use Your Own Cost Surface"),
            fileInput("userCostSurfaceFile", "Upload User Cost Surface (GeoTIFF)")
          ),
          div(
            h3("Calculate Cost Surface"),
            fileInput("demFile", "Upload DEM file (GeoTIFF)"),
            selectInput("costFunction", "Select Cost Function", choices = c(
              "tobler", "tobler offpath", "davey", "rees", "irmischer-clarke male", 
              "irmischer-clarke offpath male", "irmischer-clarke female", "irmischer-clarke offpath female",
              "modified tobler", "garmy", "kondo-saino", "wheeled transport", "herzog", 
              "llobera-sluckin", "naismith", "minetti", "campbell", "campbell 2019 1",
              "campbell 2019 5", "campbell 2019 10", "campbell 2019 15", "campbell 2019 20", 
              "campbell 2019 25", "campbell 2019 30", "campbell 2019 35", "campbell 2019 40", 
              "campbell 2019 45", "campbell 2019 50", "campbell 2019 55", "campbell 2019 60", 
              "campbell 2019 65", "campbell 2019 70", "campbell 2019 75", "campbell 2019 80", 
              "campbell 2019 85", "campbell 2019 90", "campbell 2019 95", "campbell 2019 99", 
              "sullivan 167", "sullivan 5", "sullivan 833")),
            selectInput("neighbours", "Select Neighbours", choices = c(4, 8, 16, 32)),
            numericInput("critSlope", "Critical Slope", value = 12),
            numericInput("maxSlope", "Max Slope", value = NULL),
            checkboxInput("exaggeration", "Exaggeration", value = FALSE)
          ),
          div(
          actionButton("calculateCostSurface", "Calculate Cost Surface")
          ),
          div(
          downloadButton("downloadData", "Download Cost Surface")
          ),
          div(
            h3("Cost Surface Stochasticity"),
            actionButton("addStochasticity", "Add Global Stochasticity")
          ),
          div(
            h3("Download Raster"),
            downloadButton("downloadRaster2", "Download Displayed Raster")
          )
        ), #close tabPanel
        tabPanel("Least cost path",
           fileInput("pointFile", "Upload Point File (CSV format)"),
           uiOutput("idMenu"),
           uiOutput("xMenu"),
           uiOutput("yMenu"),
           actionButton("setCoords", "Set Coordinates"),
           actionButton("showData", "Show Points Table"),
           actionButton("plotPoints", "Plot Points"),
           div(
           selectInput("lcpMethod", "Least cost path method", 
                       choices = c("Sequential", "From one to everywhere", "From everywhere to everywhere")
                       ),
           conditionalPanel(
             condition = "input.lcpMethod == 'Sequential'",
             uiOutput("sequenceFieldMenu")
           ),
           conditionalPanel(
             condition = "input.lcpMethod == 'From one to everywhere'",
             uiOutput("idValueMenu")
           ),
           selectInput("costSurface", "Select Cost Surface", 
                       choices = NULL),
           actionButton("calculatePaths", "Calculate Paths"),
           actionButton("plotPaths", "Plot Paths"),
           downloadButton("downloadPaths", "Download Paths")
           ),
           div(
            h3("Movement Corridors"),
            actionButton("calculateCorridors", "Calculate Path Corridors")
           ),
           div(
             h3("Path Density"),
             actionButton("createDensity", "Calculate Path Density"),
             downloadButton("downloadDensity", "Download Path Density Map")
           ),
           div(
             h3("Download Raster"),
             downloadButton("downloadRaster", "Download Displayed Raster")
           )
        ), #close tabPanel
        selectInput("rasterData", "Display Raster Map:", choices = NULL)
      ) #close tabsetPanel

    ), #close sidebarPanel
    
    mainPanel(
      fluidRow(
          plotOutput("demPlot", height = "100vh", width = "100%"),
          outputId = "rasterMap",
            height = "100vh" ,
            width = "100%"
      ) #close fluidRow
    ) #close mainPanel
  ) #close sidebarLayout
) #close UI fluidPage


server <- function(input, output, session) {
  options(shiny.maxRequestSize = 30*1024^2)
  # Initialize an empty list to store the raster data
  raster_data <- reactiveValues(data = list())
  values <- reactiveValues()
  cc <- list()
  
  user_cost_surface <- reactive({
    req(input$userCostSurfaceFile)
    inFile <- input$userCostSurfaceFile
    if (is.null(inFile)) return(NULL)
    terra::rast(inFile$datapath)
  })
  # observe({
  #   ucs <- user_cost_surface()
  #   if (!is.null(ucs)) {
  #     raster_data$data[["User Cost Surface"]] <- ucs
  #     updateSelectInput(session, "rasterData", choices = names(raster_data$data))
  #   }
  # })
  # When a new User Cost Surface file is uploaded, add it to the list and update the selectInput choices
  observe({
    ucs <- user_cost_surface()
    if (!is.null(ucs)) {
      raster_data$data[["User Cost Surface"]] <- ucs
      updateSelectInput(session, "rasterData", choices = names(raster_data$data))
    }
  })
  
  
  
  uploaded_dem <- reactive({
    req(input$demFile)
    inFile <- input$demFile
    if (is.null(inFile)) return(NULL)
    terra::rast(inFile$datapath)
  })
  
  # When a new DEM file is uploaded, add it to the list and update the selectInput choices
  observe({
    dem <- uploaded_dem()
    if (!is.null(dem)) {
      raster_data$data$DEM <- dem
      updateSelectInput(session, "rasterData", choices = names(raster_data$data))
    }
  })
  
  points_data <- reactive({
    req(input$pointFile)
    inFile <- input$pointFile
    if (is.null(inFile)) return(NULL)
    read.csv(inFile$datapath, stringsAsFactors = FALSE)
  })
  output$idMenu <- renderUI({
    req(points_data())
    selectInput("idField", "ID", choices = names(points_data()))
  })
  
  output$xMenu <- renderUI({
    req(points_data())
    selectInput("xField", "X Coordinate Field", choices = names(points_data()))
  })
  
  output$yMenu <- renderUI({
    req(points_data())
    selectInput("yField", "Y Coordinate Field", choices = names(points_data()))
  })
  
  observeEvent(input$setCoords, {
    req(input$idField, input$xField, input$yField, points_data())
    
    points <- points_data()
    
    points_sf <- sf::st_sf(
      geometry = sf::st_sfc(
        lapply(1:nrow(points), function(i) {
          sf::st_point(c(points[i, input$xField], points[i, input$yField]))
        })
      ),
      crs = terra::crs(raster_data$data[[input$rasterData]])
    )
    
    # Store the spatial points in a reactiveValues object for use elsewhere in the app
    values$points_sf <- points_sf
  })
  
  
  observeEvent(input$showData, {
    showModal(modalDialog(
      title = "Points Data",
      tableOutput("pointsTable")
    ))
  })
  output$pointsTable <- renderTable({
    points_data()
  })
  
  
  

  
  selected_cost_function <- reactive({
    cost_function_name <- input$costFunction
    cost_function(cost_function_name)
  })
  
  output$demPlot <- renderPlot({
    dem <- uploaded_dem()
    if (!is.null(dem)) {
      terra::plot(dem, main = "Uploaded DEM")
    }
  })
  
  
  observeEvent(input$calculateCostSurface, {
    cost_function <- input$costFunction
    neighbours <- as.numeric(input$neighbours)
    crit_slope <- input$critSlope
    max_slope <- input$maxSlope
    exaggeration <- input$exaggeration
    
    # Calculate the cost surface
    slope_cs <- create_slope_cs(
      x = uploaded_dem(), 
      cost_function = cost_function, 
      neighbours = neighbours,
      crit_slope = crit_slope,
      max_slope = max_slope,
      exaggeration = exaggeration)
    slope_cs_rast <- leastcostpath::rasterise(slope_cs)
    
    observe({
      if (!is.null(slope_cs)) {
        showNotification("Cost surface calculated")
      }
    })
    
    if (!is.null(slope_cs_rast)) {
      raster_data$data[["Cost Surface"]] <- slope_cs_rast
      updateSelectInput(session, "rasterData", choices = names(raster_data$data))
    }
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("cost-surface_", gsub(" ", "-", input$costFunction), "_", Sys.Date(), ".tif", sep="")
    },
    content = function(file) {
      terra::writeRaster(
        raster_data$data[["Cost Surface"]],
        file)
    }
  )
  
  # observeEvent(input$addStochasticity, {
  #   showModal(modalDialog(
  #     title = "Add Global Stochasticity",
  #     numericInput("percentQuantile", "Percent Quantile", value = 0.5, min = 0.00, max = 1.00),
  #     actionButton("goButton", "Go")
  #   ))
  # })
  # observeEvent(input$goButton, {
  #   req(input$percentQuantile)
  #   stochastic_cs <- leastcostpath::add_global_stochasticity(x = slope_cs, percent_quantile = input$percentQuantile)
  #   raster_data$data[["Stochastic Cost Surface"]] <- rasterise(stochastic_cs)
  #   updateSelectInput(session, "rasterData", choices = names(raster_data$data))
  # 
  #   removeModal()
  # })
  
  observeEvent(input$addStochasticity, {
    showModal(modalDialog(
      title = "Add Global Stochasticity",
      numericInput("percentQuantile", "Percent Quantile", value = 0.5, min = 0.00, max = 1.00),
      selectInput("costSurfaceSource", "Select Cost Surface Source", choices = c("Generated Cost Surface", "User Cost Surface")),
      actionButton("goButton", "Go")
    ))
  })
  observeEvent(input$goButton, {
    req(input$percentQuantile)
    
    # Determine which cost surface to use
    if (input$costSurfaceSource == "Generated Cost Surface") {
      cs <- raster_data$data[["Cost Surface"]]
    } else if (input$costSurfaceSource == "User Cost Surface") {
      cs <- create_conductance_matrix(raster_data$data[["User Cost Surface"]])
      cs$crs <- terra::crs(raster_data$data[["User Cost Surface"]])
    }
    
    stochastic_cs <- leastcostpath::add_global_stochasticity(x = cs, percent_quantile = input$percentQuantile)
    raster_data$data[["Stochastic Cost Surface"]] <- rasterise(stochastic_cs)
    updateSelectInput(session, "rasterData", choices = names(raster_data$data))
    
    removeModal()
  })
  
  
  output$demPlot <- renderPlot({
    selected_data <- input$rasterData
    if (!is.null(selected_data)) {
      terra::plot(raster_data$data[[selected_data]], main = selected_data)
    }
  })
  
  # Observe the raster_data$data list for changes
  observe({
    # Get the names of the rasters
    raster_names <- names(raster_data$data)
    
    # Filter the names to include only those that contain the word "Cost"
    cost_raster_names <- raster_names[grep("Cost", raster_names)]
    
    # Update the selectInput choices
    updateSelectInput(session, "costSurface", choices = cost_raster_names)
  })
  
  
  # observeEvent(input$plotPoints, {
  #   output$demPlot <- renderPlot({
  #     selected_data <- input$rasterData
  #     if (!is.null(selected_data)) {
  #       terra::plot(raster_data$data[[selected_data]], main = selected_data)
  #       
  #       if (!is.null(values$points_sf)) {
  #         plot(values$points_sf,
  #              add = TRUE,
  #              pch = 16)
  #         points <- points_data()
  #         coords <- sf::st_coordinates(values$points_sf)
  #         text(coords[,1], coords[,2], labels = points[[input$idField]], pos = 4, cex = 0.8)
  #         
  #       }
  #     }
  #   })
  # })
  
  observeEvent(input$plotPoints, {
    output$demPlot <- renderPlot({
      selected_data <- input$rasterData
      if (!is.null(selected_data)) {
        rast <- raster_data$data[[selected_data]]
        if (!is.null(rast)) {
          if (is.null(terra::crs(rast))) {
            showModal(modalDialog(
              title = "Enter EPSG Code",
              numericInput("epsgCode", "EPSG Code:", value = NULL),
              actionButton("submitEPSG", "Submit")
            ))
          } else {
            terra::plot(rast, main = selected_data)
            if (!is.null(values$points_sf)) {
              plot(values$points_sf,
                   add = TRUE,
                   pch = 16)
              points <- points_data()
              coords <- sf::st_coordinates(values$points_sf)
              text(coords[,1], coords[,2], labels = points[[input$idField]], pos = 4, cex = 0.8)
            }
          }
        }
      }
    })
  })
  
  observeEvent(input$submitEPSG, {
    req(input$epsgCode)
    selected_data <- input$rasterData
    rast <- raster_data$data[[selected_data]]
    terra::crs(rast) <- as.character(input$epsgCode)
    removeModal()
  })
  
  
  output$sequenceFieldMenu <- renderUI({
    req(points_data())
    selectInput("sequenceField", "Sequence field", choices = names(points_data()))
  })
  
  output$idValueMenu <- renderUI({
    req(points_data(), input$idField)
    selectInput("idValue", "Origin", choices = unique(points_data()[[input$idField]]))
  })
  
  observeEvent(input$calculatePaths, {
    req(values$points_sf, input$lcpMethod)
    
    if (input$costSurface == "Cost Surface") {
      cs <- raster_data$data[["Cost Surface"]]
    } else if (input$costSurface == "User Cost Surface") {
      cs <- create_conductance_matrix(raster_data$data[["User Cost Surface"]])
      cs$crs <- terra::crs(raster_data$data[["User Cost Surface"]])
    } else if (input$costSurface == "Stochastic Cost Surface") {
      cs <- create_conductance_matrix(raster_data$data[["Stochastic Cost Surface"]])
      cs$crs <- terra::crs(raster_data$data[["Stochastic Cost Surface"]])
    }
    
    if (input$lcpMethod == "Sequential") {
      points <- points_data()[order(points_data()[[input$sequenceField]]), ]
      #points_data()[[input$idField]] <- as.numeric(points_data()[[input$idField]])
      values$points_sf <- values$points_sf[order(points_data()[[input$sequenceField]]), ]
      
      lcps <- lapply(1:(nrow(points) - 1), function(i) {
        create_lcp(x = cs, origin = values$points_sf[i,], destination = values$points_sf[i + 1,])
      })
      
      lcps <- do.call(rbind, lcps)       # Combine all the least cost paths into one object
      
    } else if (input$lcpMethod == "From one to everywhere") {
      origin <- values$points_sf[which(points_data()[[input$idField]] == input$idValue), ]
      lcps <- create_lcp(x = cs, origin = origin, destination = values$points_sf)
      
    } else if (input$lcpMethod == "From everywhere to everywhere") {
      lcps <- create_FETE_lcps(x = cs, locations = values$points_sf)
    }
    
    values$lcps <- lcps
    
    observe(
      if(!is.null(lcps)){
    showNotification("Least cost paths calculated successfully!")
      })
  })
  
  observeEvent(input$plotPaths, {
    output$demPlot <- renderPlot({
      selected_data <- input$rasterData
      if (!is.null(selected_data)) {
        # Plot the selected raster data
        terra::plot(raster_data$data[[selected_data]], main = selected_data)
        
        # If points have been set, plot them on top of the raster data
        if (!is.null(values$points_sf)) {
          plot(values$points_sf, add = TRUE)
          coords <- sf::st_coordinates(values$points_sf)
          text(coords[,1], coords[,2], labels = points_data()[[input$idField]], pos = 4, cex = 0.8)
        }
        
        # If paths have been calculated, plot them on top of the points
        if (!is.null(values$lcps)) {
          plot(terra::vect(values$lcps), add = TRUE)
        }
      }
    })
  })
  
  output$downloadPaths <- downloadHandler(
    filename = function() {
      paste("paths_", input$lcpMethod, "_", input$costFunction, "_", Sys.Date(), ".gpkg", sep="")
    },
    content = function(file) {
      sf::st_write(values$lcps, file)
    }
  )
  
  observeEvent(input$createDensity, {
    req(values$lcps)

    # Calculate the path density
    density <- leastcostpath::create_lcp_density(uploaded_dem(), values$lcps)
    raster_data$data[["Path density map"]] <- density
    #raster_data$data[["Path density map"]] <- rasterise(density)
    updateSelectInput(session, "rasterData", choices = names(raster_data$data))

    observe(
      if(!is.null(density)){
        showNotification("Path density calculated successfully!")
      })
  })
  
  
  observeEvent(input$calculateCorridors, {
    
    # Check the selected least cost path method
    if (input$lcpMethod == "Sequential") {
      
      # Loop over the points sequentially
      for(i in 1:(nrow(values$points_sf) - 1)) {
        # Create a cost corridor for each pair of points and store it in the list
        cc[[i]] <- create_cost_corridor(slope_cs, values$points_sf[i,], values$points_sf[i + 1,])
      }
      
    } else if (input$lcpMethod == "From one to everywhere") {
      
      # Find the origin point
      origin <- values$points_sf[which(points_data()[[input$idField]] == input$idValue), ]
      
      # Loop over the other points
      for(i in 1:nrow(values$points_sf)) {
        if (i != which(points_data()[[input$idField]] == input$idValue)) {
          # Create a cost corridor from the origin to each other point and store it in the list
          cc[[i]] <- create_cost_corridor(slope_cs, origin, values$points_sf[i,])
        }
      }
      
    } else if (input$lcpMethod == "From everywhere to everywhere") {
      
      # Initialize a counter
      counter <- 1
      
      # Loop over all pairs of points
      for(i in 1:(nrow(values$points_sf) - 1)) {
        for(j in (i + 1):nrow(values$points_sf)) {
          # Create a cost corridor for each pair of points and store it in the list
          cc[[counter]] <- create_cost_corridor(slope_cs, values$points_sf[i,], values$points_sf[j,])
          counter <- counter + 1
        }
      }
      
    }
    
    # Store each cost corridor as a separate raster in your raster_data$data list and update the selectInput choices
    for(i in seq_along(cc)) {
      raster_data$data[[paste0("Cost Corridor ", i)]] <- cc[[i]]
    }
    updateSelectInput(session, "rasterData", choices = names(raster_data$data))

    observe(
      if(!is.null(cc)){
        showNotification("Cost corridors calculated successfully!")
      })
  })
  

  output$downloadRaster <- downloadHandler(
    filename = function() {
      paste(input$rasterData,"_", gsub(" ", "-", input$costFunction), "_", Sys.Date(), ".tif", sep="")
    },
    content = function(file) {
      # Get the currently selected raster data
      selected_data <- input$rasterData
        terra::writeRaster(
          raster_data$data[[selected_data]],
          file)
    }
  )
  output$downloadRaster2 <- downloadHandler(
    filename = function() {
      paste(input$rasterData,"_", gsub(" ", "-", input$costFunction), "_", Sys.Date(), ".tif", sep="")
    },
    content = function(file) {
      # Get the currently selected raster data
      selected_data <- input$rasterData
      terra::writeRaster(
        raster_data$data[[selected_data]],
        file)
    }
  )
  
  # Observe the raster_data$data list for changes
  observe({
    # Get the names of the rasters
    raster_names <- names(raster_data$data)
    
    # If there are any rasters
    if (length(raster_names) > 0) {
      # Print the names and classes of the rasters
      print("Names and classes of uploaded rasters:")
      for (raster_name in raster_names) {
        rast <- raster_data$data[[raster_name]]
        print(paste0(raster_name, " (", class(rast), ")"))
      }
    }
  })
  
  
  
  
  
} # close server

shinyApp(ui, server)
