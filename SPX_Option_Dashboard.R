# app.R
library(shiny)
library(quantmod) # For fetching financial data
library(plotly)   # For interactive plots
library(readxl)   # You already have this for local data, but we'll use quantmod
library(DT)       # For displaying data tables if needed
library(dplyr)    # For data manipulation (e.g., %>% and filter)

# --- Placeholder for subdiagCreate if subdiag package is not available --

# Define UI for application
ui <- fluidPage(
  title = "Option Pricing: Numerical Methods",
  tags$h1("Option Pricing: Numerical Methods"),
  tags$p("Author: Godfrey Nkolokosa"),
  tags$p(paste("Date:", Sys.Date())),
  
  navbarPage(
    title = "Option Pricing Dashboard",
    tabPanel("Executive Summary",
             fluidRow(
               column(12,
                      h2("Executive Summary"),
                      p("The global derivatives market is gigantic, estimated at $1 quadrillion on the high end (Maverick, 2022). Options are amongst the most useful and highly traded derivatives, mostly due to their application in portfolio hedging, speculation, and income generation among both institutional and retail investors. Thus, there has been growing demand for reliable and computationally effective techniques for valuing options. In this piece, I analyse and compare popular methods for valuing options using Yahoo Finance's S&P500 data. These methods include Monte Carlo simulation, Black-Scholes, and the Binomial Trees. Among other things, I find the Monte Carlo simulation and the Black-Scholes models to be robust and stable option pricing techniques."),
                      br(),
                      h3("Data Source:"),
                      p("For this dashboard, S&P 500 (SPX), VIX, and 13-week T-Bill yield (IRX) data are fetched dynamically using the 'quantmod' package, ensuring the latest 252 trading days are used for calculations."),
                      p("Please note: The original script used specific Excel files. This dashboard adapts to fetch live data.")
               )
             )
    ),
    # Removed FDM Tab
    tabPanel("Black-Scholes-Merton",
             sidebarLayout(
               sidebarPanel(
                 h3("Black-Scholes-Merton Parameters"),
                 selectInput("bsm_option_type", "Option Type:",
                             choices = c("Call" = FALSE, "Put" = TRUE)),
                 sliderInput("bsm_strike_price", "Strike Price (K):",
                             min = 4500, max = 6500, value = 5500, step = 50) # Updated range
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Interactive Plot", plotlyOutput("bsm_plot")),
                   tabPanel("Theory",
                            h2("Black-Scholes-Merton (BSM) Model"),
                            p("The Black-Scholes formula, also known as the Black-Scholes-Merton (BSM) model, is a mathematical formula used to calculate the theoretical price of European-style options. BSM was developed by Fischer Black, Myron Scholes, and Robert Merton in the 1970s and is largely recognised as one of the major breakthroughs in the field of finance."),
                            p("This section implements the BSM model to determine the value of European calls and puts using the SP500 index as the underlying asset and 1 year maturity.")
                   )
                 )
               )
             )
    ),
    tabPanel("Monte Carlo Simulation",
             sidebarLayout(
               sidebarPanel(
                 h3("Monte Carlo Parameters"),
                 selectInput("mc_option_type", "Option Type:",
                             choices = c("Call" = FALSE, "Put" = TRUE)),
                 sliderInput("mc_strike_price", "Strike Price (K):",
                             min = 4500, max = 6500, value = 5500, step = 50) # Updated range
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Interactive Plot", plotlyOutput("mc_plot")),
                   tabPanel("Theory",
                            h2("Monte Carlo Simulation"),
                            p("Monte Carlo simulations, also known as multiple probability simulations, is another method of pricing options especially when it is difficult or impossible to find analytical solutions. To implement MC, we estimated the future value of the SP500 index using the geometrical Brownian model based on the current index value (S[i]), time to expiration (tau), risk-free interest rate (r), and stock volatility (sigma[i]). For each index value (S[i]) at time t[i], we generated 1000 possible future stock price paths."),
                            p("For each simulated index value path, the payoff of the option was calculated based on its type (call or put) and strike price (K). If a call, the payoff was given as the maximum of zero and the difference between the simulated index value and the strike price. If it's a put option, the payoff is the maximum of zero and the difference between the strike price and the simulated index value.")
                   )
                 )
               )
             )
    ),
    tabPanel("Binomial Tree",
             sidebarLayout(
               sidebarPanel(
                 h3("Binomial Tree Parameters"),
                 selectInput("bt_option_type", "Option Type:",
                             choices = c("Call" = FALSE, "Put" = TRUE)),
                 selectInput("bt_option_style", "Option Style:",
                             choices = c("American" = "American", "European" = "European")),
                 sliderInput("bt_strike_price", "Strike Price (K):",
                             min = 4500, max = 6500, value = 5500, step = 50) # Updated range
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Interactive Plot", plotlyOutput("bt_plot")),
                   tabPanel("Theory",
                            h2("Binomial Tree Model"),
                            p("The Binomial Tree model is a discrete-time model used for pricing options. It simplifies the movement of the underlying asset's price over time into a series of upward or downward movements, forming a 'tree' of possible price paths."),
                            p("This model is particularly flexible for pricing American options because it allows for the possibility of early exercise at each node of the tree.")
                   )
                 )
               )
             )
    ),
    tabPanel("Benchmarks",
             tabsetPanel(
               tabPanel("Price Comparison",
                        sidebarLayout(
                          sidebarPanel(
                            h3("Price Comparison Parameters"),
                            selectInput("benchmark_option_type", "Option Type:",
                                        choices = c("Call" = FALSE, "Put" = TRUE)),
                            sliderInput("benchmark_strike_price", "Strike Price (K):",
                                        min = 4500, max = 6500, value = 5500, step = 50), # Updated range
                            selectInput("benchmark_option_style", "Binomial Tree Style:",
                                        choices = c("American" = "American", "European" = "European"),
                                        selected = "European") # Default for BT for easier comparison
                          ),
                          mainPanel(
                            h2("Option Price Comparison Over Time"),
                            plotlyOutput("benchmark_plot"),
                            h3("Model Insights"),
                            p("This chart compares the option prices generated by Black-Scholes, Monte Carlo, and Binomial Tree models over the last 252 trading days. Key observations:"),
                            tags$ul(
                              tags$li("Black-Scholes provides a theoretical, continuous-time solution for European options."),
                              tags$li("Monte Carlo simulation, while more flexible for complex options, can exhibit some variance due to its probabilistic nature."),
                              tags$li("Binomial Tree (European style) should generally converge towards the Black-Scholes price for European options with a sufficient number of steps. The American style Binomial Tree will always be greater than or equal to its European counterpart due to the early exercise feature.")
                            )
                          )
                        )
               ),
               tabPanel("Pricing Errors",
                        sidebarLayout(
                          sidebarPanel(
                            h3("Pricing Error Parameters"),
                            selectInput("error_option_type", "Option Type:",
                                        choices = c("Call" = FALSE, "Put" = TRUE)),
                            sliderInput("error_strike_price", "Strike Price (K):",
                                        min = 4500, max = 6500, value = 5500, step = 50), # Updated range
                            # For error comparison, Binomial Tree should be European to compare to Black-Scholes
                            selectInput("error_bt_style", "Binomial Tree Style:",
                                        choices = c("European" = "European"),
                                        selected = "European")
                          ),
                          mainPanel(
                            h2("Pricing Error Analysis (vs. Black-Scholes)"),
                            plotlyOutput("pricing_error_plot"),
                            h3("Understanding Pricing Errors"),
                            p("This plot illustrates the absolute difference between:"),
                            tags$ul(
                              tags$li("Monte Carlo option prices and Black-Scholes option prices."),
                              tags$li("Binomial Tree (European) option prices and Black-Scholes option prices."),
                              tags$li("The Black-Scholes model is used as the benchmark as it provides an analytical solution for European options. Ideally, other numerical methods for European options (like European Binomial Tree and Monte Carlo) should converge to Black-Scholes prices with sufficient iterations/steps.")
                            )
                          )
                        )
               )
             )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # Reactive data fetching and processing
  option_data <- reactive({
    # Fetch S&P500 (SPX) data
    getSymbols("^GSPC", src = "yahoo",
               from = Sys.Date() - 747, # Ensure enough historical data for 252 trading days
               to = Sys.Date(),
               auto.assign = TRUE) # Ensures ^GSPC is assigned to GSPC object
    
    # Fetch VIX data
    getSymbols("^VIX", src = "yahoo",
               from = Sys.Date() - 747, # Ensure enough historical data
               to = Sys.Date(),
               auto.assign = TRUE) # Ensures ^VIX is assigned to VIX object
    
    # Fetch 13-week T-bill yield (risk-free rate)
    getSymbols("^IRX", src = "yahoo",
               from = Sys.Date() - 747, # Ensure enough historical data
               to = Sys.Date(),
               auto.assign = TRUE) # Ensures ^IRX is assigned to IRX object
    
    # Extract adjusted close prices and clean data
    S_raw <- tryCatch(Cl(GSPC), error = function(e) {
      message("Error fetching GSPC: ", e$message)
      return(NULL)
    })
    VIX_raw <- tryCatch(Cl(VIX), error = function(e) {
      message("Error fetching VIX: ", e$message)
      return(NULL)
    })
    IRX_raw <- tryCatch(Cl(IRX), error = function(e) {
      message("Error fetching IRX: ", e$message)
      return(NULL)
    })
    
    if (is.null(S_raw) || is.null(VIX_raw) || is.null(IRX_raw)) {
      stop("Failed to fetch one or more financial symbols. Please check internet connection or symbol availability.")
    }
    
    # Combine data by date and remove NA values
    merged_data <- merge(S_raw, VIX_raw, IRX_raw, all = FALSE) %>% na.omit()
    
    # Get the last 252 trading days
    if(nrow(merged_data) < 252) {
      n_obs <- nrow(merged_data)
      showNotification(paste0("Only ", n_obs, " trading days available for calculations. Data might be insufficient."), type = "warning")
    } else {
      n_obs <- 252
    }
    
    # Extract latest 252 observations
    S_latest <- as.vector(tail(merged_data[, 1], n_obs))
    VIX_latest <- as.vector(tail(merged_data[, 2], n_obs))
    IRX_latest <- as.vector(tail(merged_data[, 3], n_obs))
    
    # Convert VIX to daily volatility (annualized percentage to daily decimal)
    sigma_daily <- (VIX_latest / 100) / sqrt(252)
    
    # Risk-free rate (use the latest IRX value, convert from annualized percentage to decimal)
    r_annualized <- as.numeric(tail(IRX_latest, 1)) / 100
    
    list(S = S_latest, sigma = sigma_daily, r = r_annualized, n = n_obs, time = 1:n_obs)
  })
  
  
  ### Black Scholes Merton
  cblackscholes1 <- function(K, isput, sigma, S, r) {
    T_mat <- 1
    n <- length(S)
    t <- seq(0, 1, length.out = n)
    call_prices <- numeric(n)
    for (i in 1:n) {
      tau <- T_mat - t[i]
      if (tau <= 0) {
        call_prices[i] <- ifelse(isput, max(0, K - S[i]), max(0, S[i] - K))
        next
      }
      
      if (sigma[i] <= 0 || is.na(sigma[i])) {
        if (!isput) {
          call_prices[i] <- max(0, S[i] - K) * exp(-r * tau)
        } else {
          call_prices[i] <- max(0, K - S[i]) * exp(-r * tau)
        }
        next
      }
      
      d1 <- 1 / (sigma[i] * sqrt(tau)) * (log(S[i] / K) + (r + sigma[i]^2 / 2) * tau)
      d2 <- d1 - sigma[i] * sqrt(tau)
      
      if (!is.finite(d1) || !is.finite(d2)) {
        call_prices[i] <- 0
        next
      }
      
      if (!isput) { # Call
        call_prices[i] <- pnorm(d1) * S[i] - pnorm(d2) * K * exp(-r * tau)
      } else { # Put
        call_prices[i] <- pnorm(-d2) * K * exp(-r * tau) - pnorm(-d1) * S[i]
      }
    }
    return(call_prices)
  }
  
  ### Monte Carlo
  montecarlo1 <- function(K, isput, sigma, S, r) {
    T_mat <- 1
    N_sim <- 1000 # Number of simulations
    n <- length(S)
    t <- seq(0, 1, length.out = n)
    option_prices <- numeric(n)
    
    for (i in 1:n) {
      tau <- T_mat - t[i]
      if (tau <= 0) {
        option_prices[i] <- ifelse(isput, max(0, K - S[i]), max(0, S[i] - K))
        next
      }
      
      if (sigma[i] <= 0 || is.na(sigma[i])) {
        st <- S[i] * exp(r * tau)
      } else {
        st <- S[i] * exp((r - 0.5 * sigma[i]^2) * tau + sigma[i] * rnorm(N_sim) * sqrt(tau))
      }
      
      if (!isput) { # Call
        payoff <- pmax(st - K, 0)
      } else { # Put
        payoff <- pmax(K - st, 0)
      }
      
      option_prices[i] <- exp(-r * tau) * mean(payoff)
    }
    return(option_prices)
  }
  
  ### Binomial Tree - American Options
  BT_am <- function(K, isput, sigma, S, r) {
    T_mat <- 1
    N_steps <- 50 # Number of tree periods
    n_days <- length(S)
    t_days <- seq(0, 1, length.out = n_days)
    option_prices <- numeric(n_days)
    
    for (o in 1:n_days) {
      tau <- T_mat - t_days[o]
      if (tau <= 0) {
        option_prices[o] <- ifelse(isput, max(0, K - S[o]), max(0, S[o] - K))
        next
      }
      
      dt <- tau / N_steps
      if (sigma[o] <= 0 || is.na(sigma[o])) {
        u <- exp(r * dt)
        d <- 1 / u
      } else {
        u <- exp(sigma[o] * sqrt(dt))
        d <- 1 / u
      }
      
      q <- (exp(r * dt) - d) / (u - d)
      
      if (q < 0 || q > 1 || is.na(q) || is.infinite(q)) {
        option_prices[o] <- 0
        next
      }
      
      M <- matrix(0, nrow = N_steps + 1, ncol = N_steps + 1)
      
      terminal_stock_prices_nodes <- S[o] * u^((N_steps):0) * d^(0:(N_steps))
      
      if (!isput) { # Call
        M[, N_steps + 1] <- pmax(terminal_stock_prices_nodes - K, 0)
      } else { # Put
        M[, N_steps + 1] <- pmax(K - terminal_stock_prices_nodes, 0)
      }
      
      for (j in N_steps:1) {
        for (i in 1:j) {
          expected_value <- exp(-r * dt) * (q * M[i, j + 1] + (1 - q) * M[i + 1, j + 1])
          
          current_stock_price_at_node <- S[o] * u^(j - i) * d^(i - 1)
          intrinsic_value <- ifelse(!isput,
                                    max(0, current_stock_price_at_node - K),
                                    max(0, K - current_stock_price_at_node))
          
          M[i, j] <- pmax(expected_value, intrinsic_value)
        }
      }
      option_prices[o] <- M[1, 1]
    }
    return(option_prices)
  }
  
  ### Binomial Tree - European Options
  BT_eu1 <- function(K, isput, sigma, S, r) {
    T_mat <- 1
    N_steps <- 50
    n_days <- length(S)
    t_days <- seq(0, 1, length.out = n_days)
    option_prices <- numeric(n_days)
    
    for (o in 1:n_days) {
      tau <- T_mat - t_days[o]
      if (tau <= 0) {
        option_prices[o] <- ifelse(isput, max(0, K - S[o]), max(0, S[o] - K))
        next
      }
      
      dt <- tau / N_steps
      if (sigma[o] <= 0 || is.na(sigma[o])) {
        u <- exp(r * dt)
        d <- 1 / u
      } else {
        u <- exp(sigma[o] * sqrt(dt))
        d <- 1 / u
      }
      q <- (exp(r * dt) - d) / (u - d)
      
      if (q < 0 || q > 1 || is.na(q) || is.infinite(q)) {
        option_prices[o] <- 0
        next
      }
      
      M <- matrix(0, nrow = N_steps + 1, ncol = N_steps + 1)
      
      terminal_stock_prices_nodes <- S[o] * u^((N_steps):0) * d^(0:(N_steps))
      
      if (!isput) { # Call
        M[, N_steps + 1] <- pmax(terminal_stock_prices_nodes - K, 0)
      } else { # Put
        M[, N_steps + 1] <- pmax(K - terminal_stock_prices_nodes, 0)
      }
      
      for (j in N_steps:1) {
        for (i in 1:j) {
          M[i, j] <- exp(-r * dt) * (q * M[i, j + 1] + (1 - q) * M[i + 1, j + 1])
        }
      }
      option_prices[o] <- M[1, 1]
    }
    return(option_prices)
  }
  
  
  # Reactive expressions for pricing each model
  bsm_prices <- reactive({
    req(option_data())
    data <- option_data()
    K_val <- input$bsm_strike_price
    is_put <- as.logical(input$bsm_option_type)
    cblackscholes1(K_val, is_put, data$sigma, data$S, data$r)
  })
  
  mc_prices <- reactive({
    req(option_data())
    data <- option_data()
    K_val <- input$mc_strike_price
    is_put <- as.logical(input$mc_option_type)
    montecarlo1(K_val, is_put, data$sigma, data$S, data$r)
  })
  
  bt_prices <- reactive({
    req(option_data())
    data <- option_data()
    K_val <- input$bt_strike_price
    is_put <- as.logical(input$bt_option_type)
    
    if (input$bt_option_style == "American") {
      BT_am(K_val, is_put, data$sigma, data$S, data$r)
    } else {
      BT_eu1(K_val, is_put, data$sigma, data$S, data$r)
    }
  })
  
  # Reactive for Benchmark Plot
  benchmark_data <- reactive({
    req(option_data())
    data <- option_data()
    K_val <- input$benchmark_strike_price
    is_put <- as.logical(input$benchmark_option_type)
    bt_style <- input$benchmark_option_style
    
    bsm_p <- cblackscholes1(K_val, is_put, data$sigma, data$S, data$r)
    mc_p <- montecarlo1(K_val, is_put, data$sigma, data$S, data$r)
    bt_p <- if (bt_style == "American") {
      BT_am(K_val, is_put, data$sigma, data$S, data$r)
    } else {
      BT_eu1(K_val, is_put, data$sigma, data$S, data$r)
    }
    
    list(
      time = data$time,
      bsm = bsm_p,
      mc = mc_p,
      bt = bt_p,
      sigma = data$sigma # Include sigma for VIX line
    )
  })
  
  # Reactive for Pricing Error Plot
  pricing_error_data <- reactive({
    req(option_data())
    data <- option_data()
    K_val <- input$error_strike_price
    is_put <- as.logical(input$error_option_type)
    # Ensure BT is European for error comparison against BSM
    bt_style <- "European"
    
    bsm_p <- cblackscholes1(K_val, is_put, data$sigma, data$S, data$r)
    mc_p <- montecarlo1(K_val, is_put, data$sigma, data$S, data$r)
    bt_p <- BT_eu1(K_val, is_put, data$sigma, data$S, data$r) # Always European for this tab
    
    # Calculate absolute differences
    mc_error <- abs(mc_p - bsm_p)
    bt_error <- abs(bt_p - bsm_p)
    
    list(
      time = data$time,
      mc_error = mc_error,
      bt_error = bt_error,
      sigma = data$sigma # Include sigma for VIX line
    )
  })
  
  
  # Output plots for BSM
  output$bsm_plot <- renderPlotly({
    req(bsm_prices())
    data <- option_data()
    prices <- bsm_prices()
    option_label <- ifelse(as.logical(input$bsm_option_type), "Put Prices (USD)", "Call Prices (USD)")
    plot_title <- paste("Black-Scholes:", ifelse(as.logical(input$bsm_option_type), "Puts", "Calls"), " (K =", input$bsm_strike_price, ")")
    
    plot_ly() %>%
      add_trace(x = data$time, y = prices, type = 'scatter', mode = 'lines',
                name = option_label, line = list(color = 'blue')) %>%
      add_trace(x = data$time, y = data$sigma * 100, type = 'scatter', mode = 'lines',
                name = "VIX (%)", yaxis = "y2", line = list(color = 'brown', dash = 'dash')) %>%
      layout(
        title = plot_title,
        xaxis = list(title = "Time Elapsed, Days"),
        yaxis = list(title = option_label),
        yaxis2 = list(title = "VIX (%)", overlaying = "y", side = "right", range = c(0, max(data$sigma * 100, na.rm = TRUE) * 1.2)),
        hovermode = "x unified",
        legend = list(orientation = "h", x = 0, y = -0.2) # Legend at bottom
      )
  })
  
  # Output plots for Monte Carlo
  output$mc_plot <- renderPlotly({
    req(mc_prices())
    data <- option_data()
    prices <- mc_prices()
    option_label <- ifelse(as.logical(input$mc_option_type), "Put Prices (USD)", "Call Prices (USD)")
    plot_title <- paste("Monte Carlo:", ifelse(as.logical(input$mc_option_type), "Puts", "Calls"), " (K =", input$mc_strike_price, ")")
    
    plot_ly() %>%
      add_trace(x = data$time, y = prices, type = 'scatter', mode = 'lines',
                name = option_label, line = list(color = 'blue')) %>%
      add_trace(x = data$time, y = data$sigma * 100, type = 'scatter', mode = 'lines',
                name = "VIX (%)", yaxis = "y2", line = list(color = 'brown', dash = 'dash')) %>%
      layout(
        title = plot_title,
        xaxis = list(title = "Time Elapsed, Days"),
        yaxis = list(title = option_label),
        yaxis2 = list(title = "VIX (%)", overlaying = "y", side = "right", range = c(0, max(data$sigma * 100, na.rm = TRUE) * 1.2)),
        hovermode = "x unified",
        legend = list(orientation = "h", x = 0, y = -0.2) # Legend at bottom
      )
  })
  
  # Output plots for Binomial Tree
  output$bt_plot <- renderPlotly({
    req(bt_prices())
    data <- option_data()
    prices <- bt_prices()
    option_label <- ifelse(as.logical(input$bt_option_type), "Put Prices (USD)", "Call Prices (USD)")
    plot_title <- paste("Binomial Tree:", input$bt_option_style, ifelse(as.logical(input$bt_option_type), "Puts", "Calls"), " (K =", input$bt_strike_price, ")")
    
    plot_ly() %>%
      add_trace(x = data$time, y = prices, type = 'scatter', mode = 'lines',
                name = option_label, line = list(color = 'blue')) %>%
      add_trace(x = data$time, y = data$sigma * 100, type = 'scatter', mode = 'lines',
                name = "VIX (%)", yaxis = "y2", line = list(color = 'brown', dash = 'dash')) %>%
      layout(
        title = plot_title,
        xaxis = list(title = "Time Elapsed, Days"),
        yaxis = list(title = option_label),
        yaxis2 = list(title = "VIX (%)", overlaying = "y", side = "right", range = c(0, max(data$sigma * 100, na.rm = TRUE) * 1.2)),
        hovermode = "x unified",
        legend = list(orientation = "h", x = 0, y = -0.2) # Legend at bottom
      )
  })
  
  # Output plot for Benchmarks (Price Comparison)
  output$benchmark_plot <- renderPlotly({
    req(benchmark_data())
    b_data <- benchmark_data()
    option_label <- ifelse(as.logical(input$benchmark_option_type), "Put Prices (USD)", "Call Prices (USD)")
    plot_title <- paste("Model Comparison:", ifelse(as.logical(input$benchmark_option_type), "Puts", "Calls"), " (K =", input$benchmark_strike_price, ")")
    
    plot_ly() %>%
      add_trace(x = b_data$time, y = b_data$bsm, type = 'scatter', mode = 'lines',
                name = "Black-Scholes", line = list(color = 'blue')) %>%
      add_trace(x = b_data$time, y = b_data$mc, type = 'scatter', mode = 'lines',
                name = "Monte Carlo", line = list(color = 'red')) %>%
      add_trace(x = b_data$time, y = b_data$bt, type = 'scatter', mode = 'lines',
                name = paste0("Binomial Tree (", input$benchmark_option_style, ")"), line = list(color = 'green')) %>%
      add_trace(x = b_data$time, y = b_data$sigma * 100, type = 'scatter', mode = 'lines',
                name = "VIX (%)", yaxis = "y2", line = list(color = 'brown', dash = 'dash')) %>%
      layout(
        title = plot_title,
        xaxis = list(title = "Time Elapsed, Days"),
        yaxis = list(title = option_label),
        yaxis2 = list(title = "VIX (%)", overlaying = "y", side = "right", range = c(0, max(b_data$sigma * 100, na.rm = TRUE) * 1.2)),
        hovermode = "x unified",
        legend = list(orientation = "h", x = 0, y = -0.2) # Legend at bottom
      )
  })
  
  # Output plot for Pricing Errors
  output$pricing_error_plot <- renderPlotly({
    req(pricing_error_data())
    error_data <- pricing_error_data()
    error_label <- "Absolute Pricing Error (USD)"
    plot_title <- paste("Pricing Error vs. Black-Scholes:", ifelse(as.logical(input$error_option_type), "Puts", "Calls"), " (K =", input$error_strike_price, ")")
    
    plot_ly() %>%
      add_trace(x = error_data$time, y = error_data$mc_error, type = 'scatter', mode = 'lines',
                name = "Monte Carlo Error", line = list(color = 'red')) %>%
      add_trace(x = error_data$time, y = error_data$bt_error, type = 'scatter', mode = 'lines',
                name = "Binomial Tree (European) Error", line = list(color = 'green')) %>%
      add_trace(x = error_data$time, y = error_data$sigma * 100, type = 'scatter', mode = 'lines',
                name = "VIX (%)", yaxis = "y2", line = list(color = 'brown', dash = 'dash')) %>%
      layout(
        title = plot_title,
        xaxis = list(title = "Time Elapsed, Days"),
        yaxis = list(title = error_label),
        yaxis2 = list(title = "VIX (%)", overlaying = "y", side = "right", range = c(0, max(error_data$sigma * 100, na.rm = TRUE) * 1.2)),
        hovermode = "x unified",
        legend = list(orientation = "h", x = 0, y = -0.2) # Legend at bottom
      )
  })
  
}

# Run the application
shinyApp(ui = ui, server = server)