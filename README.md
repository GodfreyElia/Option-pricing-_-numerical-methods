# 📈 Option Pricing Dashboard

This **Shiny application** provides an interactive dashboard for exploring and comparing different numerical methods used to price **European and American options**. It dynamically fetches **live market data** for the **S&P 500 (^GSPC)**, **VIX Index (^VIX)**, and **13-week U.S. Treasury Bill yield (^IRX)** from **Yahoo Finance**, allowing for up-to-date analysis. See the live dashboard [here](https://godfreyelia.shinyapps.io/SPX_Option_Pricing/)

---

## 🚀 Features

- **📡 Dynamic Data Fetching**  
  Automatically retrieves the latest S&P 500 prices, VIX volatility, and risk-free rates from Yahoo Finance.

- **📊 Black-Scholes-Merton (BSM) Model**  
  Calculate theoretical prices for European call and put options.

- **🎲 Monte Carlo Simulation**  
  Price options by simulating numerous future price paths of the underlying asset.

- **🌲 Binomial Tree Model**  
  Supports both European and American call and put options, demonstrating the flexibility for early exercise scenarios.

- **📈 Interactive Plots**  
  Visualize option prices and VIX movements over the last 252 trading days (~1 year).

- **🧪 Benchmarking Tab**
  - **Price Comparison**: Compare BSM, Monte Carlo, and Binomial Tree prices in a single chart.
  - **Pricing Errors**: Analyze absolute differences between Monte Carlo & Binomial (European) vs. BSM benchmark.

- **⚙️ Customizable Parameters**  
  Adjust strike prices, option types (call/put), and styles (American/European for Binomial Tree) with sliders and dropdowns.

---

## 🧩 How to Use

### 1. Install R and RStudio  
If not already installed, download from:
- [R](https://cran.r-project.org/)
- [RStudio](https://posit.co/download/rstudio-desktop/)

### 2. Install Required R Packages  
Open RStudio and run:

```R
install.packages(c("shiny", "quantmod", "plotly", "readxl", "DT", "dplyr"))
```
## Live Dashboard

See the live dashboard [here](https://godfreyelia.shinyapps.io/SPX_Option_Pricing/)
