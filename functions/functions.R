f_load_libraries <- function(){
  # ---------------------------------------------------------------------------
  # 1. Data Manipulation & Workflow 
  # ---------------------------------------------------------------------------
  library(tidyverse)   # Data manipulation, visualization, and more
  library(here)        # Simplifies file path handling
  library(reshape2)    # Data reshaping (deprecated, use tidyr instead)
  library(plan)        # Workflow management
  
  # ---------------------------------------------------------------------------
  # 2. Statistical Modeling & Bayesian Analysis 
  # ---------------------------------------------------------------------------
  library(brms)        # Bayesian regression models using Stan
  library(cmdstanr)    # Interface to CmdStan
  library(rstan)       # Interface to Stan for Bayesian modeling
  library(glmmTMB)     # Generalized linear mixed models
  library(ggmcmc)      # MCMC visualization
  library(tidybayes)   # Tidy tools for Bayesian models
  library(fdrtool)     # False Discovery Rate corrections
  library(tmvtnorm)    # Truncated multivariate normal distribution
  library(truncnorm)   # Truncated normal distribution
  library(priorsense)  # Sensitivity analysis for Bayesian priors
  library(deSolve)     # Solving differential equations (ODEs, PDEs, DDEs)
  library(easystats)   # Tools for statistical modeling and reporting
  
  # ---------------------------------------------------------------------------
  # 3. Data Visualization & Graphics
  # ---------------------------------------------------------------------------
  library(ggplot2)     # Core data visualization (part of tidyverse)
  library(ggthemes)    # Additional themes for ggplot2
  library(viridis)     # Color palettes optimized for accessibility
  library(wesanderson) # Wes Anderson-inspired color palettes
  library(RColorBrewer) # Color palettes for categorical & sequential data
  library(nord)        # Minimalist "Nord" color palettes
  library(ggdist)      # Visualization of distributions and uncertainty
  library(patchwork)   # Arrange multiple ggplots
  library(cowplot)     # Streamlined plot composition
  library(ggbreak)     # Break axes in ggplot2
  library(ggtext)      # Enhanced text rendering in ggplot2
  library(plotly)      # Interactive plots
  library(ggpubr)      # Publication-ready figures
  library(gridExtra)   # Arrange multiple plots
  library(grid)        # Low-level graphics functions
  library(extrafont)   # Additional fonts for ggplot
  
  # ---------------------------------------------------------------------------
  # 4. Model Summary & Reporting
  # ---------------------------------------------------------------------------
  library(modelsummary) # Summarize regression models
  library(car)         # Regression diagnostics and statistical tests
  library(Hmisc)       # Data summarization & missing data handling
  library(lsmeans)     # Least-squares means comparisons (replaced by emmeans)
  library(nlstools)    # Nonlinear regression diagnostics
  library(latex2exp)   # Use LaTeX expressions in plots
  library(knitr)       # Dynamic reports (Markdown/Quarto)
  library(kableExtra)  # Extended table formatting for knitr::kable()
  library(DT)         # Interactive tables (DataTables)
  
  # ---------------------------------------------------------------------------
  # 5. File Handling & External Libraries
  # ---------------------------------------------------------------------------
  library(readxl)      # Read Excel files
  library(devtools)    # Install and manage R packages
  
  # ---------------------------------------------------------------------------
  # 6. Parallel Computing & Performance
  # ---------------------------------------------------------------------------
  library(parallel)    # Multi-core processing
}


f_load_colors <- function(){
  Nord_frost <<- nord(palette = "frost")
  Nord_polar <<- nord(palette = "polarnight")
  Nord_aurora <<- nord(palette = "aurora")
  sizetitle <<- 12
  
  col_blue <<- "#5E81AC"
  col_red <<- "#f42404"
  pal_col <<- c(col_blue, col_red)
  
  pal_blue <<- c("#5E81AC", "#7F9DC4", "#A0C1D9", "#DCE9F2")
  pal_red <<- c("#f42404", "#F65E4B", "#F6876D", "#FBD3D0")  
}

f_load_libraries_colors <- function(){
  f_load_libraries()
  f_load_colors()
}