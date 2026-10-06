# ===================================================================
# INTERAKTIVNÍ SHINY APLIKACE: POPULAČNÍ GENETIKA, HWE A EVOLUČNÍ SÍLY
# (Punnettův čtverec, HWE 2 alely, HWE 3 alely / ABO, Selekce, Mutace, Migrace)
# Version 6 - Přidán dynamický Punnettův čtverec s plochami genotypů na úvod
# ===================================================================
# install.packages("shinylive")
# shinylive::export("cesta/k/aplikaci", "docs")
if (!requireNamespace("shiny", quietly = TRUE)) install.packages("shiny")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(shiny)
library(ggplot2)

# Bezpečné nastavení tématu vzhledu
app_theme <- if (requireNamespace("shinythemes", quietly = TRUE)) {
  shinythemes::shinytheme("flatly")
} else {
  NULL
}

# -------------------------------------------------------------------
# 1. MATEMATICKÉ FUNKCE PRO HWE & EVOLUČNÍ SÍLY
# -------------------------------------------------------------------

# A) Výpočet HWE pro 2 alely (n_AA, n_Aa, n_aa)
vypocti_hwe_2alely <- function(n_AA, n_Aa, n_aa) {
  n_AA <- as.numeric(n_AA)
  n_Aa <- as.numeric(n_Aa)
  n_aa <- as.numeric(n_aa)
  
  obs <- c(AA = n_AA, Aa = n_Aa, aa = n_aa)
  N <- sum(obs)
  if (is.na(N) || N == 0) return(NULL)
  
  p <- as.numeric((2 * n_AA + n_Aa) / (2 * N))
  q <- as.numeric(1 - p)
  
  exp_counts <- c(
    AA = as.numeric(N * (p^2)),
    Aa = as.numeric(N * (2 * p * q)),
    aa = as.numeric(N * (q^2))
  )
  
  chi2_val <- as.numeric(sum((obs - exp_counts)^2 / ifelse(exp_counts == 0, 1e-6, exp_counts)))
  df <- 1
  p_val <- as.numeric(pchisq(chi2_val, df = df, lower.tail = FALSE))
  
  return(list(
    N = N,
    obs = obs,
    obs_freq = obs / N,
    p = p,
    q = q,
    exp_counts = exp_counts,
    exp_freq = exp_counts / N,
    chi2 = chi2_val,
    p_val = p_val
  ))
}

# B) Výpočet HWE pro 3 alely (např. ABO: p = A, q = B, r = O)
vypocti_hwe_3alely <- function(p, q, N = 1000, obs_counts = NULL) {
  p <- as.numeric(p)
  q <- as.numeric(q)
  r <- as.numeric(1 - p - q)
  if (r < 0) r <- 0
  
  gen_freq <- c(
    AA = p^2,
    BB = q^2,
    OO = r^2,
    AB = 2 * p * q,
    AO = 2 * p * r,
    BO = 2 * q * r
  )
  
  phen_freq <- c(
    "Skupina A (AA + AO)" = p^2 + 2 * p * r,
    "Skupina B (BB + BO)" = q^2 + 2 * q * r,
    "Skupina AB (AB)"    = 2 * p * q,
    "Skupina O (OO)"      = r^2
  )
  
  exp_gen_counts <- gen_freq * N
  exp_phen_counts <- phen_freq * N
  
  chi2_gen <- NULL
  p_val_gen <- NULL
  if (!is.null(obs_counts) && length(obs_counts) == 6) {
    chi2_gen <- sum((obs_counts - exp_gen_counts)^2 / ifelse(exp_gen_counts == 0, 1e-6, exp_gen_counts))
    df_gen <- 6 - 2 - 1
    p_val_gen <- pchisq(chi2_gen, df = df_gen, lower.tail = FALSE)
  }
  
  return(list(
    p = p, q = q, r = r, N = N,
    gen_freq = gen_freq,
    phen_freq = phen_freq,
    exp_gen_counts = exp_gen_counts,
    exp_phen_counts = exp_phen_counts,
    chi2_gen = chi2_gen,
    p_val_gen = p_val_gen
  ))
}

# C) Výpočet jednoho evolučního kroku (Selekce, Mutace, Migrace pro 2 alely)
vypocti_dalsi_generaci <- function(q, W_AA, W_Aa, W_aa, u, v, m, q_m) {
  p <- 1 - q
  
  # Selekce
  W_bar <- (p^2 * W_AA) + (2 * p * q * W_Aa) + (q^2 * W_aa)
  if (W_bar > 0) {
    q_sel <- (p * q * W_Aa + q^2 * W_aa) / W_bar
  } else {
    q_sel <- q
  }
  
  # Mutace
  q_mut <- q_sel * (1 - v) + (1 - q_sel) * u
  
  # Migrace
  q_next <- (1 - m) * q_mut + m * q_m
  q_next <- max(0, min(1, q_next))
  
  return(list(
    q_next = q_next,
    delta_q = q_next - q,
    W_bar = W_bar
  ))
}

# D) Simulace evoluce v čase
simuluj_populaci <- function(q0, W_AA, W_Aa, W_aa, u, v, m, q_m, generace) {
  q_vec <- numeric(generace + 1)
  W_vec <- numeric(generace + 1)
  q_vec[1] <- q0
  
  for (t in 1:generace) {
    res <- vypocti_dalsi_generaci(q_vec[t], W_AA, W_Aa, W_aa, u, v, m, q_m)
    q_vec[t + 1] <- res$q_next
    W_vec[t] <- res$W_bar
  }
  p_last <- 1 - q_vec[generace + 1]
  W_vec[generace + 1] <- (p_last^2 * W_AA) + (2 * p_last * q_vec[generace + 1] * W_Aa) + (q_vec[generace + 1]^2 * W_aa)
  
  data.frame(
    Generace = 0:generace,
    q = q_vec,
    p = 1 - q_vec,
    W_bar = W_vec
  )
}

# -------------------------------------------------------------------
# 2. UŽIVATELSKÉ ROZHRANÍ (UI - NAVBAR)
# -------------------------------------------------------------------
ui <- navbarPage(
  theme = app_theme,
  title = "Genetika populací: HWE & Evoluční síly",
  footer = tags$div(
    class = "text-muted",
    style = "text-align: center; padding: 12px; margin-top: 20px;
             border-top: 1px solid #ddd; font-size: 0.9em;",
    HTML(paste0("Autor: <b>Tomáš Urban</b>, MENDELU &copy; ",
                format(Sys.Date(), "%Y")))
  ),  
  # =================================================================
  # ZÁLOŽKA 1: DYNAMICKÝ PUNNETTŮV ČTVEREC (PLOCHY GENOTYPŮ)
  # =================================================================
  tabPanel(
    "1. Punnettův čtverec (Plochy genotypů)",
    sidebarLayout(
      sidebarPanel(
        width = 4,
        h3("Geometrická vizualizace HWE"),
        helpText("Geometrický model HWE: Plocha každého čtverce/obdélníku v Punnettově čtverci odpovídá frekvenci daného genotypu potomků."),
        hr(),
        sliderInput("punnett_p", "Frekvence alely A (p):", min = 0.01, max = 0.99, value = 0.60, step = 0.01),
        uiOutput("punnett_q_text"),
        hr(),
        numericInput("punnett_N", "Velikost populace pro přepočet (N):", value = 1000, min = 100, max = 100000, step = 100),
        hr(),
        wellPanel(
          h4("Rozpad genotypových frekvencí:"),
          htmlOutput("punnett_summary_html")
        )
      ),
      
      mainPanel(
        width = 8,
        fluidRow(
          column(12, plotOutput("plot_punnett_square", height = "420px"))
        ),
        br(),
        fluidRow(
          column(12, plotOutput("plot_punnett_bar", height = "220px"))
        )
      )
    )
  ),
  
  # =================================================================
  # ZÁLOŽKA 2: FREKVENCE ALEL, GENOTYPŮ A HWE TESTOVÁNÍ (2 ALELY)
  # =================================================================
  tabPanel(
    "2. HWE model a testování (2 alely)",
    sidebarLayout(
      sidebarPanel(
        width = 4,
        h3("Parametry HWE a populace"),
        radioButtons("hwe_mode", "Režim zadávání dat:",
                     choices = c(
                       "Teoretický model (přesné HWE frekvence)" = "theory",
                       "Simulace náhodného vzorku (N jedinců)" = "sample",
                       "Ruční zadání pozorovaných genotypů" = "manual"
                     )),
        hr(),
        
        conditionalPanel(
          condition = "input.hwe_mode == 'theory' || input.hwe_mode == 'sample'",
          sliderInput("hwe_p", "Frekvence alely A (p):", min = 0.01, max = 0.99, value = 0.6, step = 0.01),
          helpText("Frekvence alely a je q = 1 - p.")
        ),
        
        conditionalPanel(
          condition = "input.hwe_mode == 'sample'",
          numericInput("hwe_N", "Velikost vzorku populace (N):", value = 200, min = 10, max = 10000, step = 50),
          actionButton("btn_resample", "Vygenerovat nový náhodný vzorek", class = "btn-primary"),
          br(), br()
        ),
        
        conditionalPanel(
          condition = "input.hwe_mode == 'manual'",
          h4("Pozorované počty genotypů:"),
          numericInput("manual_AA", "Počet AA (D):", value = 120, min = 0, step = 1),
          numericInput("manual_Aa", "Počet Aa (H):", value = 60, min = 0, step = 1),
          numericInput("manual_aa", "Počet aa (R):", value = 20, min = 0, step = 1)
        ),
        
        hr(),
        wellPanel(
          h4("Shrnutí modelu:"),
          htmlOutput("hwe_summary_text")
        )
      ),
      
      mainPanel(
        width = 8,
        fluidRow(
          column(6, plotOutput("plot_hwe_curves", height = "360px")),
          column(6, plotOutput("plot_hwe_counts", height = "360px"))
        ),
        br(),
        fluidRow(
          column(12,
                 wellPanel(
                   h4("Statistické vyhodnocení Chi-kvadrát (χ²) testu HWE"),
                   tableOutput("tabulka_hwe_stats"),
                   htmlOutput("hwe_test_result")
                 )
          )
        )
      )
    )
  ),
  
  # =================================================================
  # ZÁLOŽKA 3: VÍCEALELOVÝ SYSTÉM (3 ALELY - ABO)
  # =================================================================
  tabPanel(
    "3. Vícealelový systém (ABO)",
    sidebarLayout(
      sidebarPanel(
        width = 4,
        h3("Parametry 3-alelového systému"),
        helpText("Model HWE pro 3 alely: (p + q + r)² = 1"),
        hr(),
        sliderInput("abo_p", "Frekvence alely I^A (p):", min = 0.01, max = 0.90, value = 0.30, step = 0.01),
        uiOutput("ui_abo_q"),
        wellPanel(
          htmlOutput("abo_r_text")
        ),
        numericInput("abo_N", "Velikost populace (N):", value = 1000, min = 100, max = 50000, step = 100)
      ),
      
      mainPanel(
        width = 8,
        fluidRow(
          column(6, plotOutput("plot_3allele_genotypes", height = "360px")),
          column(6, plotOutput("plot_3allele_phenotypes", height = "360px"))
        ),
        br(),
        fluidRow(
          column(12,
                 wellPanel(
                   h4("Přehled genotypových a fenotypových frekvencí (ABO)"),
                   tableOutput("tabulka_3allele_summary")
                 )
          )
        )
      )
    )
  ),
  
  # =================================================================
  # ZÁLOŽKA 4: EVOLUČNÍ SÍLY (SELEKCE, MUTACE, MIGRACE)
  # =================================================================
  tabPanel(
    "4. Dynamika evolučních sil",
    sidebarLayout(
      sidebarPanel(
        width = 4,
        selectInput("preset", "Vyberte výukový scénář (Preset):",
                    choices = c(
                      "Vlastní nastavení" = "custom",
                      "Selekce: Letální recesivní alela (W_aa = 0)" = "sel_recessive",
                      "Selekce: Naddominance / Balancovaný polymorfismus" = "overdominance",
                      "Mutace: Mutační rovnováha (u vs. v)" = "mutation_eq",
                      "Mutace + Selekce: Mutačně-selekční rovnováha" = "mut_sel_eq",
                      "Migrace: Homogenizace (Pevnina -> Ostrov)" = "migration_eq"
                    )),
        hr(),
        
        h4("1. Selekce (Fitness W)"),
        fluidRow(
          column(4, numericInput("w_AA", "W_AA", value = 1.0, min = 0, max = 1, step = 0.05)),
          column(4, numericInput("w_Aa", "W_Aa", value = 1.0, min = 0, max = 1, step = 0.05)),
          column(4, numericInput("w_aa", "W_aa", value = 0.2, min = 0, max = 1, step = 0.05))
        ),
        
        h4("2. Mutace (Intenzita)"),
        fluidRow(
          column(6, numericInput("u_mut", "u (A -> a)", value = 0.0001, min = 0, max = 0.1, step = 0.0001)),
          column(6, numericInput("v_mut", "v (a -> A)", value = 0.00001, min = 0, max = 0.1, step = 0.0001))
        ),
        
        h4("3. Migrace (Pevnina -> Ostrov)"),
        fluidRow(
          column(6, numericInput("m_mig", "m (Podíl m)", value = 0.0, min = 0, max = 0.5, step = 0.01)),
          column(6, numericInput("q_m", "q_m (Pevnina)", value = 0.5, min = 0, max = 1, step = 0.05))
        ),
        
        hr(),
        h4("Parametry simulace"),
        sliderInput("q0", "Počáteční frekvence alely a (q0):", min = 0.01, max = 0.99, value = 0.5, step = 0.01),
        numericInput("generace", "Počet generací (t):", value = 100, min = 10, max = 1000, step = 10)
      ),
      
      mainPanel(
        width = 8,
        tabsetPanel(
          tabPanel("Trajektorie v čase", 
                   plotOutput("plot_trajektorie", height = "350px"),
                   plotOutput("plot_fitness", height = "250px")),
          tabPanel("Rychlost změny (Delta q)", 
                   plotOutput("plot_delta_q", height = "450px")),
          tabPanel("Tabulka výsledků", 
                   br(),
                   downloadButton("download_data", "Stáhnout data (CSV)"),
                   br(), br(),
                   tableOutput("tabulka_vystup"))
        )
      )
    )
  )
)

# -------------------------------------------------------------------
# 3. SERVEROVÁ LOGIKA
# -------------------------------------------------------------------
server <- function(input, output, session) {
  
  # -----------------------------------------------------------------
  # LOGIKA PRO ZÁLOŽKU 1: DYNAMICKÝ PUNNETTŮV ČTVEREC
  # -----------------------------------------------------------------
  
  output$punnett_q_text <- renderUI({
    q <- round(1 - input$punnett_p, 4)
    HTML(paste0("<b>Frekvence alely a (q = 1 - p):</b> <span style='color: #0072B2; font-weight: bold;'>", q, "</span>"))
  })
  
  output$punnett_summary_html <- renderUI({
    p <- input$punnett_p
    q <- 1 - p
    N <- input$punnett_N
    
    freq_AA <- p^2
    freq_Aa <- 2 * p * q
    freq_aa <- q^2
    
    cnt_AA <- round(N * freq_AA)
    cnt_Aa <- round(N * freq_Aa)
    cnt_aa <- round(N * freq_aa)
    
    HTML(paste0(
      "<b>AA (homozygot):</b> ", round(freq_AA, 4), " (", round(freq_AA * 100, 1), " %) | <b>", cnt_AA, " ind.</b><br/>",
      "<b>Aa (heterozygot):</b> ", round(freq_Aa, 4), " (", round(freq_Aa * 100, 1), " %) | <b>", cnt_Aa, " ind.</b><br/>",
      "<b>aa (homozygot):</b> ", round(freq_aa, 4), " (", round(freq_aa * 100, 1), " %) | <b>", cnt_aa, " ind.</b><br/>",
      "<hr style='margin: 5px 0;'/>",
      "<b>Soustava HWE:</b> p² + 2pq + q² = ", round(freq_AA + freq_Aa + freq_aa, 2)
    ))
  })
  
  output$plot_punnett_square <- renderPlot({
    p <- input$punnett_p
    q <- 1 - p
    N <- input$punnett_N
    
    df_rects <- data.frame(
      xmin = c(0, p, 0, p),
      xmax = c(p, 1, p, 1),
      ymin = c(1 - p, 1 - p, 0, 0),
      ymax = c(1, 1, 1 - p, 1 - p),
      genotype = factor(c("AA", "Aa", "Aa (aA)", "aa"), levels = c("AA", "Aa", "Aa (aA)", "aa")),
      label = c(
        paste0("AA\n(p²)\n", round(p^2, 4), "\n", round(p^2 * 100, 1), "% | ", round(N * p^2), " ind."),
        paste0("Aa\n(p·q)\n", round(p * q, 4), "\n", round(p * q * 100, 1), "% | ", round(N * p * q), " ind."),
        paste0("Aa\n(q·p)\n", round(q * p, 4), "\n", round(q * p * 100, 1), "% | ", round(N * q * p), " ind."),
        paste0("aa\n(q²)\n", round(q^2, 4), "\n", round(q^2 * 100, 1), "% | ", round(N * q^2), " ind.")
      )
    )
    
    df_rects$xmid <- (df_rects$xmin + df_rects$xmax) / 2
    df_rects$ymid <- (df_rects$ymin + df_rects$ymax) / 2
    
    ggplot(df_rects) +
      geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = genotype),
                color = "white", size = 1.5, alpha = 0.85) +
      geom_text(aes(x = xmid, y = ymid, label = label), size = 4.5, fontface = "bold", color = "white") +
      scale_fill_manual(values = c("AA" = "#D55E00", "Aa" = "#0072B2", "Aa (aA)" = "#56B4E9", "aa" = "#009E73")) +
      scale_x_continuous(breaks = c(p / 2, p + q / 2),
                         labels = c(paste0("A (p = ", round(p, 2), ")"), paste0("a (q = ", round(q, 2), ")")),
                         position = "top") +
      scale_y_continuous(breaks = c(1 - p / 2, (1 - p) / 2),
                         labels = c(paste0("A (p = ", round(p, 2), ")"), paste0("a (q = ", round(q, 2), ")"))) +
      labs(title = "Punnettův čtverec - Změna plochy dle frekvencí genotypů potomků",
           subtitle = paste0("Celková plocha (p + q)² = ", round(p^2 + 2*p*q + q^2, 2), " | Součet heterozygotů (2pq) = ", round(2*p*q, 4)),
           x = "Mateřské gamety (vajíčka)", y = "Ocovské gamety (spermie)") +
      theme_minimal(base_size = 14) +
      theme(
        panel.grid = element_blank(),
        axis.text = element_text(size = 13, face = "bold", color = "black"),
        axis.title = element_text(size = 12, face = "bold", color = "gray20"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "gray30"),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        legend.position = "none"
      )
  })
  
  output$plot_punnett_bar <- renderPlot({
    p <- input$punnett_p
    q <- 1 - p
    N <- input$punnett_N
    
    df_bar <- data.frame(
      Genotyp = factor(c("AA (p²)", "Aa (2pq)", "aa (q²)"), levels = c("AA (p²)", "Aa (2pq)", "aa (q²)")),
      Frekvence = c(p^2, 2 * p * q, q^2),
      Pocet = c(round(N * p^2), round(N * 2 * p * q), round(N * q^2))
    )
    
    ggplot(df_bar, aes(x = Genotyp, y = Frekvence, fill = Genotyp)) +
      geom_bar(stat = "identity", width = 0.5, alpha = 0.9) +
      geom_text(aes(label = paste0(round(Frekvence * 100, 1), "%\n(", Pocet, " ind.)")),
                vjust = -0.2, fontface = "bold", size = 4) +
      scale_fill_manual(values = c("AA (p²)" = "#D55E00", "Aa (2pq)" = "#0072B2", "aa (q²)" = "#009E73")) +
      ylim(0, max(df_bar$Frekvence) * 1.25) +
      labs(title = "Souhrnná distribuční struktura genotypů v populaci",
           x = "Genotyp", y = "Relativní frekvence") +
      theme_minimal(base_size = 13) +
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5)
      )
  })
  
  # -----------------------------------------------------------------
  # LOGIKA PRO ZÁLOŽKU 2: HWE (2 ALELY)
  # -----------------------------------------------------------------
  
  hwe_data <- reactive({
    input$btn_resample
    
    if (input$hwe_mode == "theory") {
      p <- input$hwe_p
      q <- 1 - p
      N <- 1000
      n_AA <- round(N * p^2)
      n_Aa <- round(N * 2 * p * q)
      n_aa <- N - n_AA - n_Aa
      return(vypocti_hwe_2alely(n_AA, n_Aa, n_aa))
    } else if (input$hwe_mode == "sample") {
      p <- input$hwe_p
      N <- input$hwe_N
      q <- 1 - p
      genotypes <- rbinom(N, size = 2, prob = q)
      n_AA <- sum(genotypes == 0)
      n_Aa <- sum(genotypes == 1)
      n_aa <- sum(genotypes == 2)
      return(vypocti_hwe_2alely(n_AA, n_Aa, n_aa))
    } else {
      n_AA <- ifelse(is.na(input$manual_AA), 0, max(0, input$manual_AA))
      n_Aa <- ifelse(is.na(input$manual_Aa), 0, max(0, input$manual_Aa))
      n_aa <- ifelse(is.na(input$manual_aa), 0, max(0, input$manual_aa))
      return(vypocti_hwe_2alely(n_AA, n_Aa, n_aa))
    }
  })
  
  output$hwe_summary_text <- renderUI({
    res <- hwe_data()
    if (is.null(res)) return("Zadejte platná data.")
    
    HTML(paste0(
      "<b>Počet jedinců (N):</b> ", res$N, "<br/>",
      "<b>Frekvence alely A (p):</b> ", round(res$p, 4), "<br/>",
      "<b>Frekvence alely a (q):</b> ", round(res$q, 4), "<br/>",
      "<b>Očekávaná HWE proporce:</b><br/>",
      "p² (AA) = ", round(res$p^2, 4), " | 2pq (Aa) = ", round(2*res$p*res$q, 4), " | q² (aa) = ", round(res$q^2, 4)
    ))
  })
  
  output$plot_hwe_curves <- renderPlot({
    res <- hwe_data()
    q_seq <- seq(0, 1, length.out = 200)
    p_seq <- 1 - q_seq
    
    df_curves <- data.frame(
      q = rep(q_seq, 3),
      Frekvence = c(p_seq^2, 2 * p_seq * q_seq, q_seq^2),
      Genotyp = factor(rep(c("AA (p²)", "Aa (2pq)", "aa (q²)"), each = 200),
                       levels = c("AA (p²)", "Aa (2pq)", "aa (q²)"))
    )
    
    p <- ggplot(df_curves, aes(x = q, y = Frekvence, color = Genotyp)) +
      geom_line(size = 1.2) +
      scale_color_manual(values = c("AA (p²)" = "#D55E00", "Aa (2pq)" = "#0072B2", "aa (q²)" = "#009E73")) +
      labs(title = "Teoretické křivky HWE a pozice populace",
           x = "Frekvence alely a (q)", y = "Genotypová frekvence") +
      ylim(0, 1) +
      theme_minimal(base_size = 13) +
      theme(panel.background = element_rect(fill = "white", color = NA),
            plot.background = element_rect(fill = "white", color = NA),
            legend.position = "top")
    
    if (!is.null(res)) {
      df_obs <- data.frame(
        q = rep(res$q, 3),
        Frekvence = c(unname(res$obs_freq["AA"]), unname(res$obs_freq["Aa"]), unname(res$obs_freq["aa"])),
        Genotyp = factor(c("AA (p²)", "Aa (2pq)", "aa (q²)"), levels = c("AA (p²)", "Aa (2pq)", "aa (q²)"))
      )
      p <- p + geom_point(data = df_obs, aes(x = q, y = Frekvence, color = Genotyp), size = 4, shape = 18)
    }
    
    p
  })
  
  output$plot_hwe_counts <- renderPlot({
    res <- hwe_data()
    if (is.null(res)) return(NULL)
    
    df_bar <- data.frame(
      Genotyp = factor(rep(c("AA", "Aa", "aa"), 2), levels = c("AA", "Aa", "aa")),
      Typ = factor(rep(c("Pozorované (P)", "Očekávané (O)"), each = 3), levels = c("Pozorované (P)", "Očekávané (O)")),
      Pocet = c(unname(res$obs["AA"]), unname(res$obs["Aa"]), unname(res$obs["aa"]),
                unname(res$exp_counts["AA"]), unname(res$exp_counts["Aa"]), unname(res$exp_counts["aa"]))
    )
    
    ggplot(df_bar, aes(x = Genotyp, y = Pocet, fill = Typ)) +
      geom_bar(stat = "identity", position = "dodge", width = 0.7) +
      scale_fill_manual(values = c("Pozorované (P)" = "#4C72B0", "Očekávané (O)" = "#DD8452")) +
      labs(title = "Porovnání Pozorovaných (P) a Očekávaných (O) četností",
           x = "Genotyp", y = "Počet jedinců", fill = "") +
      theme_minimal(base_size = 13) +
      theme(panel.background = element_rect(fill = "white", color = NA),
            plot.background = element_rect(fill = "white", color = NA),
            legend.position = "top")
  })
  
  output$tabulka_hwe_stats <- renderTable({
    res <- hwe_data()
    if (is.null(res)) return(NULL)
    
    obs_AA <- unname(res$obs["AA"])
    obs_Aa <- unname(res$obs["Aa"])
    obs_aa <- unname(res$obs["aa"])
    
    exp_AA <- unname(res$exp_counts["AA"])
    exp_Aa <- unname(res$exp_counts["Aa"])
    exp_aa <- unname(res$exp_counts["aa"])
    
    chi2_AA <- (obs_AA - exp_AA)^2 / ifelse(exp_AA == 0, 1e-6, exp_AA)
    chi2_Aa <- (obs_Aa - exp_Aa)^2 / ifelse(exp_Aa == 0, 1e-6, exp_Aa)
    chi2_aa <- (obs_aa - exp_aa)^2 / ifelse(exp_aa == 0, 1e-6, exp_aa)
    
    data.frame(
      Genotyp = c("AA", "Aa", "aa", "Celkem (N)"),
      Pozorovane_P = c(obs_AA, obs_Aa, obs_aa, res$N),
      Ocekavane_O = c(round(exp_AA, 2), round(exp_Aa, 2), round(exp_aa, 2), res$N),
      Podil_Chi2 = c(round(chi2_AA, 4), round(chi2_Aa, 4), round(chi2_aa, 4), round(res$chi2, 4))
    )
  }, digits = 2)
  
  output$hwe_test_result <- renderUI({
    res <- hwe_data()
    if (is.null(res)) return(NULL)
    
    sig_005 <- res$p_val < 0.05
    status_text <- if (sig_005) {
      paste0("<b>ZÁVĚR TESTU: Zamítáme hypotézu H₀ (p = ", format.pval(res$p_val, digits = 4), " < 0.05).</b><br/>",
             "Pozorované četnosti genotypů se statisticky významně liší od očekávání HWE. Populace NENÍ v Hardy-Weinbergově rovnováze.")
    } else {
      paste0("<b>ZÁVĚR TESTU: Nezamítáme hypotézu H₀ (p = ", format.pval(res$p_val, digits = 4), " ≥ 0.05).</b><br/>",
             "Pozorované a očekávané četnosti jsou v souladu. Populace JE v Hardy-Weinbergově rovnováze.")
    }
    
    HTML(paste0(
      "<b>Testová statistika (χ²):</b> ", round(res$chi2, 4), " | <b>Stupně volnosti (df):</b> 1 | <b>p-hodnota:</b> ", format.pval(res$p_val, digits = 4), "<br/><br/>",
      "<div style='padding: 10px; border-radius: 5px; background-color: ", ifelse(sig_005, "#f8d7da", "#d4edda"), "; color: ", ifelse(sig_005, "#721c24", "#155724"), ";'>",
      status_text, "</div>"
    ))
  })

  # -----------------------------------------------------------------
  # LOGIKA PRO ZÁLOŽKU 3: HWE (3 ALELY / ABO)
  # -----------------------------------------------------------------
  output$ui_abo_q <- renderUI({
    max_q <- round(1 - input$abo_p - 0.01, 2)
    val_q <- min(0.20, max_q)
    sliderInput("abo_q", "Frekvence alely I^B (q):", min = 0.01, max = max_q, value = val_q, step = 0.01)
  })
  
  output$abo_r_text <- renderUI({
    q_val <- if (!is.null(input$abo_q)) input$abo_q else 0.20
    r_val <- round(1 - input$abo_p - q_val, 4)
    HTML(paste0("<b>Frekvence alely i (r = 1 - p - q):</b> <span style='color: #D55E00; font-weight: bold;'>", r_val, "</span>"))
  })
  
  abo_data <- reactive({
    q_val <- if (!is.null(input$abo_q)) input$abo_q else 0.20
    vypocti_hwe_3alely(p = input$abo_p, q = q_val, N = input$abo_N)
  })
  
  output$plot_3allele_genotypes <- renderPlot({
    res <- abo_data()
    df_gen <- data.frame(
      Genotyp = factor(c("AA (p²)", "BB (q²)", "OO (r²)", "AB (2pq)", "AO (2pr)", "BO (2qr)"),
                       levels = c("AA (p²)", "BB (q²)", "OO (r²)", "AB (2pq)", "AO (2pr)", "BO (2qr)")),
      Frekvence = unname(res$gen_freq),
      Pocet = unname(res$exp_gen_counts)
    )
    
    ggplot(df_gen, aes(x = Genotyp, y = Pocet, fill = Genotyp)) +
      geom_bar(stat = "identity", width = 0.6) +
      scale_fill_brewer(palette = "Set2") +
      labs(title = paste0("Očekávané četnosti 6 genotypů (N = ", res$N, ")"),
           x = "Genotyp", y = "Očekávaný počet jedinců") +
      theme_minimal(base_size = 13) +
      theme(panel.background = element_rect(fill = "white", color = NA),
            plot.background = element_rect(fill = "white", color = NA),
            legend.position = "none")
  })
  
  output$plot_3allele_phenotypes <- renderPlot({
    res <- abo_data()
    df_phen <- data.frame(
      Fenotyp = factor(names(res$phen_freq), levels = names(res$phen_freq)),
      Frekvence = unname(res$phen_freq),
      Pocet = unname(res$exp_phen_counts)
    )
    
    ggplot(df_phen, aes(x = Fenotyp, y = Pocet, fill = Fenotyp)) +
      geom_bar(stat = "identity", width = 0.5) +
      scale_fill_manual(values = c("Skupina A (AA + AO)" = "#4C72B0", 
                                   "Skupina B (BB + BO)" = "#DD8452", 
                                   "Skupina AB (AB)" = "#55A868", 
                                   "Skupina O (OO)" = "#C44E52")) +
      labs(title = "Očekávané četnosti 4 fenotypů (krevní skupiny ABO)",
           x = "Fenotyp", y = "Očekávaný počet jedinců") +
      theme_minimal(base_size = 13) +
      theme(panel.background = element_rect(fill = "white", color = NA),
            plot.background = element_rect(fill = "white", color = NA),
            legend.position = "none")
  })
  
  output$tabulka_3allele_summary <- renderTable({
    res <- abo_data()
    data.frame(
      Kategorie = c("Alela A (p)", "Alela B (q)", "Alela O (r)",
                    "Genotyp AA (p²)", "Genotyp BB (q²)", "Genotyp OO (r²)",
                    "Genotyp AB (2pq)", "Genotyp AO (2pr)", "Genotyp BO (2qr)",
                    "Fenotyp Skupina A (p² + 2pr)", "Fenotyp Skupina B (q² + 2qr)",
                    "Fenotyp Skupina AB (2pq)", "Fenotyp Skupina O (r²)"),
      Frekvence = c(round(res$p, 4), round(res$q, 4), round(res$r, 4),
                    round(unname(res$gen_freq["AA"]), 4), round(unname(res$gen_freq["BB"]), 4), round(unname(res$gen_freq["OO"]), 4),
                    round(unname(res$gen_freq["AB"]), 4), round(unname(res$gen_freq["AO"]), 4), round(unname(res$gen_freq["BO"]), 4),
                    round(unname(res$phen_freq[1]), 4), round(unname(res$phen_freq[2]), 4),
                    round(unname(res$phen_freq[3]), 4), round(unname(res$phen_freq[4]), 4)),
      Ocekavany_Pocet = c("-", "-", "-",
                          round(unname(res$exp_gen_counts["AA"]), 1), round(unname(res$exp_gen_counts["BB"]), 1), round(unname(res$exp_gen_counts["OO"]), 1),
                          round(unname(res$exp_gen_counts["AB"]), 1), round(unname(res$exp_gen_counts["AO"]), 1), round(unname(res$exp_gen_counts["BO"]), 1),
                          round(unname(res$exp_phen_counts[1]), 1), round(unname(res$exp_phen_counts[2]), 1),
                          round(unname(res$exp_phen_counts[3]), 1), round(unname(res$exp_phen_counts[4]), 1))
    )
  })

  # -----------------------------------------------------------------
  # LOGIKA PRO ZÁLOŽKU 4: EVOLUČNÍ SÍLY
  # -----------------------------------------------------------------
  observeEvent(input$preset, {
    if (input$preset == "sel_recessive") {
      updateNumericInput(session, "w_AA", value = 1.0)
      updateNumericInput(session, "w_Aa", value = 1.0)
      updateNumericInput(session, "w_aa", value = 0.0)
      updateNumericInput(session, "u_mut", value = 0)
      updateNumericInput(session, "v_mut", value = 0)
      updateNumericInput(session, "m_mig", value = 0)
    } else if (input$preset == "overdominance") {
      updateNumericInput(session, "w_AA", value = 0.7)
      updateNumericInput(session, "w_Aa", value = 1.0)
      updateNumericInput(session, "w_aa", value = 0.4)
      updateNumericInput(session, "u_mut", value = 0)
      updateNumericInput(session, "v_mut", value = 0)
      updateNumericInput(session, "m_mig", value = 0)
    } else if (input$preset == "mutation_eq") {
      updateNumericInput(session, "w_AA", value = 1.0)
      updateNumericInput(session, "w_Aa", value = 1.0)
      updateNumericInput(session, "w_aa", value = 1.0)
      updateNumericInput(session, "u_mut", value = 0.001)
      updateNumericInput(session, "v_mut", value = 0.0005)
      updateNumericInput(session, "m_mig", value = 0)
    } else if (input$preset == "mut_sel_eq") {
      updateNumericInput(session, "w_AA", value = 1.0)
      updateNumericInput(session, "w_Aa", value = 1.0)
      updateNumericInput(session, "w_aa", value = 0.0)
      updateNumericInput(session, "u_mut", value = 0.0001)
      updateNumericInput(session, "v_mut", value = 0)
      updateNumericInput(session, "m_mig", value = 0)
    } else if (input$preset == "migration_eq") {
      updateNumericInput(session, "w_AA", value = 1.0)
      updateNumericInput(session, "w_Aa", value = 1.0)
      updateNumericInput(session, "w_aa", value = 1.0)
      updateNumericInput(session, "u_mut", value = 0)
      updateNumericInput(session, "v_mut", value = 0)
      updateNumericInput(session, "m_mig", value = 0.05)
      updateNumericInput(session, "q_m", value = 0.8)
    }
  })
  
  data_sim <- reactive({
    simuluj_populaci(
      q0 = input$q0,
      W_AA = input$w_AA,
      W_Aa = input$w_Aa,
      W_aa = input$w_aa,
      u = input$u_mut,
      v = input$v_mut,
      m = input$m_mig,
      q_m = input$q_m,
      generace = input$generace
    )
  })
  
  output$plot_trajektorie <- renderPlot({
    df <- data_sim()
    ggplot(df, aes(x = Generace)) +
      geom_line(aes(y = q, color = "Frekvence q (a)"), size = 1.2) +
      geom_line(aes(y = p, color = "Frekvence p (A)"), size = 1.2, linetype = "dashed") +
      scale_color_manual(values = c("Frekvence q (a)" = "#D55E00", "Frekvence p (A)" = "#0072B2")) +
      labs(title = "Trajektorie alelových frekvencí v čase",
           x = "Čas (generace t)", y = "Frekvence alel", color = "Alela") +
      ylim(0, 1) +
      theme_minimal(base_size = 14) +
      theme(panel.background = element_rect(fill = "white", color = NA),
            plot.background = element_rect(fill = "white", color = NA))
  })
  
  output$plot_fitness <- renderPlot({
    df <- data_sim()
    ggplot(df, aes(x = Generace, y = W_bar)) +
      geom_line(color = "#009E73", size = 1.2) +
      labs(title = "Průměrná adaptivní hodnota populace (W_bar)",
           x = "Čas (generace t)", y = "Průměrné W") +
      theme_minimal(base_size = 13) +
      theme(panel.background = element_rect(fill = "white", color = NA),
            plot.background = element_rect(fill = "white", color = NA))
  })
  
  output$plot_delta_q <- renderPlot({
    q_vals <- seq(0, 1, length.out = 200)
    dq_vals <- sapply(q_vals, function(qv) {
      vypocti_dalsi_generaci(qv, input$w_AA, input$w_Aa, input$w_aa, 
                             input$u_mut, input$v_mut, input$m_mig, input$q_m)$delta_q
    })
    
    df_dq <- data.frame(q = q_vals, delta_q = dq_vals)
    
    ggplot(df_dq, aes(x = q, y = delta_q)) +
      geom_line(color = "#E69F00", size = 1.3) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
      labs(title = "Rychlost evoluce: Změna Delta q v závislosti na aktuální frekvenci q",
           x = "Frekvence alely q", y = "Změna za 1 generaci (Delta q)") +
      theme_minimal(base_size = 14) +
      theme(panel.background = element_rect(fill = "white", color = NA),
            plot.background = element_rect(fill = "white", color = NA))
  })
  
  output$tabulka_vystup <- renderTable({
    df <- data_sim()
    step <- max(1, floor(nrow(df) / 15))
    df[seq(1, nrow(df), by = step), ]
  }, digits = 4)
  
  output$download_data <- downloadHandler(
    filename = function() { paste("simulace_popgen_", Sys.Date(), ".csv", sep = "") },
    content = function(file) { write.csv(data_sim(), file, row.names = FALSE) }
  )
}

# -------------------------------------------------------------------
# 4. SPUŠTĚNÍ APLIKACE
# -------------------------------------------------------------------
shinyApp(ui = ui, server = server)
