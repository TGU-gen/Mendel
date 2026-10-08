# ===================================================================
# INTERAKTIVNÍ SHINY APLIKACE: POPULAČNÍ GENETIKA, HWE A EVOLUČNÍ SÍLY
# (Punnettův čtverec, HWE 2 alely, HWE 3 alely / ABO, Evoluční síly, Genetický drift, Inbreeding, Drift vs. Selekce)
# vazba genů a X chromozom
# Autor: Tomáš Urban, urban@mendelu.cz, UMFGZ, AF  MENDELU
# ===================================================================

options(encoding = "UTF-8")

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
# 1. MATEMATICKÉ FUNKCE PRO HWE, DRIFT, INBREEDING & EVOLUČNÍ SÍLY
# -------------------------------------------------------------------

# A) Výpočet HWE pro 2 alely (n_AA, n_Aa, n_aa)
vypocti_hwe_2alely <- function(n_AA, n_Aa, n_aa, perfect_hwe = FALSE) {
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
  
  # Pokud je aktivováno perfect_hwe = TRUE, pozorovaná data nastavíme na teoretické HWE počty
  if(perfect_hwe){
    chi2_val <- 0
    p_val <- 1
  }
  
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
    IᴬIᴬ = p^2,
    IᴮIᴮ = q^2,
    ii = r^2,
    IᴬIᴮ = 2 * p * q,
    Iᴬi = 2 * p * r,
    Iᴮi = 2 * q * r
  )
  
  phen_freq <- c(
    "Skupina A (IᴬIᴬ + Iᴬi)" = p^2 + 2 * p * r,
    "Skupina B (IᴮIᴮ + Iᴮi)" = q^2 + 2 * q * r,
    "Skupina AB (IᴬIᴮ)"    = 2 * p * q,
    "Skupina O (ii)"      = r^2
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

# C) Simulace Genetického Driftu (Wright-Fisherův model)
simuluj_drift <- function(N, q0, K = 10, generace = 100) {
  N <- max(1, round(as.numeric(N)))
  q0 <- as.numeric(q0)
  K <- max(1, round(as.numeric(K)))
  generace <- max(1, round(as.numeric(generace)))
  
  alely_2N <- 2 * N
  q_mat <- matrix(NA_real_, nrow = generace + 1, ncol = K)
  q_mat[1, ] <- q0
  
  for (t in seq_len(generace)) {
    q_mat[t + 1, ] <- rbinom(K, size = alely_2N, prob = q_mat[t, ]) / alely_2N
  }
  
  df_lines <- data.frame(
    Generace = rep(0:generace, times = K),
    q = as.vector(q_mat),
    Populace = factor(rep(seq_len(K), each = generace + 1))
  )
  
  H0 <- 2 * q0 * (1 - q0)
  mean_H_obs <- rowMeans(2 * q_mat * (1 - q_mat))
  mean_q_obs <- rowMeans(q_mat)
  H_teor <- H0 * ((1 - 1 / alely_2N)^(0:generace))
  
  df_het <- data.frame(
    Generace = 0:generace,
    H_obs = mean_H_obs,
    H_teor = H_teor,
    mean_q = mean_q_obs
  )
  
  last_q <- q_mat[generace + 1, ]
  fixed_a <- sum(last_q == 1)
  lost_a <- sum(last_q == 0)
  poly <- sum(last_q > 0 & last_q < 1)
  
  return(list(
    df_lines = df_lines,
    df_het = df_het,
    last_q = last_q,
    fixed_a = fixed_a,
    lost_a = lost_a,
    poly = poly,
    N = N,
    K = K,
    generace = generace,
    q0 = q0
  ))
}

# D) Výpočet Inbreedingu (Koeficient F)
vypocti_inbreeding <- function(p, F_coeff, N = 1000) {
  p <- as.numeric(p)
  q <- 1 - p
  F_coeff <- as.numeric(F_coeff)
  N <- as.numeric(N)
  
  freq_0 <- c(AA = p^2, Aa = 2 * p * q, aa = q^2)
  freq_F <- c(AA = p^2 + F_coeff * p * q, Aa = 2 * p * q * (1 - F_coeff), aa = q^2 + F_coeff * p * q)
  
  safe_counts <- function(freq, N) {
    counts <- floor(freq * N)
    
    counts[which.max(freq)] <-
      counts[which.max(freq)] +
      (N - sum(counts))
    
    counts
  }
  
  counts_F <- safe_counts(freq_F, N)
  counts_0 <- safe_counts(freq_0, N)
  
  risk_ratio <- if (freq_0["aa"] > 0) unname(freq_F["aa"] / freq_0["aa"]) else NA_real_
  
  return(list(
    p = p, q = q, F = F_coeff, N = N,
    freq_0 = freq_0,
    freq_F = freq_F,
    counts_0 = counts_0,
    counts_F = counts_F,
    risk_ratio = risk_ratio
  ))
}

# Inbreeding v panmiktické populaci konečné velikosti N po t generacích
#   F_t = 1 - (1 - 1/(2N))^t
inbreeding_Ft <- function(N, t) {
  1 - (1 - 1 / (2 * N))^t
}

# Počet generací, za které F dosáhne cílové hodnoty
#   t = ln(1 - F) / ln(1 - 1/(2N))
inbreeding_t_do_F <- function(N, F_cil) {
  if (F_cil <= 0) return(0)
  if (F_cil >= 1) return(Inf)
  log(1 - F_cil) / log(1 - 1 / (2 * N))
}

# G) Drift vs. selekce: Wright-Fisherův model se selekcí
#    Fitness: W_AA = 1, W_Aa = 1 + h*s, W_aa = 1 + s (alela a je výhodná)
#    Každá generace: 1) selekce -> q_sel, 2) vzorkování gamet X ~ Bin(2Ne, q_sel), q' = X / 2Ne
simuluj_drift_selekce <- function(Ne, s, h, q0, K, generace) {
  Ne <- round(Ne)
  n_gamet <- 2 * Ne
  W_AA <- 1; W_Aa <- 1 + h * s; W_aa <- 1 + s
  
  selekce_krok <- function(q) {
    p <- 1 - q
    W_bar <- p^2 * W_AA + 2 * p * q * W_Aa + q^2 * W_aa
    pmin(1, pmax(0, (p * q * W_Aa + q^2 * W_aa) / W_bar))
  }
  
  # stochastické linie (K replikátů)
  Q <- matrix(NA_real_, nrow = generace + 1, ncol = K)
  Q[1, ] <- q0
  for (t in seq_len(generace)) {
    q_sel <- selekce_krok(Q[t, ])
    Q[t + 1, ] <- rbinom(K, size = n_gamet, prob = q_sel) / n_gamet
  }
  
  # deterministická trajektorie (selekce bez driftu)
  q_det <- numeric(generace + 1)
  q_det[1] <- q0
  for (t in seq_len(generace)) q_det[t + 1] <- selekce_krok(q_det[t])
  
  q_konec <- Q[generace + 1, ]
  list(
    Ne = Ne, s = s, h = h, q0 = q0, K = K, generace = generace,
    x = Ne * s,                          # Ne*s
    Q = Q, q_det = q_det, q_mean = rowMeans(Q),
    fixed = sum(q_konec == 1), lost = sum(q_konec == 0),
    poly = sum(q_konec > 0 & q_konec < 1)
  )
}

# Pravděpodobnost fixace (Kimura, aditivní selekce h = 1/2), x = Ne*s
#   P_fix = (1 - exp(-2*x*q0)) / (1 - exp(-2*x));  pro x -> 0 je P_fix = q0 (neutralita)
ds_pfix <- function(x, q0) {
  ifelse(x < 1e-8, q0, expm1(-2 * x * q0) / expm1(-2 * x))
}

# Režim podle Ne*s
ds_rezim <- function(x) {
  if (x < 1) {
    list(text = "Drift dominuje (Ne·s < 1)", barva = "#0072B2")
  } else if (x <= 10) {
    list(text = "Přechodová zóna (1 ≤ Ne·s ≤ 10)", barva = "#E69F00")
  } else {
    list(text = "Selekce dominuje (Ne·s > 10)", barva = "#009E73")
  }
}

# E) Výpočet jednoho evolučního kroku (Selekce, Mutace, Migrace pro 2 alely)
vypocti_dalsi_generaci <- function(q, W_AA, W_Aa, W_aa, u, v, m, q_m) {
  p <- 1 - q
  
  # Selekce
  W_bar <- (p^2 * W_AA) + (2 * p * q * W_Aa) + (q^2 * W_aa)
  q_sel <- if (W_bar > 0) (p * q * W_Aa + q^2 * W_aa) / W_bar else q
  
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

# F) Simulace evoluce v čase
simuluj_populaci <- function(q0, W_AA, W_Aa, W_aa, u, v, m, q_m, generace) {
  q_vec <- numeric(generace + 1)
  W_vec <- numeric(generace + 1)
  q_vec[1] <- q0
  
  for (t in seq_len(generace)) {
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
# POMOCNÉ FUNKCE PRO ZOBRAZENÍ VZORCŮ (Čisté HTML/CSS)
# -------------------------------------------------------------------
vzorce_css <- "
  .vzorce-box {
    background: #f7f9fb;
    border: 1px solid #d5dde5;
    border-left: 4px solid #2c7be5;
    border-radius: 4px;
    padding: 8px 12px;
    margin-top: 15px;
  }
  .vzorce-box summary {
    font-weight: bold;
    cursor: pointer;
    outline: none;
    color: #1a2b4c;
  }
  .vzorec-radek {
    margin: 8px 0 4px 0;
  }
  .vzorec-popis {
    display: block;
    font-size: 0.82em;
    color: #6c757d;
  }
  .vzorec {
    font-family: 'Cambria Math', 'STIX Two Math', 'Times New Roman', serif;
    font-size: 1.12em;
    line-height: 1.9;
    color: #0f172a;
  }
  .zlomek {
    display: inline-block;
    vertical-align: middle;
    text-align: center;
    margin: 0 3px;
  }
  .zlomek .cit {
    display: block;
    padding: 0 4px;
    border-bottom: 1px solid #333;
    line-height: 1.35;
  }
  .zlomek .jmen {
    display: block;
    padding: 0 4px;
    line-height: 1.35;
  }
"

zl <- function(cit, jmen) {
  paste0("<span class='zlomek'><span class='cit'>", cit, "</span><span class='jmen'>", jmen, "</span></span>")
}

vz <- function(popis, vzorec) {
  tags$div(class = "vzorec-radek",
    tags$span(class = "vzorec-popis", HTML(popis)),
    tags$div(class = "vzorec", HTML(vzorec))
  )
}

vzorce_box <- function(...) {
  tags$details(class = "vzorce-box", open = NA,
    tags$summary("Použité vzorce a teoretický základ"),
    ...
  )
}

# -------------------------------------------------------------------
# LOGO V PRAVÉM HORNÍM ROHU (s odkazem na jiný web)
# -------------------------------------------------------------------
# Soubor s logem uložte do podsložky "www" vedle tohoto skriptu (např. www/logo.png).
# Pokud soubor neexistuje, zobrazí se místo loga textový odkaz (logo_text).
logo_soubor <- "logo.png"                # název souboru ve složce www
logo_odkaz  <- "https://www.mendelu.cz"  # adresa, na kterou logo odkazuje
logo_text   <- "MENDELU"                 # záložní text + popisek při najetí myší
logo_vyska  <- 36                        # výška loga v pixelech

vytvor_logo_html <- function() {
  obsah <- if (file.exists(file.path("www", logo_soubor))) {
    tags$img(src = logo_soubor, alt = logo_text, style = paste0("height: ", logo_vyska, "px;"))
  } else {
    logo_text
  }
  as.character(tags$a(
    href = logo_odkaz, target = "_blank", rel = "noopener noreferrer",
    class = "navbar-logo", title = logo_text, obsah
  ))
}

logo_css <- "
  .navbar-logo {
    position: absolute;
    right: 15px;
    top: 0;
    height: 56px;
    display: flex;
    align-items: center;
    z-index: 1100;
    color: #ffffff;
    font-weight: bold;
    text-decoration: none;
  }
  .navbar-logo:hover, .navbar-logo:focus { opacity: 0.85; color: #ffffff; text-decoration: none; }
  @media (max-width: 767px) { .navbar-logo { right: 60px; } }
"

# JavaScript vloží logo do horní lišty (navbar) po načtení stránky
logo_js <- paste0(
  "var logoHtml = ", jsonlite::toJSON(vytvor_logo_html(), auto_unbox = TRUE), ";\n",
  "$(function() { $('.navbar > .container-fluid').first().append(logoHtml); });"
)

# -------------------------------------------------------------------
# 2. UŽIVATELSKÉ ROZHRANÍ (UI - NAVBAR)
# -------------------------------------------------------------------
ui <- navbarPage(
  theme = app_theme,
  header = tags$head(
    tags$style(HTML(vzorce_css)),
    tags$style(HTML(logo_css)),
    tags$script(HTML(logo_js))
  ),
  title = "Genetika populací: HWE & Evoluční síly",
  footer = tags$div(
    class = "text-muted",
    style = "text-align: center; padding: 12px; margin-top: 20px; border-top: 1px solid #ddd; font-size: 0.9em;",
    HTML(paste0("Autor: <b>Tomáš Urban</b>, UMFGZ, AF MENDELU © ", format(Sys.Date(), "%Y")))
  ),
  
  # =================================================================
  # ZÁLOŽKA 1: DYNAMICKÝ PUNNETTŮV ČTVEREC (PLOCHY GENOTYPŮ)
  # =================================================================
  tabPanel(
    "1. Punnettův čtverec (za HWE)",
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
        ),
        vzorce_box(
          vz("Součet frekvencí alel", "<i>p</i> + <i>q</i> = 1"),
          vz("Plochy v Punnettově čtverci = genotypové frekvence", "AA = <i>p</i><sup>2</sup> &nbsp;&nbsp; Aa = 2<i>p</i><i>q</i> &nbsp;&nbsp; aa = <i>q</i><sup>2</sup>"),
          vz("Kontrola součtu", "(<i>p</i> + <i>q</i>)<sup>2</sup> = <i>p</i><sup>2</sup> + 2<i>p</i><i>q</i> + <i>q</i><sup>2</sup> = 1"),
          vz("Očekávané počty jedinců", "<i>N</i>·<i>p</i><sup>2</sup> &nbsp;|&nbsp; <i>N</i>·2<i>p</i><i>q</i> &nbsp;|&nbsp; <i>N</i>·<i>q</i><sup>2</sup>")
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
    "2. HWE model a testování",
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
        ),
        vzorce_box(
          vz("Frekvence alel z pozorovaných počtů",
             paste0("<i>p</i> = ", zl("2·<i>n</i><sub>AA</sub> + <i>n</i><sub>Aa</sub>", "2<i>N</i>"), " &nbsp;&nbsp; <i>q</i> = 1 − <i>p</i>")),
          vz("Očekávané počty podle HWE",
             "<i>E</i><sub>AA</sub> = <i>N</i><i>p</i><sup>2</sup> &nbsp; <i>E</i><sub>Aa</sub> = 2<i>N</i><i>p</i><i>q</i> &nbsp; <i>E</i><sub>aa</sub> = <i>N</i><i>q</i><sup>2</sup>"),
          vz("Test dobré shody (df = 1)",
             paste0("χ<sup>2</sup> = Σ ", zl("(<i>O</i> − <i>E</i>)<sup>2</sup>", "<i>E</i>")))
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
        sliderInput("abo_p", "Frekvence alely Iᴬ (p):", min = 0.01, max = 0.90, value = 0.30, step = 0.01),
        uiOutput("ui_abo_q"),
        wellPanel(
          htmlOutput("abo_r_text")
        ),
        numericInput("abo_N", "Velikost populace (N):", value = 1000, min = 100, max = 50000, step = 100),
        vzorce_box(
          vz("Frekvence třetí alely (i)", "<i>r</i> = 1 − <i>p</i> − <i>q</i>"),
          vz("Genotypové frekvence", "IᴬIᴬ = <i>p</i><sup>2</sup> &nbsp; IᴮIᴮ = <i>q</i><sup>2</sup> &nbsp; ii = <i>r</i><sup>2</sup><br>IᴬIᴮ = 2<i>p</i><i>q</i> &nbsp; Iᴬi = 2<i>p</i><i>r</i> &nbsp; Iᴮi = 2<i>q</i><i>r</i>"),
          vz("Fenotypové frekvence (krevní skupiny)", "A = <i>p</i><sup>2</sup> + 2<i>p</i><i>r</i><br>B = <i>q</i><sup>2</sup> + 2<i>q</i><i>r</i><br>AB = 2<i>p</i><i>q</i><br>O = <i>r</i><sup>2</sup>"),
          vz("Kontrola součtu", "(<i>p</i> + <i>q</i> + <i>r</i>)<sup>2</sup> = 1")
        )
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
        numericInput("generace", "Počet generací (t):", value = 100, min = 10, max = 1000, step = 10),
        vzorce_box(
          vz("Průměrná zdatnost populace",
             "<i>W&#772;</i> = <i>p</i><sup>2</sup><i>W</i><sub>AA</sub> + 2<i>p</i><i>q</i><i>W</i><sub>Aa</sub> + <i>q</i><sup>2</sup><i>W</i><sub>aa</sub>"),
          vz("1. Selekce",
             paste0("<i>q</i>′ = ", zl("<i>p</i><i>q</i><i>W</i><sub>Aa</sub> + <i>q</i><sup>2</sup><i>W</i><sub>aa</sub>", "<i>W&#772;</i>"))),
          vz("2. Mutace (u: A → a, v: a → A)",
             "<i>q</i>″ = <i>q</i>′(1 − <i>v</i>) + (1 − <i>q</i>′)<i>u</i>"),
          vz("3. Migrace (podíl migrantů m)",
             "<i>q</i><sub><i>t</i>+1</sub> = (1 − <i>m</i>)<i>q</i>″ + <i>m</i>·<i>q</i><sub><i>m</i></sub>"),
          vz("Změna frekvence za generaci", "Δ<i>q</i> = <i>q</i><sub><i>t</i>+1</sub> − <i>q</i><sub><i>t</i></sub>"),
          vz("Rovnovážné stavy (při působení jediné síly)",
             paste0("mutace: &nbsp;<i>q&#770;</i> = ", zl("<i>u</i>", "<i>u</i> + <i>v</i>"),
                    "<br>migrace: &nbsp;<i>q&#770;</i> = <i>q</i><sub><i>m</i></sub>",
                    "<br>naddominance: &nbsp;<i>q&#770;</i> = ", zl("<i>W</i><sub>Aa</sub> − <i>W</i><sub>AA</sub>", "2<i>W</i><sub>Aa</sub> − <i>W</i><sub>AA</sub> − <i>W</i><sub>aa</sub>"),
                    "<br>mutačně-selekční (aproximace pro plně recesivní škodlivou alelu a malé hodnoty mutační intenzity): &nbsp;<i>q&#770;</i> ≈ √(", zl("<i>u</i>", "<i>s</i>"), "), &nbsp;<i>s</i> = 1 − <i>W</i><sub>aa</sub>"))
        )
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
  ),
  
  # =================================================================
  # ZÁLOŽKA 5: GENETICKÝ DRIFT (WRIGHT-FISHERŮV MODEL)
  # =================================================================
  tabPanel(
    "5. Genetický drift",
    sidebarLayout(
      sidebarPanel(
        width = 4,
        h3("Genetický drift (Wright-Fisher)"),
        helpText("Stochastické kolísání alelových frekvencí v konečných populacích v důsledku náhodného výběru gamet."),
        hr(),
        numericInput("drift_N", "Velikost populace N (diploidi):", value = 50, min = 5, max = 2000, step = 5),
        sliderInput("drift_q0", "Počáteční frekvence alely a (q₀):", min = 0.01, max = 0.99, value = 0.50, step = 0.01),
        sliderInput("drift_K", "Počet replikátních populací (linií):", min = 1, max = 30, value = 10, step = 1),
        numericInput("drift_G", "Počet generací (t):", value = 100, min = 10, max = 500, step = 10),
        hr(),
        actionButton("btn_drift_sim", "Spustit novou simulaci driftu", class = "btn-primary"),
        br(), br(),
        wellPanel(
          h4("Výsledek po t generacích:"),
          htmlOutput("drift_summary_html")
        ),
        vzorce_box(
          vz("Vzorkování gamet v generaci <i>t</i> → <i>t</i>+1", "<i>X</i> ~ Bin(2<i>N</i>, <i>q</i><sub><i>t</i></sub>), &nbsp; <i>q</i><sub><i>t</i>+1</sub> = <i>X</i> / 2<i>N</i>"),
          vz("Rozptyl změny frekvence za jednu generaci", paste0("σ<sup>2</sup><sub>Δ<i>p,q</i></sub> = ", zl("<i>pq</i>", "2<i>N</i>"))),
          vz("Směrodatná odchylka změny frekvence jako míra genetického driftu", paste0("σ<sub>Δ<i>p,q</i></sub> = √(", zl("<i>pq</i>", "2<i>N</i>"), ")")),          
          vz("Pokles heterozygotnosti", paste0("<i>H</i><sub><i>t</i></sub> = <i>H</i><sub>0</sub>(1 − ", zl("1", "2<i>N</i>"), ")<sup><i>t</i></sup>, &nbsp; <i>H</i><sub>0</sub> = 2<i>q</i><sub>0</sub>(1 − <i>q</i><sub>0</sub>)")),
          vz("Pravděpodobnost fixace alely a", "<i>P</i><sub>fix</sub> = <i>q</i><sub>0</sub>")
        )
      ),
      
      mainPanel(
        width = 8,
        fluidRow(
          column(12, plotOutput("plot_drift_lines", height = "380px"))
        ),
        br(),
        fluidRow(
          column(12, plotOutput("plot_drift_het", height = "300px"))
        )
      )
    )
  ),
  
  # =================================================================
  # ZÁLOŽKA 6: INBREEDING (PŘÍBUZENSKÉ KŘÍŽENÍ F)
  # =================================================================
  tabPanel(
    "6. Inbreeding",
    sidebarLayout(
      sidebarPanel(
        width = 4,
        h3("Příbuzenské křížení (Inbreeding)"),
        helpText("Nenáhodné páření příbuzných jedinců zvyšuje homozygotnost v populaci měřené koeficientem inbreedingu F. Inbreeding ale vzniká i v panmiktické populaci konečné velikosti: čím je populace menší, tím rychleji F roste (viz graf akumulace)."),
        hr(),
        sliderInput("inbreed_p", "Frekvence alely A (p):", min = 0.01, max = 0.99, value = 0.50, step = 0.01),
        radioButtons("inbreed_zdroj", "Zdroj koeficientu inbreedingu F:",
                     choices = c(
                       "Zadat F ručně" = "manual",
                       "Spočítat F z velikosti populace N a počtu generací t" = "drift"
                     )),
        conditionalPanel(
          condition = "input.inbreed_zdroj == 'manual'",
          sliderInput("inbreed_F", "Koeficient inbreedingu (F):", min = 0.00, max = 1.00, value = 0.25, step = 0.01),
          helpText("Příklady F: F = 0.25 (plní sourozenci / rodič-potomek), F = 0.125 (polorodí sourozenci, strýc-neteř), F = 0.0625 (vlastní bratranci/sestřenice), F = 0 (panmixie).")
        ),
        conditionalPanel(
          condition = "input.inbreed_zdroj == 'drift'",
          sliderInput("inbreed_t", "Počet generací (t):", min = 1, max = 500, value = 20, step = 1),
          helpText("Při náhodném páření v konečné populaci mají i náhodně vybraní partneři společné předky.Na počátku roste F přibližně o 1/(2N) za generaci, ale skutečný přírůstek se s rostoucím F postupně zpomaluje. N chápejte zjednodušeně jako efektivní velikost populace Ne.")
        ),
        hr(),
        numericInput("inbreed_N", "Velikost populace (N):", value = 100, min = 10, max = 50000, step = 10),
        numericInput("inbreed_tmax", "Horizont grafu akumulace (generací):", value = 100, min = 10, max = 2000, step = 10),
        hr(),
        wellPanel(
          h4("Analýza inbrední zátěže:"),
          htmlOutput("inbreed_summary_html")
        ),
        vzorce_box(
          vz("Genotypové frekvence při inbreedingu F", "AA = <i>p</i><sup>2</sup> + <i>F</i><i>p</i><i>q</i><br>Aa = 2<i>p</i><i>q</i>(1 − <i>F</i>)<br>aa = <i>q</i><sup>2</sup> + <i>F</i><i>p</i><i>q</i>"),
          vz("Definice koeficientu F", paste0("<i>F</i> = 1 − ", zl("<i>H</i><sub>pozor.</sub>", "<i>H</i><sub>oček.</sub>"))),
          vz("Přebytek homozygotů / úbytek heterozygotů", "+<i>F</i><i>p</i><i>q</i> každý homozygot, &nbsp; −2<i>F</i><i>p</i><i>q</i> heterozygoti"),
          vz("Nárůst rizika recesivních homozygotů", paste0("RR = ", zl("<i>q</i><sup>2</sup> + <i>F</i><i>p</i><i>q</i>", "<i>q</i><sup>2</sup>"))),
          vz("Akumulace inbreedingu v konečné populaci (panmixie)",
             paste0("<i>F</i><sub><i>t</i></sub> = 1 − (1 − ", zl("1", "2<i>N</i>"), ")<sup><i>t</i></sup>")),
          vz("Přesný přírůstek inbreedingu za generaci", paste0("ΔF = (1 − F)", zl("1", "2N"))),
          vz("Počet generací do dosažení zadaného F",
             paste0("<i>t</i> = ", zl("ln(1 − <i>F</i>)", "ln(1 − 1/2<i>N</i>)"))),
          vz("Souvislost s driftem (úbytek heterozygotnosti)", "<i>H</i><sub><i>t</i></sub> = <i>H</i><sub>0</sub>(1 − <i>F</i><sub><i>t</i></sub>)")
        )
      ),
      
      mainPanel(
        width = 8,
        fluidRow(
          column(6, plotOutput("plot_inbreed_bar", height = "360px")),
          column(6, plotOutput("plot_inbreed_curve", height = "360px"))
        ),
        br(),
        fluidRow(
          column(12, plotOutput("plot_inbreed_time", height = "380px"))
        ),
        br(),
        fluidRow(
          column(12,
                 wellPanel(
                   h4("Srovnání genotypových četností: HWE (F = 0) vs. Inbreeding (F > 0)"),
                   tableOutput("tabulka_inbreed_comp")
                 )
          )
        )
      )
    )
  ),
  
  # =================================================================
  # ZÁLOŽKA 7: DRIFT VS. SELEKCE (KRITÉRIUM Ne*s)
  # =================================================================
  tabPanel(
    "7. Drift vs. Selekce",
    sidebarLayout(
      sidebarPanel(
        width = 4,
        h3("Drift vs. selekce (Ne · s)"),
        helpText("Wright-Fisherův model, ve kterém se v každé generaci nejprve uplatní selekce (výhodná alela a) a poté náhodné vzorkování gamet. O tom, zda převládne náhoda, nebo selekce, rozhoduje součin efektivní velikosti populace a selekčního koeficientu, Ne · s."),
        selectInput("ds_preset", "Výukový scénář (Preset):",
                    choices = c(
                      "Vlastní nastavení" = "custom",
                      "Drift dominuje (Ne·s = 0,5)" = "drift",
                      "Přechodová zóna (Ne·s = 3)" = "prechod",
                      "Selekce dominuje (Ne·s = 15)" = "selekce"
                    ), selected = "custom"),
        hr(),
        numericInput("ds_Ne", "Efektivní velikost populace Ne:", value = 100, min = 5, max = 5000, step = 10),
        numericInput("ds_s", "Selekční koeficient s (W_aa = 1 + s):", value = 0.03, min = 0, max = 0.5, step = 0.001),
        sliderInput("ds_h", "Dominance h (W_Aa = 1 + h·s):", min = 0, max = 1, value = 0.5, step = 0.05),
        sliderInput("ds_q0", "Počáteční frekvence výhodné alely a (q₀):", min = 0.01, max = 0.99, value = 0.10, step = 0.01),
        uiOutput("ds_ns_text"),
        hr(),
        sliderInput("ds_K", "Počet replikátních populací (linií):", min = 10, max = 200, value = 50, step = 10),
        numericInput("ds_G", "Počet generací (t):", value = 300, min = 10, max = 1000, step = 10),
        actionButton("btn_ds_sim", "Spustit novou simulaci", class = "btn-primary"),
        br(), br(),
        wellPanel(
          h4("Výsledek simulace:"),
          htmlOutput("ds_summary_html")
        ),
        vzorce_box(
          vz("Selekce (fitness: AA = 1, Aa = 1 + <i>hs</i>, aa = 1 + <i>s</i>)",
             paste0("<i>q</i><sub>sel</sub> = ", zl("<i>pq</i>(1 + <i>hs</i>) + <i>q</i><sup>2</sup>(1 + <i>s</i>)", "<i>p</i><sup>2</sup> + 2<i>pq</i>(1 + <i>hs</i>) + <i>q</i><sup>2</sup>(1 + <i>s</i>)"))),
          vz("Vzorkování gamet (drift) po selekci",
             "<i>X</i> ~ Bin(2<i>N</i><sub>e</sub>, <i>q</i><sub>sel</sub>), &nbsp; <i>q</i><sub><i>t</i>+1</sub> = <i>X</i> / 2<i>N</i><sub>e</sub>"),
          vz("Deterministická změna (selekce), aditivní případ <i>h</i> = ½",
             paste0("Δ<i>q</i><sub>sel</sub> ≈ ", zl("<i>s</i>", "2"), " <i>q</i>(1 − <i>q</i>)")),
          vz("Náhodná změna (drift)",
             paste0("Var(Δ<i>q</i><sub>drift</sub>) = ", zl("<i>q</i>(1 − <i>q</i>)", "2<i>N</i><sub>e</sub>"))),
          vz("Kritérium režimu",
             "<i>N</i><sub>e</sub> · <i>s</i> &lt; 1: drift dominuje; &nbsp; 1 ≤ <i>N</i><sub>e</sub> · <i>s</i> ≤ 10: přechodová zóna; &nbsp; <i>N</i><sub>e</sub> · <i>s</i> &gt; 10: selekce dominuje"),
          vz("Pravděpodobnost fixace (Kimura, aditivní selekce <i>h</i> = ½)",
             paste0("<i>P</i><sub>fix</sub> = ", zl("1 − e<sup>−2<i>N</i><sub>e</sub><i>s</i><i>q</i><sub>0</sub></sup>", "1 − e<sup>−2<i>N</i><sub>e</sub><i>s</i></sup>"), ", &nbsp; neutrálně <i>P</i><sub>fix</sub> = <i>q</i><sub>0</sub>")),
          helpText("Poznámka: přesná hodnota, na které selekce a drift soutěží, je v difuzní aproximaci přibližně 2·Ne·s (při aditivní selekci). Meze 1 a 10 jsou řádové a slouží k orientaci.")
        )
      ),
      
      mainPanel(
        width = 8,
        fluidRow(
          column(12, plotOutput("plot_ds_traj", height = "400px"))
        ),
        br(),
        fluidRow(
          column(12, plotOutput("plot_ds_regime", height = "340px"))
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
           x = "Mateřské gamety (vajíčka)", y = "Otcovské gamety (spermie)") +
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
           x = "Genotypy", y = "Relativní frekvence") +
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
      return(vypocti_hwe_2alely(n_AA, n_Aa, n_aa, perfect_hwe = TRUE))
    } else if (input$hwe_mode == "sample") {
      p <- input$hwe_p
      N <- input$hwe_N
      q <- 1 - p
      genotypes <- rbinom(N, size = 2, prob = q)
      n_AA <- sum(genotypes == 0)
      n_Aa <- sum(genotypes == 1)
      n_aa <- sum(genotypes == 2)
      return(vypocti_hwe_2alely(n_AA, n_Aa, n_aa, perfect_hwe = TRUE))
    } else {
      n_AA <- ifelse(is.na(input$manual_AA), 0, max(0, input$manual_AA))
      n_Aa <- ifelse(is.na(input$manual_Aa), 0, max(0, input$manual_Aa))
      n_aa <- ifelse(is.na(input$manual_aa), 0, max(0, input$manual_aa))
      return(vypocti_hwe_2alely(n_AA, n_Aa, n_aa, perfect_hwe = TRUE))
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
    sliderInput("abo_q", "Frekvence alely Iᴮ (q):", min = 0.01, max = max_q, value = val_q, step = 0.01)
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
      Genotyp = factor(c("IᴬIᴬ (p²)", "IᴮIᴮ (q²)", "ii (r²)", "IᴬIᴮ (2pq)", "Iᴬi (2pr)", "Iᴮi (2qr)"),
                       levels = c("IᴬIᴬ (p²)", "IᴮIᴮ (q²)", "ii (r²)", "IᴬIᴮ (2pq)", "Iᴬi (2pr)", "Iᴮi (2qr)")),
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
      scale_fill_manual(values = c("Skupina A (IᴬIᴬ + Iᴬi)" = "#4C72B0", 
                                   "Skupina B (IᴮIᴮ + Iᴮi)" = "#DD8452", 
                                   "Skupina AB (IᴬIᴮ)" = "#55A868", 
                                   "Skupina O (ii)" = "#C44E52")) +
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
      Kategorie = c("Alela Iᴬ (p)", "Alela Iᴮ (q)", "Alela i (r)",
                    "Genotyp IᴬIᴬ (p²)", "Genotyp IᴮIᴮ (q²)", "Genotyp ii (r²)",
                    "Genotyp IᴬIᴮ (2pq)", "Genotyp Iᴬi (2pr)", "Genotyp Iᴮi (2qr)",
                    "Fenotyp Skupina A (p² + 2pr)", "Fenotyp Skupina B (q² + 2qr)",
                    "Fenotyp Skupina AB (2pq)", "Fenotyp Skupina O (r²)"),
      Frekvence = c(round(res$p, 4), round(res$q, 4), round(res$r, 4),
                    round(unname(res$gen_freq["IᴬIᴬ"]), 4), round(unname(res$gen_freq["IᴮIᴮ"]), 4), round(unname(res$gen_freq["ii"]), 4),
                    round(unname(res$gen_freq["IᴬIᴮ"]), 4), round(unname(res$gen_freq["Iᴬi"]), 4), round(unname(res$gen_freq["Iᴮi"]), 4),
                    round(unname(res$phen_freq[1]), 4), round(unname(res$phen_freq[2]), 4),
                    round(unname(res$phen_freq[3]), 4), round(unname(res$phen_freq[4]), 4)),
      Ocekavany_Pocet = c("-", "-", "-",
                          round(unname(res$exp_gen_counts["IᴬIᴬ"]), 1), round(unname(res$exp_gen_counts["IᴮIᴮ"]), 1), round(unname(res$exp_gen_counts["ii"]), 1),
                          round(unname(res$exp_gen_counts["IᴬIᴮ"]), 1), round(unname(res$exp_gen_counts["Iᴬi"]), 1), round(unname(res$exp_gen_counts["Iᴮi"]), 1),
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

  # -----------------------------------------------------------------
  # LOGIKA PRO ZÁLOŽKU 5: GENETICKÝ DRIFT
  # -----------------------------------------------------------------
  drift_data <- reactive({
    input$btn_drift_sim
    isolate({
      simuluj_drift(
        N = input$drift_N,
        q0 = input$drift_q0,
        K = input$drift_K,
        generace = input$drift_G
      )
    })
  })
  
  output$plot_drift_lines <- renderPlot({
    res <- drift_data()
    ggplot(res$df_lines, aes(x = Generace, y = q, color = Populace)) +
      geom_line(size = 0.9, alpha = 0.8) +
      geom_hline(yintercept = c(0, 1), linetype = "dashed", color = "gray40") +
      labs(title = paste0("Trajektorie alely a u ", res$K, " linií (N = ", res$N, " diploidních jedinců)"),
           x = "Čas (generace t)", y = "Frekvence alely a (q)") +
      ylim(-0.02, 1.02) +
      theme_minimal(base_size = 13) +
      theme(
        legend.position = "none",
        plot.title = element_text(face = "bold"),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      )
  })
  
  output$plot_drift_het <- renderPlot({
    res <- drift_data()
    ggplot(res$df_het, aes(x = Generace)) +
      geom_line(aes(y = H_obs, color = "Pozorovaná průměrná heterozygotnost (2pq)"), size = 1.2) +
      geom_line(aes(y = H_teor, color = "Teoretická očakávaná H_t = H₀(1 - 1/2N)^t"), size = 1.2, linetype = "dashed") +
      scale_color_manual(values = c("Pozorovaná průměrná heterozygotnost (2pq)" = "#D55E00", 
                                    "Teoretická očakávaná H_t = H₀(1 - 1/2N)^t" = "#0072B2")) +
      labs(title = "Rychlost ztráty heterozygotnosti v důsledku genetického driftu",
           x = "Čas (generace t)", y = "Průměrná heterozygotnost (H)", color = "") +
      ylim(0, max(0.6, max(res$df_het$H_obs) * 1.1)) +
      theme_minimal(base_size = 13) +
      theme(
        legend.position = "top",
        plot.title = element_text(face = "bold"),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      )
  })
  
  output$drift_summary_html <- renderUI({
    res <- drift_data()
    mean_q_final <- round(mean(res$last_q), 4)
    
    HTML(paste0(
      "<b>Fixované linie (q = 1):</b> <span style='color: #D55E00; font-weight: bold;'>", res$fixed_a, "</span> (", round(res$fixed_a / res$K * 100, 1), "%)<br/>",
      "<b>Ztracené linie (q = 0):</b> <span style='color: #0072B2; font-weight: bold;'>", res$lost_a, "</span> (", round(res$lost_a / res$K * 100, 1), "%)<br/>",
      "<b>Polymorfní linie (0 < q < 1):</b> ", res$poly, " (", round(res$poly / res$K * 100, 1), "%)<br/>",
      "<b>Průměrná q v generaci ", res$generace, ":</b> ", mean_q_final, " (počáteční q₀ = ", res$q0, ")"
    ))
  })

  # -----------------------------------------------------------------
  # LOGIKA PRO ZÁLOŽKU 6: INBREEDING
  # -----------------------------------------------------------------
  # Efektivní F: buď zadané ručně, nebo vypočtené z velikosti populace a počtu generací
  inbreed_F_eff <- reactive({
    req(is.finite(input$inbreed_N), input$inbreed_N >= 2)
    if (identical(input$inbreed_zdroj, "drift")) {
      req(is.finite(input$inbreed_t))
      round(inbreeding_Ft(input$inbreed_N, input$inbreed_t), 3)
    } else {
      req(is.finite(input$inbreed_F))
      input$inbreed_F
    }
  })
  
  inbreed_data <- reactive({
    vypocti_inbreeding(p = input$inbreed_p, F_coeff = inbreed_F_eff(), N = input$inbreed_N)
  })
  
  output$plot_inbreed_bar <- renderPlot({
    res <- inbreed_data()
    df_bar <- data.frame(
      Genotyp = factor(rep(c("AA", "Aa", "aa"), 2), levels = c("AA", "Aa", "aa")),
      Stav = factor(rep(c("Panmixie (F = 0)", paste0("Inbreeding (F = ", res$F, ")")), each = 3),
                    levels = c("Panmixie (F = 0)", paste0("Inbreeding (F = ", res$F, ")"))),
      Frekvence = c(unname(res$freq_0), unname(res$freq_F)),
      Pocet = c(unname(res$counts_0), unname(res$counts_F))
    )
    
    barvy_stav <- setNames(c("#4C72B0", "#C44E52"), c("Panmixie (F = 0)", paste0("Inbreeding (F = ", res$F, ")")))
    
    ggplot(df_bar, aes(x = Genotyp, y = Frekvence, fill = Stav)) +
      geom_bar(stat = "identity", position = "dodge", width = 0.7) +
      geom_text(aes(label = paste0(round(Frekvence * 100, 1), "%")),
                position = position_dodge(width = 0.7), vjust = -0.3, fontface = "bold", size = 3.8) +
      scale_fill_manual(values = barvy_stav) +
      labs(title = "Srovnání genotypových frekvencí: Panmixie vs. Inbreeding",
           x = "Genotyp", y = "Relativní frekvence", fill = "") +
      ylim(0, max(df_bar$Frekvence) * 1.2) +
      theme_minimal(base_size = 13) +
      theme(
        legend.position = "top",
        plot.title = element_text(face = "bold"),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      )
  })
  
  output$plot_inbreed_curve <- renderPlot({
    res <- inbreed_data()
    p_val <- res$p
    q_val <- res$q
    F_seq <- seq(0, 1, length.out = 100)
    
    f_AA <- p_val^2 + F_seq * p_val * q_val
    f_Aa <- 2 * p_val * q_val * (1 - F_seq)
    f_aa <- q_val^2 + F_seq * p_val * q_val
    
    lev_AA <- "AA (p² + Fpq)"
    lev_Aa <- "Aa (2pq(1-F))"
    lev_aa <- "aa (q² + Fpq)"
    
    # Pořadí úrovní = pořadí vykreslování: AA se kreslí poslední (navrch),
    # takže je vidět i tehdy, když se kryje s křivkou aa (p = q = 0,5).
    df_curve <- data.frame(
      F_val = rep(F_seq, 3),
      Frekvence = c(f_aa, f_Aa, f_AA),
      Genotyp = factor(rep(c(lev_aa, lev_Aa, lev_AA), each = 100),
                       levels = c(lev_aa, lev_Aa, lev_AA))
    )
    
    podtitul <- if (abs(p_val - q_val) < 1e-9) {
      "Při p = q se křivky AA (čárkovaně) a aa překrývají, protože p² + Fpq = q² + Fpq."
    } else {
      "Křivka AA je kreslena čárkovaně."
    }
    
    ggplot(df_curve, aes(x = F_val, y = Frekvence, color = Genotyp, linetype = Genotyp)) +
      geom_line(linewidth = 1.2) +
      geom_vline(xintercept = res$F, linetype = "dashed", color = "gray30") +
      scale_color_manual(values = c("AA (p² + Fpq)" = "#D55E00", "Aa (2pq(1-F))" = "#0072B2", "aa (q² + Fpq)" = "#009E73"),
                         breaks = c(lev_AA, lev_Aa, lev_aa)) +
      scale_linetype_manual(values = c("AA (p² + Fpq)" = "dashed", "Aa (2pq(1-F))" = "solid", "aa (q² + Fpq)" = "solid"),
                            breaks = c(lev_AA, lev_Aa, lev_aa)) +
      labs(title = paste0("Změna genotypových frekvencí v závislosti na F (p = ", p_val, ")"),
           subtitle = podtitul,
           x = "Koeficient inbreedingu (F)", y = "Frekvence genotypu") +
      ylim(0, 1) +
      theme_minimal(base_size = 13) +
      theme(
        legend.position = "top",
        plot.title = element_text(face = "bold"),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      )
  })
  
  output$inbreed_summary_html <- renderUI({
    res <- inbreed_data()
    excess_homo <- round(res$F * res$p * res$q, 4)
    deficit_het <- round(2 * res$F * res$p * res$q, 4)
    risk_txt <- if (!is.na(res$risk_ratio)) paste0(round(res$risk_ratio, 2), "x") else "N/A"
    
    if (identical(input$inbreed_zdroj, "drift")) {
      radek_zdroj <- paste0("<b>Inbreeding po ", input$inbreed_t, " generacích (N = ", res$N, "):</b> ",
                            "<span style='font-weight: bold;'>F = ", res$F, "</span><br/>")
    } else {
      t_cil <- inbreeding_t_do_F(res$N, res$F)
      t_txt <- if (is.finite(t_cil)) paste0(ceiling(t_cil), " generací") else "nikdy (F = 1 nelze dosáhnout)"
      radek_zdroj <- paste0("<b>Panmixie o velikosti N = ", res$N, " dosáhne F = ", res$F, " za:</b> ",
                            "<span style='font-weight: bold;'>", t_txt, "</span><br/>")
    }
    
    HTML(paste0(
      "<b>Přebytek homozygotů (+Fpq):</b> <span style='color: #D55E00; font-weight: bold;'>+", excess_homo, "</span> na každý homozygotný genotyp<br/>",
      "<b>Úbytek heterozygotů (-2Fpq):</b> <span style='color: #C44E52; font-weight: bold;'>-", deficit_het, "</span><br/>",
      "<b>Nárůst rizika recesivních homozygotů (aa):</b> <span style='color: #721c24; font-weight: bold;'>", risk_txt, "</span> vyšší výskyt oproti panmixii<br/>",
      "<small>U vzácných patogenních alel (malé q) inbreeding dramaticky zvyšuje incidenci onemocnění!</small>",
      "<hr style='margin: 6px 0;'/>",
      "<b>Přírůstek inbreedingu za generaci (ΔF ≈ 1/2N):</b> ", round(1 / (2 * res$N), 4), "<br/>",
      radek_zdroj
    ))
  })
  
  output$plot_inbreed_time <- renderPlot({
    req(is.finite(input$inbreed_tmax), input$inbreed_tmax >= 1)
    res <- inbreed_data()
    N_user <- res$N
    drift_mode <- identical(input$inbreed_zdroj, "drift")
    tmax <- max(input$inbreed_tmax, if (drift_mode) input$inbreed_t else 0)
    t_seq <- 0:tmax
    
    N_vals <- sort(unique(c(10, 50, 500, 5000, N_user)))
    lab_N <- function(n) ifelse(n == N_user, paste0("N = ", n, " (vybrané)"), paste0("N = ", n))
    
    df_t <- do.call(rbind, lapply(N_vals, function(n) {
      data.frame(t = t_seq, F_t = inbreeding_Ft(n, t_seq), N = n)
    }))
    df_t$N_lab <- factor(lab_N(df_t$N), levels = lab_N(N_vals))
    df_t$vybrane <- df_t$N == N_user
    
    g <- ggplot(df_t, aes(x = t, y = F_t, color = N_lab, linewidth = vybrane)) +
      geom_hline(yintercept = res$F, linetype = "dashed", color = "gray30") +
      annotate("text", x = tmax, y = res$F, label = paste0("F = ", res$F),
               hjust = 1, vjust = -0.6, size = 3.8, color = "gray30") +
      geom_line() +
      scale_linewidth_manual(values = c("FALSE" = 0.8, "TRUE" = 1.8), guide = "none") +
      scale_color_viridis_d(end = 0.8) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(title = "Akumulace inbreedingu v konečné populaci",
           subtitle = "F(t) = 1 − (1 − 1/2N)ᵗ: v malých populacích roste inbreeding mnohem rychleji",
           x = "Generace (t)", y = "Koeficient inbreedingu F(t)", color = "")
    
    if (drift_mode) {
      g <- g + geom_vline(xintercept = input$inbreed_t, linetype = "dotted", color = "gray30") +
        geom_point(data = data.frame(t = input$inbreed_t, F_t = res$F),
                   aes(x = t, y = F_t), inherit.aes = FALSE, size = 3.5, color = "black")
    } else {
      t_cil <- inbreeding_t_do_F(N_user, res$F)
      if (is.finite(t_cil) && t_cil <= tmax) {
        g <- g + geom_vline(xintercept = t_cil, linetype = "dotted", color = "gray30") +
          geom_point(data = data.frame(t = t_cil, F_t = res$F),
                     aes(x = t, y = F_t), inherit.aes = FALSE, size = 3.5, color = "black")
      }
    }
    
    g +
      theme_minimal(base_size = 13) +
      theme(
        legend.position = "right",
        plot.title = element_text(face = "bold"),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      )
  })
  
  output$tabulka_inbreed_comp <- renderTable({
    res <- inbreed_data()
    fmt <- function(x) formatC(x, format = "f", digits = 4)
    gen <- c("AA", "Aa", "aa")
    data.frame(
      Genotyp = c(gen, "Celkem"),
      Frekvence_Panmixie = c(fmt(unname(res$freq_0[gen])), fmt(1)),
      Pocet_Panmixie = c(as.integer(res$counts_0[gen]), as.integer(res$N)),
      Frekvence_Inbreeding = c(fmt(unname(res$freq_F[gen])), fmt(1)),
      Pocet_Inbreeding = c(as.integer(res$counts_F[gen]), as.integer(res$N)),
      Rozdil_Pocet = c(as.integer(res$counts_F[gen] - res$counts_0[gen]), 0L)
    )
  }, digits = 0)

  # -----------------------------------------------------------------
  # LOGIKA PRO ZÁLOŽKU 7: DRIFT VS. SELEKCE
  # -----------------------------------------------------------------
  observeEvent(input$ds_preset, {
    if (input$ds_preset == "drift") {
      updateNumericInput(session, "ds_Ne", value = 100)
      updateNumericInput(session, "ds_s", value = 0.005)
      updateNumericInput(session, "ds_G", value = 400)
    } else if (input$ds_preset == "prechod") {
      updateNumericInput(session, "ds_Ne", value = 100)
      updateNumericInput(session, "ds_s", value = 0.03)
      updateNumericInput(session, "ds_G", value = 600)
    } else if (input$ds_preset == "selekce") {
      updateNumericInput(session, "ds_Ne", value = 100)
      updateNumericInput(session, "ds_s", value = 0.15)
      updateNumericInput(session, "ds_G", value = 250)
    }
    if (input$ds_preset != "custom") {
      updateSliderInput(session, "ds_h", value = 0.5)
      updateSliderInput(session, "ds_q0", value = 0.10)
    }
  }, ignoreInit = TRUE)
  
  output$ds_ns_text <- renderUI({
    req(is.finite(input$ds_Ne), is.finite(input$ds_s))
    x <- input$ds_Ne * input$ds_s
    r <- ds_rezim(x)
    HTML(paste0(
      "<b>Ne · s =</b> <span style='font-weight: bold;'>", round(x, 3), "</span><br/>",
      "<span style='color: ", r$barva, "; font-weight: bold;'>", r$text, "</span>"
    ))
  })
  
  ds_data <- reactive({
    input$btn_ds_sim   # tlačítko = nový náhodný běh
    req(is.finite(input$ds_Ne), input$ds_Ne >= 5,
        is.finite(input$ds_s), input$ds_s >= 0,
        is.finite(input$ds_G), input$ds_G >= 10)
    simuluj_drift_selekce(
      Ne = input$ds_Ne, s = input$ds_s, h = input$ds_h,
      q0 = input$ds_q0, K = input$ds_K, generace = input$ds_G
    )
  })
  
  output$plot_ds_traj <- renderPlot({
    res <- ds_data()
    G <- res$generace
    K <- res$K
    r <- ds_rezim(res$x)
    
    df_lines <- data.frame(
      Generace = rep(0:G, times = K),
      Linie = rep(seq_len(K), each = G + 1),
      q = as.vector(res$Q)
    )
    lab_det <- "Deterministicky (jen selekce, bez driftu)"
    lab_mean <- "Průměr simulovaných linií"
    lab_neu <- "Neutrální očekávání (E[q] = q₀)"
    df_ref <- data.frame(
      Generace = rep(0:G, 3),
      q = c(res$q_det, res$q_mean, rep(res$q0, G + 1)),
      Typ = factor(rep(c(lab_det, lab_mean, lab_neu), each = G + 1),
                   levels = c(lab_det, lab_mean, lab_neu))
    )
    
    ggplot() +
      geom_line(data = df_lines, aes(x = Generace, y = q, group = Linie),
                color = "#7f8fa6", alpha = 0.45, linewidth = 0.5) +
      geom_hline(yintercept = c(0, 1), linetype = "dashed", color = "gray40") +
      geom_line(data = df_ref, aes(x = Generace, y = q, color = Typ, linetype = Typ), linewidth = 1.3) +
      scale_color_manual(values = setNames(c("black", "#D55E00", "gray45"), c(lab_det, lab_mean, lab_neu))) +
      scale_linetype_manual(values = setNames(c("solid", "solid", "dotted"), c(lab_det, lab_mean, lab_neu))) +
      labs(title = paste0("Drift vs. selekce: Ne = ", res$Ne, ", s = ", res$s,
                          "  (Ne·s = ", round(res$x, 2), ")"),
           subtitle = r$text,
           x = "Čas (generace t)", y = "Frekvence výhodné alely a (q)", color = "", linetype = "") +
      ylim(-0.02, 1.02) +
      theme_minimal(base_size = 13) +
      theme(
        legend.position = "top",
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(color = r$barva, face = "bold"),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      ) +
      guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))
  })
  
  output$plot_ds_regime <- renderPlot({
    res <- ds_data()
    x_seq <- 10^seq(-2, 3, length.out = 300)
    df_c <- data.frame(x = x_seq, P = ds_pfix(x_seq, res$q0))
    
    p_sim <- res$fixed / res$K
    p_teor <- ds_pfix(res$x, res$q0)
    x_ok <- res$x >= 0.01 && res$x <= 1000
    
    g <- ggplot(df_c, aes(x = x, y = P)) +
      annotate("rect", xmin = 0.01, xmax = 1, ymin = 0, ymax = 1, fill = "#0072B2", alpha = 0.10) +
      annotate("rect", xmin = 1, xmax = 10, ymin = 0, ymax = 1, fill = "#E69F00", alpha = 0.14) +
      annotate("rect", xmin = 10, xmax = 1000, ymin = 0, ymax = 1, fill = "#009E73", alpha = 0.10) +
      annotate("text", x = 0.1, y = 1.04, label = "DRIFT", fontface = "bold", color = "#0072B2") +
      annotate("text", x = sqrt(10), y = 1.04, label = "PŘECHOD", fontface = "bold", color = "#B07800") +
      annotate("text", x = 100, y = 1.04, label = "SELEKCE", fontface = "bold", color = "#009E73") +
      geom_hline(yintercept = res$q0, linetype = "dashed", color = "gray40") +
      annotate("text", x = 0.011, y = res$q0, label = paste0("neutrálně P = q₀ = ", res$q0),
               hjust = 0, vjust = -0.6, size = 3.6, color = "gray30") +
      geom_line(linewidth = 1.3, color = "#1a2b4c")
    
    if (x_ok) {
      g <- g +
        geom_vline(xintercept = res$x, linetype = "dotted", color = "gray30") +
        geom_point(data = data.frame(x = res$x, P = p_teor), size = 4, color = "#1a2b4c") +
        geom_point(data = data.frame(x = res$x, P = p_sim), size = 4, shape = 18, color = "#D55E00")
    }
    
    g +
      scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100, 1000),
                    labels = c("0,01", "0,1", "1", "10", "100", "1000")) +
      coord_cartesian(ylim = c(0, 1.07)) +
      labs(title = "Pravděpodobnost fixace výhodné alely v závislosti na Ne·s",
           subtitle = if (abs(res$h - 0.5) > 1e-9) "Teoretická křivka platí pro aditivní selekci (h = 0,5)." else NULL,
           caption = paste0("● teorie (Kimura, h = ½)     ◆ simulace: podíl fixovaných linií v generaci ", res$generace,
                            if (!x_ok) "     (aktuální Ne·s je mimo zobrazený rozsah)" else ""),
           x = "Ne · s  (logaritmická osa)", y = "Pravděpodobnost fixace") +
      theme_minimal(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0, size = 10),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      )
  })
  
  output$ds_summary_html <- renderUI({
    res <- ds_data()
    r <- ds_rezim(res$x)
    K <- res$K
    p_teor <- ds_pfix(res$x, res$q0)
    teor_txt <- if (abs(res$h - 0.5) < 1e-9) {
      paste0("<b>Teoretická pravděpodobnost fixace (Kimura):</b> ", round(p_teor, 3), "<br/>")
    } else {
      "<small>Kimurova teorie platí pro h = 0,5.</small><br/>"
    }
    
    HTML(paste0(
      "<b>Režim:</b> <span style='color: ", r$barva, "; font-weight: bold;'>", r$text, "</span><br/>",
      "<b>Fixované linie (q = 1):</b> <span style='color: #D55E00; font-weight: bold;'>", res$fixed, "</span> (", round(res$fixed / K * 100, 1), "%)<br/>",
      "<b>Ztracené linie (q = 0):</b> <span style='color: #0072B2; font-weight: bold;'>", res$lost, "</span> (", round(res$lost / K * 100, 1), "%)<br/>",
      "<b>Polymorfní linie (0 < q < 1):</b> ", res$poly, " (", round(res$poly / K * 100, 1), "%)<br/>",
      teor_txt,
      "<b>Neutrální očekávání (drift):</b> ", res$q0, "<br/>",
      "<b>Deterministické q v generaci ", res$generace, ":</b> ", round(res$q_det[res$generace + 1], 3), "<br/>",
      "<b>Průměrné q simulovaných linií:</b> ", round(res$q_mean[res$generace + 1], 3),
      if (res$poly > 0) paste0("<br/><small>Pozn.: ", res$poly, " linií ještě segreguje, simulovaný podíl fixace je proto v čase t dolním odhadem. Prodlužte počet generací.</small>") else ""
    ))
  })
}

# -------------------------------------------------------------------
# 4. SPUŠTĚNÍ APLIKACE
# -------------------------------------------------------------------
shinyApp(ui = ui, server = server)
