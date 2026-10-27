#### main function for regressions ####

# run_regressions: Automatizza l'esecuzione dei modelli panel (Fixed Effects) in base ai parametri
# Argomenti:
#   - data: Il dataset contenente tutte le variabili.
#   - emission_var: Vettore con gli indici o nomi delle colonne della variabile dipendente (CO2).
#   - ine_var: Vettore con gli indici o nomi delle colonne della variabile indipendente principale (disuguaglianza).
#   - formulas: Finestra di formule di regressione (es. dl1, dl2, dl5) da stimare.
#   - weight: Booleano (TRUE/FALSE). Se FALSE usa il pacchetto 'plm' (non pesato, output stargazer). 
#             Se TRUE usa il pacchetto 'fixest' (pesato su popolazione, output etable).
#   - start_year: L'anno di inizio da cui filtrare il campione (es: 1980).
#   - min_obs: Numero minimo di osservazioni valide per paese (di default 1).
#   - oecd_only
#   - output_type: Formato di output per le tabelle ("text" di default, oppure "latex")

run_regressions <- function(data, emission_var, ine_var, formulas, weight = FALSE, start_year = 1980, min_obs = 1, oecd_only = F, output_type = "text") {
  
  # Assicurati che "fixest" sia installato se si utilizza l'approccio pesato
  if (weight) {
    if (!requireNamespace("fixest", quietly = TRUE)) stop("Package 'fixest' is required when weight=TRUE")
  }
  
  # Inizializzo vettori vuoti per salvare statistiche aggregate estratte durante i loop
  all_thresholds <- c()
  all_c_a_tps <- c()
  all_beta_0s <- c()
  all_beta_ints <- c()
  all_tables <- list()
  
  f_idx <- 1
  # Loop 1: Itera su tutte le formule fornite dall'utente
  for (f in formulas) { 
    formula <- f
    
    # Liste vuote per salvare ciascun modello temporaneo per poi farne una tabella unica
    fee <- list()
    see <- list()
    i <- 1 
    
    # Statistiche specifiche del modello corrente (formula f), una lista per ogni formuala
    xy <- c() #main variables used
    thresholds <- c()
    p_used <- c() #numero di paesi usati
    beta_0s <- c()
    beta_ints <- c()
    c_a_tps <- c()
    
    # Loop 2: Itera su ogni variante del termine per la CO2
    for (co2_col in emission_var) { 
      
      # Loop 3: Itera su ogni variante del termine per l'ineguaglianza
      for (ineq_col in ine_var) {
        
        # Step 1: Preparazione dataset temporaneo
        # Le colonne target vengono unificate nei nomi "co2" e "inequality" così che
        # le formule in R capiscano automaticamente quali variabili assegnare
        panel_data <- data
        panel_data$co2 <- data[,co2_col]   
        panel_data$inequality <- data[,ineq_col]
        
        # Rimozione dei missing values dalle principali x e y per pulizia primaria
        panel_data <- subset(panel_data, (!is.na(panel_data$co2)) & (!is.na(panel_data$inequality)))    
        
        # Step 2: Filtro il campione
        # Conto i record rimanenti: mantengo solo i paesi che raggiungono o superano min_obs.
        tutti <- plyr::count(panel_data$Country.Code)$x
        vf <- plyr::count(panel_data$Country.Code)$freq >= min_obs
        paesi <- tutti[vf]
        panel_data <- panel_data %>% filter(Country.Code %in% paesi)  
        
        # Filtro dall'anno di inizio
        panel_data <- panel_data %>% filter(Time >= start_year)
        
        # Filtro eventualmente OECD
        if(oecd_only){panel_data <- panel_data %>% filter(oecd == "OECD")}
        
        # Creazione della struttura panel dataset formale
        panel <- pdata.frame(panel_data, index=c("Country.Code","Time"))
        
        # Creo una legenda pulita (LaTeX-friendly) combinando i veri nomi delle variabili
        legenda <- colnames(data[ineq_col])
        
        # Inizializzo le metriche di test
        threshold <- NA; beta_0 <- NA; beta_int <- NA; c_a_tp <- NA; cn <- NA;
        
        if (!weight) {
          # ========== MODELLO SENZA PESI (plm) ==========
          # Modello fixed effetti "twoways" (paese e temporalità)
          fe <- plm(formula, model="within", effect = "twoways", data=panel)
          
          # Compute Small Sample Correction tipico di Stata: N / (N - 1)
          n_countries <- length(unique(panel_data$Country.Code))
          small_sample_correction <- n_countries/(n_countries - 1)
          
          # Estrazione degli errori standard robusti clusterizzati per gruppo (paese)
          coef_test <- coeftest(fe, small_sample_correction * vcovHC(fe, type = "HC2", cluster = "group"))
          se <- sqrt(diag(small_sample_correction * vcovHC(fe, type="HC2", cluster = "group")))
          
          # Calcolo numero di paesi intercettati effettivamente dal dataset residuale
          c_used <- as.data.frame(fe$residual)
          c_used <- tibble::rownames_to_column(c_used, "cc")
          c_used$cc <- substr(c_used$cc, 1, 3) 
          cn <- length(unique(c_used$cc))
          
          # Utilizzo logiche dinamiche per rintracciare i coefficienti dell'inequality main effect
          # e dell'interazione (inequality : variabile) che servono per calcolare il turning point
          c_names <- rownames(coef_test)
          int_name <- c_names[grepl(":", c_names) & grepl("inequality", c_names)][1]
          main_name <- c_names[grepl("inequality", c_names) & !grepl(":", c_names)][1]
          if (is.na(main_name)) main_name <- c_names[1]
          
          if (!is.na(int_name) && !is.na(main_name)) {
            beta_0 <- as.numeric(coef_test[main_name, 1])
            beta_int <- as.numeric(coef_test[int_name, 1])
            
            # Calcolo Turning Point: exp(-beta_linear / beta_interaz) - 1
            threshold <- round(exp(-(beta_0)/beta_int)-1, 2)
          }
          
          if (!is.na(threshold)) {
            # Conto quanti dataset superano il TP alla fine del periodo osservato (2018 default max time)
            c_a_tp <- ((panel_data %>% filter(Time==max(Time)) %>% filter(Country.Code %in% c_used$cc) %>%  filter(serv_managr>threshold) %>% tally()) / cn)[1,1]
          }
        } else {
          # ========== MODELLO PESATO (fixest) ==========
          # Se weight = T, si usa il pacchetto performante "feols" aggiungendo i fissi per stringa
          f_str <- paste(deparse(formula), collapse = " ")
          if (!grepl("\\|", f_str)) {
            f_str <- paste(f_str, "| Country.Code + Time")
            formula <- as.formula(f_str)
          }
          # Pesato per popolazione (pop) -> "weights = ~ pop"
          fe <- fixest::feols(formula, data = panel, weights = ~ pop, cluster = ~ Country.Code)
          # Cluster standard error incorporato automaticamente
          se <- fixest::se(fe)
          
          c_test <- fe$coeftable
          used_rows <- fixest::obs(fe)
          countries_used <- unique(as.character(panel$Country.Code[used_rows]))
          cn <- length(countries_used)
          
          # Stessa procedura per individuare i coefficienti e derivare il turning point dinamicamente
          c_names <- rownames(c_test)
          int_name <- c_names[grepl(":", c_names) & grepl("inequality", c_names)][1]
          main_name <- c_names[grepl("inequality", c_names) & !grepl(":", c_names)][1]
          if (is.na(main_name)) main_name <- c_names[1]
          
          if (!is.na(int_name) && !is.na(main_name)) {
            beta_0 <- as.numeric(c_test[main_name, 1]) # Column 1 corresponds to "Estimate"
            beta_int <- as.numeric(c_test[int_name, 1])
            threshold <- round(exp(-(beta_0)/beta_int)-1, 2)
          }
          
          if (!is.na(threshold)) {
            c_a_tp <- ((panel_data %>% filter(Time==max(Time)) %>% filter(serv_managr>threshold) %>% tally()) / cn)[1,1]
          }
        }
        
        # Immagazzinamento del loop singolo
        fee[[i]] <- fe
        see[[i]] <- se
        i <- i + 1
        
        xy <- c(xy, legenda)
        thresholds <- c(thresholds, threshold)
        p_used <- c(p_used, cn)
        beta_0s <- c(beta_0s, beta_0)
        beta_ints <- c(beta_ints, beta_int)
        c_a_tps <- c(c_a_tps, c_a_tp)
      }
    }
    
    var_labels <- c(
      "inequality"         = "Inequality (ln)",    
      "l_serv_managr"      = "Tertiarization (ln)",
      "l_gdp_pc"           = "GDP per capita (ln)",
      "l_pop"              = "Population (ln)",
      "l_urban_pop"        = "Urban population (ln)",
      "l_services"         = "Services, value-added (ln)",
      "l_prim_en"          = "Primary enrollment (ln)",
      "l_over_65"          = "Population ages 65 and above (ln)",
      "l_final_hous_cons"  = "Final households consumption per capita (ln)",
      "l_trade_gdp"        = "Trade (ln)",
      "l_fdi_i"            = "FDI inflows (ln)",
      "l_manufacturing"    = "Manufacturing, value-added (ln)",
      "l_industry"         = "Industry, value-added (ln)",
      "str_index"          = "Environmental Policy Stringency Index",
      "l_empl_serv_managr" = "Tertiarization (employment) (ln)",
      "l_serv_empl"        = "Employment in services (ln)",
      "l_serv_indagr"      = "Tertiarization (with industry)",
      "l_serv_managrserv"  = "(Services/ All sectors) (with manufacturing)",
      "l_serv_indagrserv"  = "(Services/ All sectors) (with industry)",
      "inequality:l_serv_managr" = "Inequality*Tertiarization",
      "inequality:l_empl_serv_managr" = "Inequality*Tertiarization (employment)",
      "inequality:l_serv_indagr" = "Inequality*Tertiarization (with industry)",
      "inequality:l_serv_managrserv" = "Inequality*(Services/ All sectors) (with manufacturing)",
      "inequality:l_serv_indagrserv" = "Inequality*(Services/ All sectors) (with industry)"
    )
    
    if (output_type == "latex") {
      coef_names <- unique(unlist(lapply(fee, function(m) names(coef(m)))))
      coef_labels <- sapply(coef_names, function(x) {
        if (x %in% names(var_labels)) return(var_labels[[x]])
        return(x)
      }, USE.NAMES = FALSE)
    } else {
      coef_labels <- NULL
    }
    
    # Label lookups
    ine_labels <- c(
      "l_gini_wid"      = "Gini pre-tax",
      "l_palma"          = "Palma pre-tax",
      "l_r20"            = "20:20 ratio pre-tax",
      "l_p80p100"        = "Top 20\\% pre-tax",
      "l_p90p100"        = "Top 10\\% pre-tax",
      "l_gini_wid_disp"  = "Gini post-tax",
      "l_palma_d"        = "Palma post-tax",
      "l_r20_d"          = "20:20 ratio post-tax",
      "l_p80p100d"       = "Top 20\\% post-tax",
      "l_p90p100d"       = "Top 10\\% post-tax"
    )
    xy_labels <- ifelse(xy %in% names(ine_labels), ine_labels[xy], sapply(xy, function(x) xtable::sanitize(x, type = "latex")))
    
    # Before the stargazer call, build the caption dynamically
    dep_var_name <- colnames(data[co2_col])
    dep_var_label <- if (dep_var_name == "l_co2_pc_wb") {
      "CO$_{2}$ emissions from production (ln)"
    } else if (dep_var_name == "l_eora_cons") {
      "CO$_{2}$ emissions from consumption (ln)"
    } else {
      dep_var_name
    }
    
    # Step 3: Assemblaggio output a schermo/oggetto list.
    if (!weight) {
      tb <- capture.output(
        stargazer::stargazer(fee, se = see,
                             column.labels = xy_labels,
                             dep.var.labels.include = F,
                             dep.var.caption = paste0("\\textit{Dependent variable: ", dep_var_label, "}"),
                             covariate.labels = coef_labels,
                             add.lines=list(
                               c("Country fixed effects", rep("\\checkmark", length(thresholds))),
                               c("Year dummies", rep("\\checkmark", length(thresholds))),
                               c("Turning point", thresholds),
                               c("Countries", p_used),
                               c("\\% of obs. $>$ turning point", round(c_a_tps,2)*100)
                             ),
                             no.space = TRUE,
                             type = output_type)
      )
      if (output_type == "latex") {
        # Add \\[-1.8ex] before column labels
        col_labels_idx <- grep("\\\\cline\\{", tb)
        if (length(col_labels_idx) > 0) {
          for (idx in col_labels_idx) {
            if (idx < length(tb)) {
              tb[idx + 1] <- paste0("\\\\[-1.8ex] ", tb[idx + 1])
            }
          }
        }
        # Add adjustbox
        caption_idx <- grep("\\\\caption\\{", tb)
        if (length(caption_idx) > 0) {
          tb <- append(tb, "\\begin{adjustbox}{width=0.97\\textwidth,center}", after = caption_idx)
        }
        endtab_idx <- grep("\\\\end\\{tabular\\}", tb)
        if (length(endtab_idx) > 0) {
          tb <- append(tb, "\\end{adjustbox}", after = endtab_idx)
        }
      }
      cat(tb, sep = "\n")
      
      all_tables[[f_idx]] <- tb
    } else {
      names(fee) <- xy_labels
      is_tex <- if (output_type == "latex") TRUE else FALSE
      tb <- fixest::etable(fee,
                           extralines = list("Country fixed effects" = rep("\\checkmark", length(thresholds)),
                                             "Year dummies" = rep("\\checkmark", length(thresholds)),
                                             "Turning point" = thresholds,
                                             "Countries" = p_used,
                                             "\\% of obs. $>$ turning point" = round(c_a_tps,2)*100),
                           tex = is_tex)
      all_tables[[f_idx]] <- tb
      print(tb)
    }
    
    f_idx <- f_idx + 1
    all_thresholds <- c(all_thresholds,thresholds)
    all_beta_0s <- c(all_beta_0s,beta_0s)
    all_beta_ints <- c(all_beta_ints,beta_ints)
    all_c_a_tps <- c(all_c_a_tps,c_a_tps)
  }
  
  # Restituisci la lista contenente metriche e reports tabellari
  return(list(
    all_thresholds = all_thresholds,
    all_c_a_tps = all_c_a_tps,
    tables = all_tables
  ))
}

#### function for regressions to use to create Figures 3 and 4 ####
run_regressions2 <- function(data, emission_var, ine_var, formulas, weight = FALSE, start_year = 1980, min_obs = 1, oecd_only = F, output_type = "text") {
  
  # Assicurati che "fixest" sia installato se si utilizza l'approccio pesato
  if (weight) {
    if (!requireNamespace("fixest", quietly = TRUE)) stop("Package 'fixest' is required when weight=TRUE")
  }
  
  # Inizializzo vettori vuoti per salvare statistiche aggregate estratte durante i loop
  all_thresholds <- c()
  all_c_a_tps <- c()
  all_beta_0s <- c()
  all_beta_ints <- c()
  all_tables <- list()
  
  # New outputs used to construct the marginal-effects graph
  all_models <- list()
  all_vcov_matrices <- list()
  all_model_info <- list()
  
  model_id <- 1
  
  f_idx <- 1
  # Loop 1: Itera su tutte le formule fornite dall'utente
  for (f in formulas) { 
    formula <- f
    
    # Liste vuote per salvare ciascun modello temporaneo per poi farne una tabella unica
    fee <- list()
    see <- list()
    i <- 1 
    
    # Statistiche specifiche del modello corrente (formula f), una lista per ogni formuala
    xy <- c() #main variables used
    thresholds <- c()
    p_used <- c() #numero di paesi usati
    beta_0s <- c()
    beta_ints <- c()
    c_a_tps <- c()
    
    # Loop 2: Itera su ogni variante del termine per la CO2
    for (co2_col in emission_var) { 
      
      # Loop 3: Itera su ogni variante del termine per l'ineguaglianza
      for (ineq_col in ine_var) {
        
        # Step 1: Preparazione dataset temporaneo
        # Le colonne target vengono unificate nei nomi "co2" e "inequality" così che
        # le formule in R capiscano automaticamente quali variabili assegnare
        panel_data <- data
        panel_data$co2 <- data[,co2_col]   
        panel_data$inequality <- data[,ineq_col]
        
        # Rimozione dei missing values dalle principali x e y per pulizia primaria
        panel_data <- subset(panel_data, (!is.na(panel_data$co2)) & (!is.na(panel_data$inequality)))    
        
        # Step 2: Filtro il campione
        # Conto i record rimanenti: mantengo solo i paesi che raggiungono o superano min_obs.
        tutti <- plyr::count(panel_data$Country.Code)$x
        vf <- plyr::count(panel_data$Country.Code)$freq >= min_obs
        paesi <- tutti[vf]
        panel_data <- panel_data %>% filter(Country.Code %in% paesi)  
        
        # Filtro dall'anno di inizio
        panel_data <- panel_data %>% filter(Time >= start_year)
        
        # Filtro eventualmente OECD
        if(oecd_only){panel_data <- panel_data %>% filter(oecd == "OECD")}
        
        # Creazione della struttura panel dataset formale
        panel <- pdata.frame(panel_data, index=c("Country.Code","Time"))
        
        # Creo una legenda pulita (LaTeX-friendly) combinando i veri nomi delle variabili
        legenda <- colnames(data[ineq_col])
        
        # Inizializzo le metriche di test
        threshold <- NA; beta_0 <- NA; beta_int <- NA; c_a_tp <- NA; cn <- NA;
        
        if (!weight) {
          # ========== MODELLO SENZA PESI (plm) ==========
          # Modello fixed effetti "twoways" (paese e temporalità)
          fe <- plm(formula, model="within", effect = "twoways", data=panel)
          
          # Compute Small Sample Correction tipico di Stata: N / (N - 1)
          n_countries <- length(unique(panel_data$Country.Code))
          small_sample_correction <- n_countries/(n_countries - 1)
          
          # Estrazione degli errori standard robusti clusterizzati per gruppo (paese)
          
          # Calculate the clustered covariance matrix once.
          # This is exactly the same matrix used previously.
          V_cluster <- small_sample_correction *
            vcovHC(
              fe,
              type = "HC2",
              cluster = "group"
            )
          
          coef_test <- coeftest(
            fe,
            V_cluster
          )
          
          se <- sqrt(
            diag(V_cluster)
          )
          
          # Calcolo numero di paesi intercettati effettivamente dal dataset residuale
          c_used <- as.data.frame(fe$residual)
          c_used <- tibble::rownames_to_column(c_used, "cc")
          c_used$cc <- substr(c_used$cc, 1, 3) 
          cn <- length(unique(c_used$cc))
          
          # Utilizzo logiche dinamiche per rintracciare i coefficienti dell'inequality main effect
          # e dell'interazione (inequality : variabile) che servono per calcolare il turning point
          c_names <- rownames(coef_test)
          int_name <- c_names[grepl(":", c_names) & grepl("inequality", c_names)][1]
          main_name <- c_names[grepl("inequality", c_names) & !grepl(":", c_names)][1]
          if (is.na(main_name)) main_name <- c_names[1]
          
          if (!is.na(int_name) && !is.na(main_name)) {
            beta_0 <- as.numeric(coef_test[main_name, 1])
            beta_int <- as.numeric(coef_test[int_name, 1])
            
            # Calcolo Turning Point: exp(-beta_linear / beta_interaz) - 1
            threshold <- round(exp(-(beta_0)/beta_int)-1, 2)
          }
          
          if (!is.na(threshold)) {
            # Conto quanti dataset superano il TP alla fine del periodo osservato (2018 default max time)
            c_a_tp <- ((panel_data %>% filter(Time==max(Time)) %>% filter(Country.Code %in% c_used$cc) %>%  filter(serv_managr>threshold) %>% tally()) / cn)[1,1]
          }
          
          # Save the fitted model
          all_models[[model_id]] <- fe
          
          # Save the same covariance matrix used for the regression table
          all_vcov_matrices[[model_id]] <- V_cluster
          
          # Save information needed to identify the model later
          all_model_info[[model_id]] <- data.frame(
            model_id = model_id,
            formula_id = f_idx,
            formula_text = paste(
              deparse(formula),
              collapse = " "
            ),
            emission_variable = colnames(data[co2_col]),
            inequality_variable = colnames(data[ineq_col]),
            main_coefficient = main_name,
            interaction_coefficient = int_name,
            beta_inequality = beta_0,
            beta_interaction = beta_int,
            turning_point = threshold,
            countries = cn,
            stringsAsFactors = FALSE
          )
          
        } else {
          # ========== MODELLO PESATO (fixest) ==========
          # Se weight = T, si usa il pacchetto performante "feols" aggiungendo i fissi per stringa
          f_str <- paste(deparse(formula), collapse = " ")
          if (!grepl("\\|", f_str)) {
            f_str <- paste(f_str, "| Country.Code + Time")
            formula <- as.formula(f_str)
          }
          # Pesato per popolazione (pop) -> "weights = ~ pop"
          fe <- fixest::feols(formula, data = panel, weights = ~ pop, cluster = ~ Country.Code)
          # Cluster standard error incorporato automaticamente
          se <- fixest::se(fe)
          
          c_test <- fe$coeftable
          cn <- length(unique(panel$Country.Code))
          
          # Stessa procedura per individuare i coefficienti e derivare il turning point dinamicamente
          c_names <- rownames(c_test)
          int_name <- c_names[grepl(":", c_names) & grepl("inequality", c_names)][1]
          main_name <- c_names[grepl("inequality", c_names) & !grepl(":", c_names)][1]
          if (is.na(main_name)) main_name <- c_names[1]
          
          if (!is.na(int_name) && !is.na(main_name)) {
            beta_0 <- as.numeric(c_test[main_name, 1]) # Column 1 corresponds to "Estimate"
            beta_int <- as.numeric(c_test[int_name, 1])
            threshold <- round(exp(-(beta_0)/beta_int)-1, 2)
          }
          
          if (!is.na(threshold)) {
            c_a_tp <- ((panel_data %>% filter(Time==max(Time)) %>% filter(serv_managr>threshold) %>% tally()) / cn)[1,1]
          }
        }
        
        # Immagazzinamento del loop singolo
        fee[[i]] <- fe
        see[[i]] <- se
        i <- i + 1
        
        xy <- c(xy, legenda)
        thresholds <- c(thresholds, threshold)
        p_used <- c(p_used, cn)
        beta_0s <- c(beta_0s, beta_0)
        beta_ints <- c(beta_ints, beta_int)
        c_a_tps <- c(c_a_tps, c_a_tp)
        
        model_id <- model_id + 1
      }
    }
    
    var_labels <- c(
      "inequality"         = "Inequality (ln)",    
      "l_serv_managr"      = "Tertiarization (ln)",
      "l_gdp_pc"           = "GDP per capita (ln)",
      "l_pop"              = "Population (ln)",
      "l_urban_pop"        = "Urban population (ln)",
      "l_services"         = "Services, value-added (ln)",
      "l_prim_en"          = "Primary enrollment (ln)",
      "l_over_65"          = "Population ages 65 and above (ln)",
      "l_final_hous_cons"  = "Final households consumption per capita (ln)",
      "l_trade_gdp"        = "Trade (ln)",
      "l_fdi_i"            = "FDI inflows (ln)",
      "l_manufacturing"    = "Manufacturing, value-added (ln)",
      "l_industry"         = "Industry, value-added (ln)",
      "str_index"          = "Environmental Policy Stringency Index",
      "l_empl_serv_managr" = "Tertiarization (employment) (ln)",
      "l_serv_empl"        = "Employment in services (ln)",
      "l_serv_indagr"      = "Tertiarization (with industry)",
      "l_serv_managrserv"  = "(Services/ All sectors) (with manufacturing)",
      "l_serv_indagrserv"  = "(Services/ All sectors) (with industry)",
      "inequality:l_serv_managr" = "Inequality*Tertiarization",
      "inequality:l_empl_serv_managr" = "Inequality*Tertiarization (employment)",
      "inequality:l_serv_indagr" = "Inequality*Tertiarization (with industry)",
      "inequality:l_serv_managrserv" = "Inequality*(Services/ All sectors) (with manufacturing)",
      "inequality:l_serv_indagrserv" = "Inequality*(Services/ All sectors) (with industry)"
    )
    
    if (output_type == "latex") {
      coef_names <- unique(unlist(lapply(fee, function(m) names(coef(m)))))
      coef_labels <- sapply(coef_names, function(x) {
        if (x %in% names(var_labels)) return(var_labels[[x]])
        return(x)
      }, USE.NAMES = FALSE)
    } else {
      coef_labels <- NULL
    }
    
    # Label lookups
    ine_labels <- c(
      "l_gini_wid"      = "Gini pre-tax",
      "l_palma"          = "Palma pre-tax",
      "l_r20"            = "20:20 ratio pre-tax",
      "l_p80p100"        = "Top 20\\% pre-tax",
      "l_p90p100"        = "Top 10\\% pre-tax",
      "l_gini_wid_disp"  = "Gini post-tax",
      "l_palma_d"        = "Palma post-tax",
      "l_r20_d"          = "20:20 ratio post-tax",
      "l_p80p100d"       = "Top 20\\% post-tax",
      "l_p90p100d"       = "Top 10\\% post-tax"
    )
    xy_labels <- ifelse(xy %in% names(ine_labels), ine_labels[xy], sapply(xy, function(x) xtable::sanitize(x, type = "latex")))
    
    # Before the stargazer call, build the caption dynamically
    dep_var_name <- colnames(data[co2_col])
    dep_var_label <- if (dep_var_name == "l_co2_pc_wb") {
      "CO$_{2}$ emissions from production (ln)"
    } else if (dep_var_name == "l_eora_cons") {
      "CO$_{2}$ emissions from consumption (ln)"
    } else {
      dep_var_name
    }
    
    # Step 3: Assemblaggio output a schermo/oggetto list.
    if (!weight) {
      tb <- capture.output(
        stargazer::stargazer(fee, se = see,
                             column.labels = xy_labels,
                             dep.var.labels.include = F,
                             dep.var.caption = paste0("\\textit{Dependent variable: ", dep_var_label, "}"),
                             covariate.labels = coef_labels,
                             add.lines=list(
                               c("Country fixed effects", rep("\\checkmark", length(thresholds))),
                               c("Year dummies", rep("\\checkmark", length(thresholds))),
                               c("Turning point", thresholds),
                               c("Countries", p_used),
                               c("\\% of obs. $>$ turning point", round(c_a_tps,2)*100)
                             ),
                             no.space = TRUE,
                             type = output_type)
      )
      if (output_type == "latex") {
        # Add \\[-1.8ex] before column labels
        col_labels_idx <- grep("\\\\cline\\{", tb)
        if (length(col_labels_idx) > 0) {
          for (idx in col_labels_idx) {
            if (idx < length(tb)) {
              tb[idx + 1] <- paste0("\\\\[-1.8ex] ", tb[idx + 1])
            }
          }
        }
        # Add adjustbox
        caption_idx <- grep("\\\\caption\\{", tb)
        if (length(caption_idx) > 0) {
          tb <- append(tb, "\\begin{adjustbox}{width=0.97\\textwidth,center}", after = caption_idx)
        }
        endtab_idx <- grep("\\\\end\\{tabular\\}", tb)
        if (length(endtab_idx) > 0) {
          tb <- append(tb, "\\end{adjustbox}", after = endtab_idx)
        }
      }
      cat(tb, sep = "\n")
      
      all_tables[[f_idx]] <- tb
    } else {
      names(fee) <- xy_labels
      is_tex <- if (output_type == "latex") TRUE else FALSE
      tb <- fixest::etable(fee,
                           extralines = list("Country fixed effects" = rep("\\checkmark", length(thresholds)),
                                             "Year dummies" = rep("\\checkmark", length(thresholds)),
                                             "Turning point" = thresholds,
                                             "Countries" = p_used,
                                             "\\% of obs. $>$ turning point" = round(c_a_tps,2)*100),
                           tex = is_tex)
      all_tables[[f_idx]] <- tb
      print(tb)
    }
    
    f_idx <- f_idx + 1
    all_thresholds <- c(all_thresholds,thresholds)
    all_beta_0s <- c(all_beta_0s,beta_0s)
    all_beta_ints <- c(all_beta_ints,beta_ints)
    all_c_a_tps <- c(all_c_a_tps,c_a_tps)
  }
  
  if (length(all_model_info) > 0) {
    model_info <- do.call(
      rbind,
      all_model_info
    )
    
    rownames(model_info) <- NULL
  } else {
    model_info <- data.frame()
  }
  
  # Restituisci la lista contenente metriche e reports tabellari
  return(list(
    all_thresholds = all_thresholds,
    all_c_a_tps = all_c_a_tps,
    tables = all_tables,
    
    # New outputs    
    models = all_models,    
    vcov_matrices = all_vcov_matrices,    
    model_info = model_info
  ))
}

# calculate_marginal_effects
#
# Calculates the marginal effect of log inequality on log emissions:
#
#   ME(T) = beta_inequality + beta_interaction * T
#
# where T is logged tertiarization, l_serv_managr.
#
# The associated variance is:
#
#   Var[ME(T)] =
#       Var(beta_inequality)
#       + T^2 Var(beta_interaction)
#       + 2 T Cov(beta_inequality, beta_interaction)
#
# The function evaluates this expression over the full observed support of
# logged tertiarization in the estimation sample used by the fitted model.

calculate_marginal_effects <- function(
    model,
    V_cluster,
    model_information,
    moderator = "l_serv_managr",
    grid_points = 200,
    confidence_level = 0.95
) {
  
  # ---------------------------------------------------------------------------
  # 1. Retrieve the exact coefficient names saved by run_regressions2()
  # ---------------------------------------------------------------------------
  
  main_name <- model_information$main_coefficient
  interaction_name <- model_information$interaction_coefficient
  
  
  # Check that the required coefficient names exist.
  if (
    is.na(main_name) ||
    is.na(interaction_name)
  ) {
    stop(
      paste0(
        "The main or interaction coefficient could not be identified ",
        "for model ",
        model_information$model_id,
        "."
      )
    )
  }
  
  
  if (
    !(main_name %in% names(coef(model))) ||
    !(interaction_name %in% names(coef(model)))
  ) {
    stop(
      paste0(
        "The saved coefficient names do not match the coefficients ",
        "in model ",
        model_information$model_id,
        "."
      )
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # 2. Extract the actual estimation sample
  # ---------------------------------------------------------------------------
  #
  # model.frame(model) should contain only observations actually used by plm,
  # including the complete-case selection arising from all model controls.
  # ---------------------------------------------------------------------------
  
  estimation_data <- tryCatch(
    model.frame(model),
    error = function(e) NULL
  )
  
  
  # Fallback for differences between plm versions.
  if (
    is.null(estimation_data) ||
    !(moderator %in% names(estimation_data))
  ) {
    estimation_data <- tryCatch(
      model$model,
      error = function(e) NULL
    )
  }
  
  
  if (
    is.null(estimation_data) ||
    !(moderator %in% names(estimation_data))
  ) {
    stop(
      paste0(
        "Could not find '",
        moderator,
        "' in the estimation sample for model ",
        model_information$model_id,
        "."
      )
    )
  }
  
  
  # Extract observed logged tertiarization values.
  tertiarization_observed_log <- as.numeric(
    estimation_data[[moderator]]
  )
  
  tertiarization_observed_log <-
    tertiarization_observed_log[
      is.finite(tertiarization_observed_log)
    ]
  
  
  if (length(tertiarization_observed_log) < 2) {
    stop(
      paste0(
        "Insufficient valid tertiarization observations for model ",
        model_information$model_id,
        "."
      )
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # 3. Create a grid covering the full observed support
  # ---------------------------------------------------------------------------
  
  support_min_log <- min(
    tertiarization_observed_log
  )
  
  support_max_log <- max(
    tertiarization_observed_log
  )
  
  tertiarization_grid_log <- seq(
    from = support_min_log,
    to = support_max_log,
    length.out = grid_points
  )
  
  
  # ---------------------------------------------------------------------------
  # 4. Retrieve the two relevant coefficients
  # ---------------------------------------------------------------------------
  
  beta_inequality <- as.numeric(
    coef(model)[main_name]
  )
  
  beta_interaction <- as.numeric(
    coef(model)[interaction_name]
  )
  
  
  # ---------------------------------------------------------------------------
  # 5. Calculate the marginal effect
  # ---------------------------------------------------------------------------
  
  marginal_effect <-
    beta_inequality +
    beta_interaction *
    tertiarization_grid_log
  
  
  # ---------------------------------------------------------------------------
  # 6. Calculate the variance and standard error
  # ---------------------------------------------------------------------------
  
  marginal_variance <-
    V_cluster[main_name, main_name] +
    tertiarization_grid_log^2 *
    V_cluster[interaction_name, interaction_name] +
    2 *
    tertiarization_grid_log *
    V_cluster[main_name, interaction_name]
  
  
  # Protect against negligible negative values resulting from numerical
  # rounding, such as -1e-16.
  marginal_variance <- pmax(
    marginal_variance,
    0
  )
  
  marginal_standard_error <- sqrt(
    marginal_variance
  )
  
  
  # ---------------------------------------------------------------------------
  # 7. Calculate pointwise confidence intervals
  # ---------------------------------------------------------------------------
  
  critical_value <- qnorm(
    1 - (1 - confidence_level) / 2
  )
  
  confidence_lower <-
    marginal_effect -
    critical_value *
    marginal_standard_error
  
  confidence_upper <-
    marginal_effect +
    critical_value *
    marginal_standard_error
  
  
  # ---------------------------------------------------------------------------
  # 8. Calculate the turning point
  # ---------------------------------------------------------------------------
  
  turning_point_log <-
    -beta_inequality /
    beta_interaction
  
  turning_point_level <-
    exp(turning_point_log) - 1
  
  
  # ---------------------------------------------------------------------------
  # 9. Convert tertiarization to its original scale
  # ---------------------------------------------------------------------------
  #
  # Your variable is:
  #
  #   l_serv_managr = log(1 + serv_managr)
  #
  # Therefore:
  #
  #   serv_managr = exp(l_serv_managr) - 1
  #
  # The marginal effect must be calculated using the logged value, but the
  # horizontal axis can display the more interpretable original-scale value.
  # ---------------------------------------------------------------------------
  
  tertiarization_grid_level <-
    exp(tertiarization_grid_log) - 1
  
  
  # ---------------------------------------------------------------------------
  # 10. Assign the panel label
  # ---------------------------------------------------------------------------
  
  inequality_variable <-
    model_information$inequality_variable
  
  panel_label <- if (
    inequality_variable == "l_gini_wid"
  ) {
    
    "A. Pre-tax Gini"
    
  } else if (
    inequality_variable == "l_gini_wid_disp"
  ) {
    
    "B. Post-tax Gini"
    
  } else {
    
    inequality_variable
  }
  
  
  # ---------------------------------------------------------------------------
  # 11. Return the graph data
  # ---------------------------------------------------------------------------
  
  data.frame(
    model_id =
      model_information$model_id,
    
    inequality_variable =
      inequality_variable,
    
    panel_label =
      panel_label,
    
    tertiarization_log =
      tertiarization_grid_log,
    
    tertiarization_level =
      tertiarization_grid_level,
    
    marginal_effect =
      marginal_effect,
    
    standard_error =
      marginal_standard_error,
    
    conf_low =
      confidence_lower,
    
    conf_high =
      confidence_upper,
    
    beta_inequality =
      beta_inequality,
    
    beta_interaction =
      beta_interaction,
    
    turning_point_log =
      turning_point_log,
    
    turning_point_level =
      turning_point_level,
    
    support_min_log =
      support_min_log,
    
    support_max_log =
      support_max_log,
    
    support_min_level =
      exp(support_min_log) - 1,
    
    support_max_level =
      exp(support_max_log) - 1,
    
    stringsAsFactors = FALSE
  )
}

#### how many countries above the tp ####

# conta quanti paesi sono sopra il tp (dalla variable usata come interazione, default é serv_managr), in ogni anno

countries_above_tp_per_year <- function(input_data, tp_variable = "serv_managr", tp = th) {
  
  an<-c()
  pa<-c()
  for (i in c(1980:2018)) {
    serv_man18<-input_data %>% filter(Time==i)
    serv_man18<-serv_man18[c("Country.Code","Country.Name",tp_variable)]
    
    paesi_more4 <-serv_man18 %>% filter(.data[[tp_variable]]>=tp) 
    paesi_more4 <- c(paesi_more4$Country.Code)
    
    paesi_less4 <-serv_man18 %>% filter(.data[[tp_variable]]<tp)
    paesi_less4 <- c(paesi_less4$Country.Code)
    
    p <- length(paesi_more4)
    
    an<-c(an,i)
    pa<-c(pa,p)
    
  }
  graph <- data_frame(an,pa) #year - number of country above th
  return(graph)
}

# example
# graph <- countries_above_tp_per_year(input_data = data, tp_variable = "serv_managr", tp = th)

#### frequency countries above the tp ####

# conta per ogni paese, quante volte é al di sopra del tp  (dalla variable usata come interazione, default é serv_managr)
frequancy_country_above_tp <- function(input_data = data, tp_variable = "serv_managr", tp = th) {
  pa<-c()
  ccc<-c()
  input_data <- input_data[!is.na(input_data[tp_variable]),]
  for (c in levels(as.factor(input_data$Country.Code)) ) { 
    dada<-input_data %>% filter(Country.Code==c) 
    cc<-0
    for (i in levels(as.factor(dada$Time)) ) {
      dadada<-dada %>% filter(Time==i)
      if (dadada[tp_variable] >= tp) {cc = cc+1}
    }
    cp<-cc/nrow(dada)
    pa<-c(pa,c)
    ccc<-c(ccc,cp)
  }
  jkl <- data_frame(pa,ccc)
  return(jkl)
}

# example
# geo <- frequancy_country_above_tp(input_data = data, tp = 2.99, tp_variable = "serv_managr") 

#### split the countries in above/below the tp ####

# paesi sopra th in almeno over_in_at_least% of years 
contries_splitted <- function (input_data = data, over_in_at_least, tp_variable = "serv_managr", tp = th, output = "data") {
  
  jkl <- frequancy_country_above_tp(input_data, tp_variable = tp_variable, tp = tp)
  jkl <- jkl %>% filter(ccc>over_in_at_least) #x%
  if (output == "list") {
    paesi_over_list <- jkl$pa
    return (paesi_over_list)
  }
  else {
    paesi_more4 <- input_data %>% filter(Country.Code %in% jkl$pa)
    paesi_less4 <- input_data %>% filter(!Country.Code %in% jkl$pa)
    return(list(paesi_over = paesi_more4, paesi_under = paesi_less4))
  } 
} 

# Example usage:
# paesi_over <- contries_splitted(over_in_at_least = 0.7, input_data = data, tp_variable = "serv_managr", tp = 2.99)$paesi_over
# list_paesi_over <-contries_splitted(over_in_at_least = 0.7, input_data = data, tp_variable = "serv_managr", tp = 2.99, output = "list")

#### mechanism regressions ####

run_mechanism_regressions <- function(data = data, dep_vars, ineq_vars, interaction_var = "over", output_type = "text", vcov = F) {
  
  panel_data <- pdata.frame(data, index=c("Country.Code","Time"))
  
  # Base controls used in all models
  controls <- "l_gdp_pc + l_pop + l_services + l_urban_pop + l_prim_en + l_over_65 + l_final_hous_cons + l_trade_gdp + l_fdi_i"
  
  # Calculate small-sample correction factor for Stata-like robust standard errors
  G <- length(unique(attr(panel_data, "index")[[1]]))
  g <- G / (G - 1)
  
  # Store robust vcov matrices
  vcov_list <- list()
  
  for (y in dep_vars) {
    
    fee <- list()
    see <- list()
    
    for (i in seq_along(ineq_vars)) {
      
      f <- as.formula(paste(y, "~", ineq_vars[i], "*", interaction_var, "+", controls))
      
      mecd <- plm(f, model = "within", effect = "twoways", data = panel_data)
      
      robust_vcov <- g * vcovHC(mecd, type = "HC2", cluster = "group")
      se <- sqrt(diag(robust_vcov))
      
      fee[[i]] <- mecd
      see[[i]] <- se
      vcov_list[[paste(y, ineq_vars[i], sep = "_")]] <- robust_vcov
    }
    
    stargazer::stargazer(
      fee, 
      se = see,
      title = paste("Results for:", y),
      column.labels = c("Gini", "Top 10\\%"),
      add.lines=list(
        c("Country fixed effects", rep("\\checkmark", 2)),
        c("Year dummies", rep("\\checkmark", 2))
      ),
      type = output_type, 
      no.space = TRUE
    )
  }
  if (vcov==T) {
    return(vcov_list)
  }
}

# Example usage:
# ineq_vars_list <- c("gini_wid", "palma", "r20", "p80p100", "p90p100")
# y_list <- c("cons_cf_p0p50", "cons_cf_p90p100", "cons_cf_p99p100", "inv_cf_p90p100", "inv_cf_p99p100")
# run_mechanism_regressions(panel, y_list, ineq_vars_list, output_type = "latex")



plot_residualized_mechanisms <- function(data, dep_vars, ineq_vars, interaction_var = "over", output_dir = ".") {
  
  controls_vars <- c("l_gdp_pc", "l_pop", "l_services", "l_urban_pop", "l_prim_en",
                     "l_over_65", "l_final_hous_cons", "l_trade_gdp", "l_fdi_i")
  controls <- paste(controls_vars, collapse = " + ")
  idx <- c("Country.Code", "Time")
  
  for (y in dep_vars) {
    for (x in ineq_vars) {
      
      vars_needed <- c(idx, y, x, interaction_var, controls_vars)
      
      missing_cols <- vars_needed[!vars_needed %in% colnames(data)]
      if (length(missing_cols) > 0) {
        stop("Colonne non trovate: ", paste(missing_cols, collapse = ", "))
      }
      
      sub_data <- data[complete.cases(data[, vars_needed]), ]
      rownames(sub_data) <- as.character(seq_len(nrow(sub_data)))  # reset a 1:n
      
      f_y <- as.formula(paste(y, "~", controls))
      resid_y <- residuals(plm(f_y, model = "within", effect = "twoways",
                               data = sub_data, index = idx))
      
      f_x <- as.formula(paste(x, "~", controls))
      resid_x <- residuals(plm(f_x, model = "within", effect = "twoways",
                               data = sub_data, index = idx))
      
      common_idx <- intersect(names(resid_y), names(resid_x))
      
      plot_df <- data.frame(
        resid_y = as.numeric(resid_y[common_idx]),
        resid_x = as.numeric(resid_x[common_idx]),
        group   = factor(sub_data[[interaction_var]][as.integer(common_idx)])
      )
      
      var_labels <- c(
        "gini_wid"        = "Gini (pre-tax)",
        "gini_wid_disp"   = "Gini (post-tax)",
        "man_empl"        = "Employment in industry",
        "l_man_empl"      = "Employment in industry",
        "top10_em"          = "Carbon footprint of Top 10%"
      )
      x_label <- ifelse(x %in% names(var_labels), var_labels[x], x)
      y_label <- ifelse(y %in% names(var_labels), var_labels[y], y)
      
      p <- ggplot(plot_df, aes(x = resid_x, y = resid_y, color = group)) +
        geom_point(alpha = 0.25, size = 1.2) +
        geom_smooth(method = "lm", se = TRUE, linewidth = 0.9) +
        scale_color_manual(values = c("0" = "#2166ac", "1" = "#d6604d"),
                           labels = c("0" = "Below threshold", "1" = "Above threshold")) +
        labs(
          x     = x_label,
          y     = y_label,
          color = NULL
        ) +
        theme_minimal(base_size = 12) +
        theme(legend.position = "bottom",
              axis.text = element_text(size = 14))
      
      fname <- file.path(output_dir, paste0("residplot_", y, "_", x, ".png"))
      ggsave(fname, plot = p, width = 6, height = 5, dpi = 300)
      message("Saved: ", fname)
    }
  }
}
