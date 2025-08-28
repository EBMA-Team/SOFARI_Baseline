```{r}

#| label: setup-packages

# Enable parallel processing
future::plan(multisession)
```

```{r}
#| label: configuration-setup

# Configuration object for all patterns and rules
create_config <- function() {
  list(
    # Text standardization patterns
    text_patterns = list(
      punctuation_normalize = c("\\.|\\. |\\: |\\, |w\\/" = ";"),
      hash_replace = c("\\#" = "fracture"),
      multiple_semicolon = c(";+" = ";"),
      conjunction_replace = c("\\bwith\\b|\\band\\b|\\bas well as\\b" = ";"),
      whitespace_normalize = c("\\s+" = " ")
    ),
    
    # Classification patterns for pathology
    anatomical_patterns = list(
      ankle = "ankle|tibiotalar|\\bplafond\\b|dome|malleol*|weber|achilles|tendo-achilles|fibula|\\btibia\\b|gutter|perone*|syndesmo|gastrocnemius|talo-fibular|talofibular|calcaneofibular|calcaneo-fibular|gutter|(lateral|medial)\\s+ligament|tibia|deltoid|compartment.+syndrome",
      rearfoot = "\\btarsal\\b|\\btalar(?!\\s+dome)|talus|talonavic*|plantar|rearfoot|hindfoot|trigonum|hindfeet|tarsi|calcaneus|\\bcalcaneal\\b|heel|subtalar",
      midfoot = "metatarsus|tarsometatarsal|tarso-metatarsal|\\bmetatarsal\\b|talonavicular|navicular|cuneiform|lisfranc|cuboid|midfoot|jones|chopart",
      forefoot = "digit*|morton*|metatarsophalangeal|phalange*|phalanx|hallux|nail|forefoot|forefeet|bunionette|hallucis|onychomycosis|paronychia|bunion|hammertoe|claw|sesamoid",
      foot = "\\bfeet\\b|\\bfoot\\b|cavovarus|equinovarus|equinus|pes|charcot|footdrop|neuropathy"
    ),
    
    pathological_patterns = list(
      arthritis = "psoria|arthritis|osteoarthritis|rheumatoid|gout|erosion|arthropathy",
      injury = "injury|injuries|axial|impact|crush|rotation|inversion|forced|accident|tear|torn|ruptur|avulsion|fracture|defect|(osteochondral|chondral|cartilage).+lesion|sprain|haemarthrosis|disruption|wound|laceration|penetrating|hernia|maisonneuve",
      deformity = "malalignment|deformit|angulation|contracture|contraction|\\bvalgus\\b|varus|planovalgus|valgoplanus|dysfunction|extension|adductus|crossover|hammer|claw|bunionette|bunion|interphalangeus|(relatively|significantly).+long|relative.+long",
      metatarsalgia = "metatarsalgia|forefoot.+overload",
      soft_tissue_disorder = "tenosynovitis|enthesopathy|teno-synovitis|tendinopathy|tendinitis|tendonitis|tendinosis|fasciosis|fasciitis|sesamoiditis|arthrofibro|scar|tibialis posterior.+dysfunction|dysfunction tibialis posterior",
      growth = "cyst|ganglion|neuroma|malformation|fibroma|tumour|accessory|ingrown|in_grown|coalition|(?<!(?:chondral|osteochondral|cartilage)\\s)\\blesion\\b|xanthomas|osteoma|gioma|schwannoma|chondroma|lump|villonodular|callosity|corn|mass|\\b(non|mal|delayed)[-]?union|pseudo-articulation|bone.+loss|exostosis|spur|osteophyte|onychogryphosis|bipartite|neuroma|chondromatosis",
      neural = "foot.+drop|footdrop|nerve|neuropathy|neural|sensory|charcot|motor|pain.+syndrome|tunnel.+syndrome|neuropathic|denervation",
      infection = "infect|osteomyelitis|cellulitis|onychomycosis|ulcer|paronychia",
      impingement = "impingement|stiffness|os.+trigonum",
      instability = "disloc|unstable|sublux|instability|talar.+shift|widening|maisonneuve"
    ),
    
    # Exclusion patterns
    exclusion_patterns = list(
      ankle_exclusions = "(lateral|medial)\\s+collateral\\s+ligament",
      numeric_exclusions = "\\d+(?!(?:st|nd|rd|th)\\b)|(left|right)",
      negation_patterns = "(?<!ab?)normal|(?<!(non|mal)-?)|nil.+pathology|\\bheal\\b|reduced|non-tender|unremarkable"
    )
  )
}
```

```{r}
#| label: load-target-terms

# Memoized function to load target terms from Google Sheets
load_target_terms <- memoise::memoise(function(sheet_url, sheet_name = "DiagTerm", range = "A1:C") {
  message(glue::glue("Loading target terms from {sheet_name}..."))
  
  terms <- googlesheets4::range_read(
    ss = sheet_url,
    sheet = sheet_name,
    range = range,
    col_names = TRUE,
    trim_ws = TRUE
  ) |> 
    mutate(
      target_term = paste0("\\b", str_escape(Term), "\\b")
    )
  
  list(
    terms = terms,
    pattern = str_c(terms$target_term, collapse = "|")
  )
})
```

```{r}
#| label: text-standardization

# Unified text standardization function
standardize_text <- function(df, config) {
  message(glue::glue("Standardizing text for {nrow(df)} records..."))
  
  # Input validation
  required_cols <- c("TreatmentID", "DiagnosisRawFinal", "DiagnosisRawPrelim")
  missing_cols <- setdiff(required_cols, names(df))
  
  if (length(missing_cols) > 0) {
    stop(glue::glue("Missing required columns: {paste(missing_cols, collapse = ', ')}"))
  }
  
  df |> 
    select(TreatmentID, DiagnosisRawFinal, DiagnosisRawPrelim) |> 
    unite("diagnosis_raw", c(DiagnosisRawFinal, DiagnosisRawPrelim), 
          na.rm = TRUE, remove = FALSE, sep = "; ") |> 
    filter(str_count(str_to_lower(diagnosis_raw), "") > 1) |> 
    mutate(
      # Apply all text patterns sequentially
      diagnosis_clean = reduce(
        config$text_patterns,
        ~ str_replace_all(.x, names(.y), .y),
        .init = str_squish(diagnosis_raw)
      ),
      # Final cleanup
      diagnosis_clean = str_trim(str_remove_all(diagnosis_clean, "^;|;$"))
    ) |> arrange(TreatmentID)
}
```

```{r}
#| label: tokenize-and-classify

# Combined tokenization and term replacement
tokenize_and_classify <- function(df, target_terms, config, batch_size = 1000) {
  message(glue::glue("Tokenizing and classifying {nrow(df)} diagnosis records..."))
  
  # Create replacement function
  replace_function <- function(term) {
    match <- filter(target_terms$terms, Term == term)
    if (nrow(match) == 1) match$ReplaceTerm else term
  }
  
  # Add batch and global IDs
  df_with_ids <- df |> 
    mutate(
      batch_id = ceiling(row_number() / batch_size),
      global_id = row_number()
    )
  
  # Process in batches
  result_list <- df_with_ids |> 
    group_split(batch_id) |> 
    map_dfr(~ {
      batch_num <- unique(.x$batch_id)
      message(glue::glue("Processing batch {batch_num}..."))
      
      .x |> 
        # Split diagnosis into sequences
        mutate(
          diagnosis_lower = str_replace_all(
            str_to_lower(diagnosis_clean),
            "\\bwith\\b|\\band\\b|\\bas well as\\b", ";"
          )
        ) |> 
        separate_rows(diagnosis_lower, sep = ";") |> 
        mutate(
          diagnosis_lower = str_trim(diagnosis_lower),
          SequenceID = row_number()
        ) |> 
        filter(nchar(diagnosis_lower) > 0) |> 
        
        # Tokenize
        unnest_tokens(
          output = term,
          input = diagnosis_lower,
          token = "regex",
          pattern = "\\s+",
          to_lower = TRUE,
          drop = FALSE
        ) |> 
        
        # Remove stop words
        anti_join(get_stopwords(), by = c("term" = "word")) |> 
        
        # Apply target term replacements
        mutate(
          term_replaced = map_chr(term, replace_function)
        ) |> 
        
        # Filter unwanted patterns
        filter(
          !str_detect(term_replaced, "\\d+(?!(?:st|nd|rd|th)\\b)|(left|right)")
        ) |> 
        
        # Add sequence information
        group_by(TreatmentID, SequenceID) |> 
        mutate(
          term_position = row_number(),
          term_count = n()
        ) |> 
        ungroup() |> dplyr::arrange(TreatmentID)
    })
  
  # Return with proper ordering
  result_list |> 
    arrange(global_id, term_position)
}
```

```{r}
#| label: categorize-pathology

# Streamlined pathology categorization
categorize_pathology <- function(df, config) {
  message(glue::glue("Categorizing pathology for {nrow(df)} terms..."))
  
  # Input validation
  if (!"term_replaced" %in% names(df)) {
    stop("Missing required column 'term_replaced'")
  }
  
  df |> 
    mutate(
      # Anatomical classifications
      Ankle = as.numeric(str_detect(term_replaced, config$anatomical_patterns$ankle) & 
                           !str_detect(term_replaced, "(lateral|medial)\\s+collateral\\s+ligament")),
      Rearfoot = as.numeric(str_detect(term_replaced, config$anatomical_patterns$rearfoot)),
      Midfoot = as.numeric(str_detect(term_replaced, config$anatomical_patterns$midfoot)),
      Forefoot = as.numeric(str_detect(term_replaced, config$anatomical_patterns$forefoot)),
      Foot = as.numeric(str_detect(term_replaced, config$anatomical_patterns$foot)),
      
      # Pathological classifications  
      Arthritis = as.numeric(str_detect(term_replaced, config$pathological_patterns$arthritis)),
      Injury = as.numeric(str_detect(term_replaced, config$pathological_patterns$injury)),
      Deformity = as.numeric(str_detect(term_replaced, config$pathological_patterns$deformity)),
      Metatarsalgia = as.numeric(str_detect(term_replaced, config$pathological_patterns$metatarsalgia)),
      Soft_tissue_disorder = as.numeric(str_detect(term_replaced, config$pathological_patterns$soft_tissue_disorder)),
      Growth = as.numeric(str_detect(term_replaced, config$pathological_patterns$growth)),
      Neural = as.numeric(str_detect(term_replaced, config$pathological_patterns$neural)),
      Infection = as.numeric(str_detect(term_replaced, config$pathological_patterns$infection)),
      Impingement = as.numeric(str_detect(term_replaced, config$pathological_patterns$impingement)),
      Instability = as.numeric(str_detect(term_replaced, config$pathological_patterns$instability)),
      
      # Calculate composite scores
      PathologySum = Arthritis + Injury + Deformity + Metatarsalgia + Soft_tissue_disorder + 
        Growth + Neural + Infection + Impingement + Instability,
      AnatomySum = Ankle + Rearfoot + Midfoot + Forefoot + Foot,
      AnatPathSum = PathologySum + AnatomySum,
      
      # Final classifications
      Other = as.numeric(str_detect(term_replaced, "foreign|body|material|object")),
      # Define negative pathology indicators (things that suggest no problems)  
      NegatePathology = as.numeric(str_detect(term_replaced, "\\b(normal(?<!abnormal)|nil.*pathology|heal|healed|healing|reduced|non-tender|unremarkable)\\b"))
    )
}
```

```{r}
#| label: concatenate-results

# Simplified result concatenation - preserving all classification columns
concatenate_results <- function(df) {
  message(glue::glue("Concatenating results for {length(unique(df$TreatmentID))} treatments..."))
  
  # First, aggregate by sequence to get term combinations and max classifications
  sequence_level <- df |> 
    group_by(TreatmentID, SequenceID) |> 
    summarise(
      CombinedTerm = str_c(term_replaced, collapse = " "),
      # Preserve all classification columns by taking the maximum value for each sequence
      across(c(Ankle:Other, NegatePathology), ~ max(.x, na.rm = TRUE)),
      PathologySum = sum(c_across(Arthritis:Instability)) + sum(c_across(Other), na.rm = TRUE),
      # Sum the anatomy categories (Ankle:Foot) 
      AnatomySum = sum(c_across(Ankle:Foot), na.rm = TRUE),
      .groups = "keep"
    ) |>
    arrange(TreatmentID, SequenceID)
  
  # Filter out nopathology sequences
  
  # Then create the final output with row-level data (similar to your DiagTerm structure)
  sequence_level |> 
    group_by(TreatmentID) |> 
    mutate(
      CumulativeDiagnosis = accumulate(
        CombinedTerm,
        ~ if_else(is.na(.x), .y, paste(.y, .x, sep = "; "))
      ) |> last(),
      SequenceRow = row_number()
    ) |> 
    ungroup() |> 
    select(
      TreatmentID,
      SequenceID,
      CombinedTerm,
      Ankle:NegatePathology,
      PathologySum,
      AnatomySum,
      AnatPathSum,
      CumulativeDiagnosis,
      SequenceRow
    )
}

```






```{r}
#| label: main-pipeline

# Main processing pipeline
process_clinical_diagnoses <- function(
    data, 
    target_terms_url,
    batch_size = 1000,
    use_parallel = TRUE
) {
  
  # Setup
  config <- create_config()
  
  if (use_parallel) {
    plan(multisession)
  }
  
  # Input validation
  message(glue::glue("Starting clinical diagnosis processing pipeline..."))
  message(glue::glue("Input data: {nrow(data)} records"))
  
  tryCatch({
    # Load target terms
    target_terms <- load_target_terms(target_terms_url)
    message(glue::glue("Loaded {nrow(target_terms$terms)} target terms"))
    
    # Main processing pipeline
    result <- data |> 
      standardize_text(config) |> 
      tokenize_and_classify(target_terms, config, batch_size) |> 
      categorize_pathology(config) |> 
      concatenate_results()
    
    message("Processing completed successfully!")
    message(glue::glue("Output: {nrow(result)} processed records"))
    
    return(result)
    
  }, error = function(e) {
    message(glue::glue("Processing failed: {e$message}"))
    stop(e)
  })
}
```