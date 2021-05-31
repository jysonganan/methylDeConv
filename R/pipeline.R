#'Deconvolution Pipeline
#'
#'Deconvolution Pipeline
#'@param input_methyl Methylation array for mixture samples. The rows correpsond to the probes/CpGs;
#'the columns correspond to the samples.
#'@param normalized If TRUE, the methylation array is the beta value matrix. Default value: FALSE.
#'@param array_type The array type can be 450k or EPIC. Default value: "450k".
#'@param feature_selection Different feature selection methods such as "oneVsAllttest" (Default value), "oneVsAllLimma",
#'"pairwiseLimma", "pairwiseGlmnet", "multiGlmnet", "glmnetpreselect", "RFpreselect", "OptVariables".
#'@param deconv_algorithm Different deconvolution algorithms such as "Houseman" (default value), "RPC", "CBS", "MethylResolver",
#'"MCP-counter", "ssGSEA", "ESTIMATE".
#'@param tissue Different tissue types such as "general" (default value) which correspond to blood tissue and general epithelial
#'tissue; "brain" correpsonds to brain tissues.
#'@param extend_reference Default: TRUE. If TRUE, extend the reference library.
#'@param custom_probes Default: NULL. The user-defined probe lists.
#'@return The deconvolution results (cell type proportions) with rows as the samples and the columns as cell types.
#'@export


MethylDeconv_pipeline <- function(input_methyl, normalized = FALSE, array_type = "450k", feature_selection = "oneVsAllttest",
                         deconv_algorithm = "Houseman", tissue = "general", extend_reference = TRUE,
                         custom_probes = NULL){
  if(normalized){
    # beta value input
    MethylDeconv_normalized(input_methyl = input_methyl, array_type = array_type, feature_selection = feature_selection,
                            deconv_algorithm = deconv_algorithm, tissue = tissue, extend_reference = extend_reference,
                            custom_probes = custom_probes)
  }

  if (is(object, "RGChannelSet")){
    GRset.normalized <- getBeta(preprocessNoob(input_methyl, dyeMethod = "single"))
    MethylDeconv_normalized(input_methyl = GRset.normalized, array_type = array_type, feature_selection = feature_selection,
                            deconv_algorithm = deconv_algorithm, tissue = tissue, extend_reference = extend_reference,
                            custom_probes = custom_probes)
  }
  }




MethylDeconv_normalized <- function(input_methyl, array_type, feature_selection, deconv_algorithm,
                                    tissue, extend_reference, custom_probes){
  if (tissue == "general" && array_type == "EPIC"){
    ref <- build_reference_EPIC(extend = extend_reference)
    ref_betamatrix <- ref$ref_betamatrix
    ref_phenotype <- ref$ref_phenotype
  }

  if (tissue == "general" && array_type == "450k"){
    ref <- build_reference_450k(extend = extend_reference)
    ref_betamatrix <- ref$ref_betamatrix
    ref_phenotype <- ref$ref_phenotype
  }

  if (tissue == "brain" && array_type == "EPIC"){
    ref <- build_reference_EPIC_neuron()
    ref_betamatrix <- ref$ref_betamatrix
    ref_phenotype <- ref$ref_phenotype
  }

  if (tissue == "brain" && array_type == "450k"){
    ref <- build_reference_450k_neuron()
    ref_betamatrix <- ref$ref_betamatrix
    ref_phenotype <- ref$ref_phenotype
  }

  ## average reference profiles
  compTable <- ref_compTable(ref_betamatrix, ref_phenotype)
  ntypes <- length(unique(ref_phenotype))
  compTable <- compTable[,3:(2 + ntypes)]

  ## feature selection
  if (is.null(custom_probes)){
    if (feature_selection == "oneVsAllttest"){
      probes_select <- ref_probe_selection_oneVsAllttest(ref_betamatrix, ref_phenotype)
    }
    else if (feature_selection == "oneVsAllLimma"){
      probes_select <- ref_probe_selection_oneVsAllLimma(ref_betamatrix, ref_phenotype)
    }
    else if (feature_selection == "pairwiseLimma"){
      probes_select <- ref_probe_selection_pairwiseLimma(ref_betamatrix, ref_phenotype)
    }
    else if (feature_selection == "pairwiseGlmnet"){
      probes_select <- ref_probe_selection_pairwiseGlmnet_cv(ref_betamatrix, ref_phenotype)
    }
    else if (feature_selection == "multiGlmnet"){
      probes_select <- ref_probe_selection_multiclassGlmnet_cv(ref_betamatrix, ref_phenotype)
    }
    else if (feature_selection == "glmnetpreselect"){
      probes_select <- ref_probe_selection_twoStage(ref_betamatrix, ref_phenotype, ml_model = "elastic net")
    }
    else if (feature_selection == "RFpreselect"){
      probes_select <- ref_probe_selection_twoStage(ref_betamatrix, ref_phenotype, ml_model = "RF")
    }
    else if (feature_selection == "OptVariables"){
      probes_select <- ref_probe_selection_twoStage(ref_betamatrix, ref_phenotype, ml_model = "rfe")
    }
  }
  else{
    probes_select = custom_probes
  }

  ## deconvolution
  if (deconv_algorithm == "Houseman"){
    deconv_res = Houseman_project(input_methyl, compTable, probes_select)
  }
  else if (deconv_algorithm == "RPC"){
    deconv_res = RPC(input_methyl, compTable, probes_select)
  }
  else if (deconv_algorithm == "CBS"){
    deconv_res = CBS(input_methyl, compTable, probes_select)
  }
  else if (deconv_algorithm == "MethylResolver"){
    deconv_res = MethylResolver(input_methyl, compTable)
  }
  else if (deconv_algorithm == "MCP-counter"){
    deconv_res = enrichment_score(input_methyl, ref_betamatrix, ref_phenotype, method = "MCP-counter")
  }
  else if (deconv_algorithm == "ssGSEA"){
    deconv_res = enrichment_score(input_methyl, ref_betamatrix, ref_phenotype, method = "ssGSEA")
  }
  else if (deconv_algorithm == "ESTIMATE"){
    deconv_res = enrichment_score(input_methyl, ref_betamatrix, ref_phenotype, method = "ESTIMATE")
  }

  return(deconv_res)
}










