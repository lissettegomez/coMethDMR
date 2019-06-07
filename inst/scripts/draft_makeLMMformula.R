.MakeLmmFormula <- function(pred1, pred2 = NULL,
                            modelType = c("randCoeffMixed", "mixed")){

  modelType <- match.arg(modelType)
  baseMod_char <- "Mvals ~ (1|Sample)"


  ###  pred2 Formula String  ###
  if(modelType == "randCoeffMixed"){
    # If we need the predictor for the random coefficient component, i.e.
    #   (<pred2>|ProbeID), then we check that it is a length-1 string
    if(is.character(pred2) && length(pred2) == 1L){

      # If so, then we create the random coefficient component
      pred2Form_char <- paste0("(", pred2, "|ProbeID)")
      # Combine with the base model
      rcMod_char <- paste(baseMod_char, pred2Form_char, sep = " + ")

    } else {
      stop("pred2 must be a length-1 character vector")
    }

  } else {
    # If we aren't using the random coefficient model, then the base model is
    #   fine by itself. We will print a warning if the user supplied a value to
    #   pred2

    if(!is.null(pred2)){
      warning(
        paste0(
  "You have supplied a <FEATURE TYPE HERE> but specified a basic mixed model.\n  ",
  pred2, " will not be included in the model unless you select the random
  coefficient mixed model type."
        )
      )
    }
    rcMod_char <- baseMod_char

  }


  ###  pred1 Formula string  ###
  if(is.character(pred1)){
    pred1Form_char <- paste(pred1, collapse = " + ")
  } else {
    stop("pred1 must be a character vector")
  }


  ###  Combine the Formula Strings  ###
  fullMod_char <- paste(
    rcMod_char,
    pred1Form_char,
    sep = " + "
  )
  fullMod_char

}


###  Test  ###
# mixed model
paste0("x", 1:5)
.MakeLmmFormula(
  pred1 = paste0("x", 1:5),
  modelType = "mixed"
)
.MakeLmmFormula(
  pred1 = c(paste0("x", 1:5), "x2 * x4"),
  modelType = "mixed"
)
.MakeLmmFormula(
  pred1 = lapply(1:5, function(i) paste0("x", i)),
  modelType = "mixed"
)

# RC mixed model
.MakeLmmFormula(
  pred1 = paste0("x", 1:5),
  pred2 = "age",
  modelType = "rand"
)
.MakeLmmFormula(
  pred1 = paste0("x", 1:5),
  pred2 = "age",
  modelType = "mixed"
)
.MakeLmmFormula(
  pred1 = paste0("x", 1:5),
  pred2 = c("age", "sex"),
  modelType = "rand"
)


###  Test with formula() wrapping  ###
formula(
  .MakeLmmFormula(
    pred1 = paste0("x", 1:5),
    modelType = "mixed"
  )
)

formula(
  .MakeLmmFormula(
    pred1 = paste0("x", 1:5),
    pred2 = "age",
    modelType = "rand"
  )
)
