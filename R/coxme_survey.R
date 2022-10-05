
#' functions needed to make svycoxme go.
#'
#' svycoxme
#'
#'
#'
#'


#' Need a function that takes the random effect information from coxme and converts it into an offset.
#' This function will need to handle interaction terms from the random effects.


#' Generic method
#'
#' rescale should always be true in my opinion. Coxme often fails with unscaled weights. When it runs, the estimates it produces are really wrong.
#'

svycoxme <- function(formula, design, subset=NULL, ...){
  .svycheck(design)
  UseMethod("svycoxme", design)
}

svycoxme.DBIsvydesign <- function(formula, design, subset, ...){

  warning("coxme has not been implemented for class = \"DBIsvydesign\"")

  # return something?

}

svycoxme.survey.design <- function(formula, design, subset, ...){

  warning("coxme has not been implemented for class = \"survey.design\"")

  # return something?

}

svycoxme.svyrep.design <- function(formula, design, subset, ...){

  warning("coxme has not been implemented for \"class = svyrep.design\" ")

  # return something?

}

svycontrast <- function(stat, contrasts, add = FALSE, ...){

  UseMethod("svycontrast")

}

svycontrast.svycoxme <- function(){

  warning("svycontrast has not been implemented for \"class = svycoxme\" ")

  # return something?

}

svycontrast.svrepcoxme <- function(){

  warning("svycontrast has not been implemented for \"class = svrepcoxme\" ")

  # return something?

}

#' check next release for the fix to the rescale issues.
svycoxme.survey.design<-function(formula,design, subset=NULL, rescale=TRUE, ...){
  subset<-substitute(subset)
  subset<-eval(subset, model.frame(design),parent.frame())
  if (!is.null(subset))
    design<-design[subset,]

  if(any(weights(design)<0)) stop("weights must be non-negative")

  data<-model.frame(design)

  g<-match.call()
  g$formula<-eval.parent(g$formula)
  g$design<-NULL
  g$var<-NULL
  g$rescale <- NULL
  if (is.null(g$weights))
    g$weights<-quote(.survey.prob.weights)
  else
    g$weights<-bquote(.survey.prob.weights*.(g$weights))
  ## g[[1]]<-quote(coxph)
  g[[1]] <- quote(coxme)
  g$data<-quote(data)
  g$subset<-quote(.survey.prob.weights>0)
  g$model <- TRUE

  ##need to rescale weights for stability
  ## unless the user doesn't want to
  if (rescale)
    data$.survey.prob.weights<-(1/design$prob)/mean(1/design$prob)
  if (!all(all.vars(formula) %in% names(data)))
    stop("all variables must be in design= argument")

  g_ph <- g

  g<-with(list(data=data), eval(g))



  g_ph[[1]] <- quote(coxph)
  g_ph$offset <-




  ## not needed.
  # if (inherits(g, "coxph.penal"))
  #   warning("svycoxph does not support penalised terms")

  g$call<-match.call()
  g$call[[1]]<-as.name(.Generic)
  g$printcall<-sys.call(-1)
  g$printcall[[1]]<-as.name(.Generic)
  # class(g)<-c("svycoxph", class(g))
  class(g)<-c("svycoxme", class(g))
  g$survey.design<-design

  nas<-g$na.action
  if (length(nas))
    design<-design[-nas,]

  # Problem: there is no residuals.coxme method.
  # Solution: fit a coxph model with the random effects as offsets.

  dbeta.subset<-resid(g,"dfbeta",weighted=TRUE)
  if (nrow(design)==NROW(dbeta.subset)){
    dbeta<-as.matrix(dbeta.subset)
  } else {
    dbeta<-matrix(0,ncol=NCOL(dbeta.subset),nrow=nrow(design))
    dbeta[is.finite(design$prob),]<-dbeta.subset
  }
  g$inv.info<-g$var

  if (inherits(design,"survey.design2"))
    g$var<-svyrecvar(dbeta, design$cluster,
                     design$strata, design$fpc,
                     postStrata=design$postStrata)
  else if (inherits(design, "twophase"))
    g$var<-twophasevar(dbeta, design)
  else if(inherits(design, "twophase2"))
    g$var<-twophase2var(dbeta, design)
  else if(inherits(design, "pps"))
    g$var<-ppsvar(dbeta,design)
  else
    g$var<-svyCprod(dbeta, design$strata,
                    design$cluster[[1]], design$fpc,design$nPSU,
                    design$certainty,design$postStrata)

  g$wald.test<-coef(g)%*%solve(g$var,coef(g))
  g$ll<-g$loglik
  g$loglik<-NULL
  g$rscore<-NULL
  g$score<-NA
  g$degf.resid<-degf(design)-length(coef(g)[!is.na(coef(g))])+1

  g
}




