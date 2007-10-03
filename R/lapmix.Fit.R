lapmix.Fit <- function(Y, asym=FALSE, fast=TRUE, two.step=TRUE, w=0.1, V=10, beta=0, gamma=1, alpha=0.1)
{
	if(is(Y,"exprSet") || is(Y,"ExpressionSet"))
		Y <- exprs(Y)
	if(is.matrix(Y))
	{
		n_vect <- apply(!is.nan(Y), 1, sum)
		y_bar <- apply(Y, 1, mean, na.rm=TRUE)
		s_sq <- apply(Y, 1, var, na.rm=TRUE)
		sum_sq <- apply(Y^2, 1, sum, na.rm=TRUE)		
	}
	else if(is.list(Y))
	{
		G <- length(Y)
		y_bar <- NULL
		sum_sq <- NULL
		s_sq <- NULL
		n_vect <- NULL
		for(g in 1:G)
		{
			if(sum(is.nan(Y[[g]])))
				stop('No NaN allowed when data is stored in a list')
			s_sq <- c(s_sq, var(Y[[g]]))
			n_vect <- c(n_vect, length(Y[[g]]))
			y_bar <- c(y_bar, mean(Y[[g]]))
			sum_sq <- c(sum_sq, sum(Y[[g]]^2))
		}
	}
	opt <- lap.maxlike(y_bar=y_bar, s_sq=s_sq, sum_sq=sum_sq, n_vect=n_vect, asym=asym, fast=fast, two.step=two.step, w=w, V=V, beta=beta, gamma=gamma, alpha=alpha)
	estimates <- NULL
	estimates$w <- opt$w
	estimates$V <- opt$V
	estimates$beta <- opt$beta
	estimates$gamma <- opt$gamma
	estimates$alpha <- opt$alpha
	diff <- post_odds(w=estimates$w, V=estimates$V, beta=estimates$beta, gamma=estimates$gamma, alpha=estimates$alpha, y_bar=y_bar, s_sq=s_sq, sum_sq=sum_sq, n_vect=n_vect)
	list(lap.lods=log(diff$odds), prob=diff$prob, med.number=diff$med.number, estimates=estimates, code=opt$code, M=y_bar, s_sq=s_sq, nb.rep=n_vect)
}
