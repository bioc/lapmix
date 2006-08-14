lap.maxlike <- function(y_bar, s_sq, sum_sq, n_vect, w=0.1, V=10, beta=0, gamma=1, alpha=0.1, asym=FALSE, fast=!asym, two.step=TRUE)
{
	out <- list(NULL)
	if(two.step)
	{
		res1 <- nlm(var.loglike, p=c(log(gamma), log(alpha)), s_sq=s_sq, n_vect=n_vect)
		phi_hat <- exp(res1$estimate)
		if(res1$code >= 4)
			print('convergence failed in fisrt stage of hyperparameter estimation')
		if(!asym)
		{
			res <- nlm(marginal.loglike, p=c(log(w/(1-w)), log(V)), phi_hat=phi_hat, y_bar=y_bar, s_sq=s_sq, sum_sq=sum_sq, n_vect=n_vect, fast=fast)
			l_Sigma_hat <- res$estimate
			if(res$code >= 4)
				print('convergence failed in second stage of hyperparameter estimation')
			out <- list(w=exp(l_Sigma_hat[1])/(1+exp(l_Sigma_hat[1])), V=exp(l_Sigma_hat[2]), beta=0, gamma=phi_hat[1], alpha=phi_hat[2], code=c(res1$code, res$code))
		}
		else
		{
			res <- nlm(marginal.loglike, p=c(log(w/(1-w)), log(V), log((1+beta)/(1-beta))), phi_hat=phi_hat, y_bar=y_bar, s_sq=s_sq, sum_sq=sum_sq, n_vect=n_vect, fast=fast)
			l_Sigma_hat <- res$estimate
			if(res$code >= 4)
				print('convergence failed in second stage of hyperparameter estimation')
			out <- list(w=exp(l_Sigma_hat[1])/(1+exp(l_Sigma_hat[1])), V=exp(l_Sigma_hat[2]), beta=2*exp(l_Sigma_hat[3])/(1+exp(l_Sigma_hat[3]))-1, gamma=phi_hat[1], alpha=phi_hat[2], code=c(res1$code, res$code))
		}
	}
	else
	{
		if(!asym)
		{
			res <- nlm(marginal.loglike, p=c(log(w/(1-w)), log(V), log(gamma), log(alpha)), y_bar=y_bar, s_sq=s_sq, sum_sq=sum_sq, n_vect=n_vect, fast=fast)
			l_Sigma_hat <- res$estimate
			if(res$code >= 4)
				print('convergence failed in hyperparameter estimation')
			out <- list(w=exp(l_Sigma_hat[1])/(1+exp(l_Sigma_hat[1])), V=exp(l_Sigma_hat[2]), beta=0, gamma=exp(l_Sigma_hat[3]), alpha=exp(l_Sigma_hat[4]), code=res$code)
		}
		else
		{
			res <- nlm(marginal.loglike, p=c(log(w/(1-w)), log(V), log((1+beta)/(1-beta)), log(gamma), log(alpha)), y_bar=y_bar, s_sq=s_sq, sum_sq=sum_sq, n_vect=n_vect, fast=fast)
			l_Sigma_hat <- res$estimate
			if(res$code >= 4)
				print('convergence failed in hyperparameter estimation')
			out <- list(w=exp(l_Sigma_hat[1])/(1+exp(l_Sigma_hat[1])), V=exp(l_Sigma_hat[2]), beta=2*exp(l_Sigma_hat[3])/(1+exp(l_Sigma_hat[3]))-1, gamma=exp(l_Sigma_hat[4]), alpha=exp(l_Sigma_hat[5]), code=res$code)
		}
	}
	out
}

