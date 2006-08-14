var.loglike <- function(l_phi, s_sq, n_vect)
{
	gamma <- exp(l_phi[1])
	alpha <- exp(l_phi[2])
	if(gamma < 10e-20)
		loglike <- Inf
	else
		loglike <- -sum(log(alpha*gamma*df(alpha*gamma*s_sq, df1=n_vect-1, df2=2*gamma)))
	loglike
}

marginal.loglike <- function(l_Sigma, phi_hat=NULL, y_bar, s_sq, sum_sq, n_vect, fast=TRUE)
{
	full <- FALSE
	if(length(l_Sigma) == 2)   ## symmetric model,  fixed gamma and alpha
	{	
		w <- exp(l_Sigma[1])/(1+exp(l_Sigma[1]))
		V <- exp(l_Sigma[2])
		beta <- 0
		gamma <- phi_hat[1]
		alpha <- phi_hat[2]
	}
	else if(length(l_Sigma) == 3)  ## asymmetric model,  fixed gamma and alpha
	{
		w <- exp(l_Sigma[1])/(1+exp(l_Sigma[1]))
		V <- exp(l_Sigma[2])
		beta <- 2*exp(l_Sigma[3])/(1+exp(l_Sigma[3]))-1
		gamma <- phi_hat[1]
		alpha <- phi_hat[2]
	}
	else if(length(l_Sigma) == 4)  ## symmetric model
	{
		w <- exp(l_Sigma[1])/(1+exp(l_Sigma[1]))
		V <- exp(l_Sigma[2])
		beta <- 0
		gamma <- exp(l_Sigma[3])
		alpha <- exp(l_Sigma[4])
		full <- TRUE
	}
	else if(length(l_Sigma) == 5)  ## asymmetric model
	{
		w <- exp(l_Sigma[1])/(1+exp(l_Sigma[1]))
		V <- exp(l_Sigma[2])
		beta <- 2*exp(l_Sigma[3])/(1+exp(l_Sigma[3]))-1
		gamma <- exp(l_Sigma[4])
		alpha <- exp(l_Sigma[5])
		full <- TRUE
	}
	nu <- 2*gamma+n_vect+1
	Vp <- V*(1+beta)
	Vm <- V*(1-beta)
	Bp <- y_bar-1/(Vp*n_vect)
	Bm <- y_bar+1/(Vm*n_vect)
	Cp <- ((n_vect-1)*s_sq+2*y_bar/Vp-1/(n_vect*Vp^2)+2/alpha)/n_vect
	Cm <- ((n_vect-1)*s_sq-2*y_bar/Vm-1/(n_vect*Vm^2)+2/alpha)/n_vect
	Dp <- sqrt(abs(nu/Cp))
	Dm <- sqrt(abs(nu/Cm))
	E <- sqrt(nu*pi)*gamma(nu/2)/gamma((nu+1)/2)
	Ip <- E*(Cp*n_vect/2)^(-(nu+1)/2)*pt(Dp*Bp, df=nu)/Dp
	Im <- E*(Cm*n_vect/2)^(-(nu+1)/2)*pt(-Dm*Bm, df=nu)/Dm
	negp <- which(Cp < 0)
	negm <- which(Cm < 0)
	if(fast)
	{
		Ip[negp] <- 0
		Im[negm] <- 0
	}
	else
	{
		for(i in negp)
			Ip[i] <- integrate(integrand, y_bar=y_bar[i], s_sq=s_sq[i], n=n_vect[i], V=Vp, gamma=gamma, alpha=alpha, lower=0, upper=Inf)$value
		for(i in negm)
			Im[i] <- integrate(integrand, y_bar=y_bar[i], s_sq=s_sq[i], n=n_vect[i], V=Vm, gamma=gamma, alpha=alpha, lower=-Inf, upper=0)$value
	}
	if(!full)
		loglike <- -mean(log((w/(2*V))*gamma((nu+1)/2)*(Im+Ip)+(1-w)*gamma((nu-1)/2)*(sum_sq/2+1/alpha)^(-(nu-1)/2)))
	else
		loglike <- -mean(log((w/(2*V))*gamma((nu+1)/2)*(Im+Ip)+(1-w)*gamma((nu-1)/2)*(sum_sq/2+1/alpha)^(-(nu-1)/2))-log(alpha^gamma*gamma(gamma)))
	loglike
}

integrand <- function(mu, y_bar, s_sq, n, V, gamma, alpha)
{
(((n-1)*s_sq+n*(y_bar-mu)^2)/2+1/alpha+abs(mu)/V)^(-gamma-n/2-1)
}
