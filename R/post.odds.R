post_odds <- function(w, V, beta, gamma, alpha, y_bar, s_sq, sum_sq, n_vect)
{
	G <- length(n_vect)
	nu <- n_vect+2*gamma+1
	Vp <- V*(1+beta)
	Vm <- V*(1-beta)
	Bp <- y_bar-1/(Vp*n_vect)
	Bm <- y_bar+1/(Vm*n_vect)
	Cp <- ((n_vect-1)*s_sq+2*y_bar/Vp-1/(n_vect*Vp^2)+2/alpha)/n_vect
	Cm <- ((n_vect-1)*s_sq-2*y_bar/Vm-1/(n_vect*Vm^2)+2/alpha)/n_vect
	Dp <- sqrt(abs(nu/Cp))
	Dm <- sqrt(abs(nu/Cm))
	negp <- which(Cp<0)
	negm <- which(Cm<0)
	E <- sqrt(nu*pi)*(gamma(nu/2)/gamma((nu+1)/2))*(n_vect/2)^(-(nu+1)/2)
	Ip <- E*(Cp^(-(nu+1)/2)*pt(Dp*Bp, df=nu)/Dp)
	Im <- E*(Cm^(-(nu+1)/2)*pt(-Dm*Bm, df=nu)/Dm)
	for(i in negp)
		Ip[i] <- integrate(integrand, y_bar=y_bar[i], s_sq=s_sq[i], n=n_vect[i], V=Vp, gamma=gamma, alpha=alpha, lower=0, upper=Inf)$value
	for(i in negm)
		Im[i] <- integrate(integrand, y_bar=y_bar[i], s_sq=s_sq[i], n=n_vect[i], V=Vm, gamma=gamma, alpha=alpha, lower=-Inf, upper=0)$value
	odds <- w*((nu-1)/2)*(Ip+Im)*(sum_sq/2+1/alpha)^((nu-1)/2)/((1-w)*2*V)
	prob <- odds/(1+odds)
	l_threshold <- Im*(Im+Ip+2*V*(2*(1-w)/(w*(nu-1)))*(sum_sq/2+1/alpha)^(-(nu-1)/2))^(-1)
	h_threshold <- l_threshold+1-prob
	med.number <- sum(!((l_threshold<1/2)*(h_threshold>1/2)))
	list(odds=odds, prob=prob, med.number=med.number)
}
