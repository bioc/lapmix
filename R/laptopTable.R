laptopTable <- function(res, number=res$med.number, sort.by='L')
{
	lap.lods <- res$lap.lods
	prob <- res$prob
	nondiff_prob <- 1-prob
	Im <- res$Im
	Ip <- res$Ip
	V <- res$estimates$V
	w <- res$estimates$w
	gamma <- res$estimates$gamma
	alpha <- res$estimates$alpha
	if(number <= 0)
		stop('Number of displayed genes must be a positive integer!')
	if(sort.by == 'L')
	{
		ix <- sort(lap.lods, decreasing=TRUE, index.return=TRUE)$ix
		diff.genes <- ix[1:number]
	}
	else if(sort.by == 'M')
	{
		ix <- sort(res$M, decreasing=TRUE, index.return=TRUE)$ix
		diff.genes <- ix[1:number]	
	}
	else
		stop('invalid value given to argument sort.by')
	table <- data.frame(gene=diff.genes, M=res$M[diff.genes], log.odds=lap.lods[diff.genes], row.names=NULL)
	table
}

