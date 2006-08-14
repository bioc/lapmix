lap.volcanoplot <- function(res, highlight=0, ...)
{
	M <- res$M
	odds <- res$lap.lods
	plot(M, odds, xlab='average log fold change', ylab='log posterior odds', pch=16, cex=0.8, ...)
	if(highlight > 0)
	{
		table <- laptopTable(res, number=highlight)
		points(table$M, table$log_odds, pch=3, col='red')
	}
	invisible()
}
