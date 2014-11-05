


#' @export
my.write.table <- function(x,file="",append=FALSE,quote=FALSE,dec=".",
                           row.names=FALSE,col.names=TRUE,width=12,digits=2,
                           justify="right",scientific=FALSE)
{   # Combines "format" and "write.table" to produce a prettier ascii table

    header = format(names(x),justify=justify,width=width)
    out = format(x,justify=justify,scientific=scientific,width=width,digits=digits)
    if (col.names) out = rbind(header,out)
    write.table(out,file=file,quote=quote,row.names=row.names,col.names=FALSE)
    return(head(out))
}
