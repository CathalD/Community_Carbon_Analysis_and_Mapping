lines <- readLines("docs/workshop.md",warn=FALSE)
in_code <- FALSE;code <- character();n <- 0L;figure <- NULL;checked <- character()
for(line in lines) {
  if(grepl("^!\\[",line)) figure <- sub(".*\\]\\(\\.\\./(.*)\\)","\\1",line)
  if(line=="```r") {in_code<-TRUE;code<-character();next}
  if(line=="```" && in_code) {
    parse(text=code);n<-n+1L;in_code<-FALSE
    for(cl in code) {
      if(grepl('^source\\("',cl)) {
        script <- sub('^source\\("([^" ]+)".*','\\1',cl)
        stopifnot(file.exists(script));parse(script)
      }
      if(grepl('^view_figure\\("',cl)) {
        path <- sub('^view_figure\\("([^"]+)".*','\\1',cl)
        stopifnot(identical(path,figure),file.exists(path));checked<-c(checked,path)
      }
    }
  } else if(in_code) code<-c(code,line)
}
stopifnot(!in_code,length(unique(checked))==20)
message("Workshop checks passed: ",n," parseable code blocks; 20 figure paths match the immediately preceding figures and maintained scripts.")
