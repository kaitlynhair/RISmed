# MAKING RNCBI PACKAGE
collapse <- function(...){paste(...,sep="",collapse="")}

EUtilsURL <- function(type="esearch",db="pubmed"){
	
	# CONSTRUCT ANY SERVICE TYPE AND DATABASE
	url <- 	"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/type.fcgi?db=DB&"
	
sub("(.*)(type)(.*)(DB)(.*)",collapse("\\1",type,"\\3",db,"\\5"),url)
}


EUtilsQuery <- function(query,type="esearch",db="pubmed",...){

	# CREATE URL WITH REQUEST
	PubMedURL <- EUtilsURL(type,db)
	Query <- gsub(" ","+",query)

	# PARTIAL MATCH OF POSSIBLE LIMITS
	OPTIONS <- c(
		"retstart",
		"retmax",
		"rettype",
		"field",
		"datetype",
		"reldate",
		"mindate",
		"maxdate"
	)
	
	ArgList <- list(...) # NEED FUNCTION PARTIAL MATCH HERE
	
	if(length(ArgList)==0){
		ArgList$retmax <- 1000 # DEFAULT 1000 RECORDS RETURNED
	}
	else{
		WhichArgs <- pmatch(names(ArgList),OPTIONS)	
		if(any(is.na(WhichArgs))||sapply(WhichArgs,length)>1)
			stop("Error in specified limits.")
		names(ArgList) <- OPTIONS[WhichArgs]
		if(all(names(ArgList)!="retmax"))
			ArgList$retmax <- 1000 # DEFAULT 1000 RECORDS RETURNED
	}
	
	ArgList$tool <- "RISmed"
	ArgList$email <- "s.a.kovalchik@gmail.com"
	
	# REPLACE RETMAX IF NOT USED
	ArgStr <- paste(names(ArgList),unlist(ArgList),sep="=")
	ArgStr <- paste(ArgStr,collapse="&")

paste(PubMedURL,"term=",Query,"&",ArgStr,sep="",collapse="")	
}

ParseTags <- function(lines){

	StripLines <- sapply(lines, function(x)gsub("> +",">",x)) # REPLACE WHITE SPACE
	
	Fields <- c("Count",
				"RetMax",
				"RetStart",
				"Id",
				"QueryTranslation")
			
	Patterns <- paste("(.*<",Fields,">)(.*)(<\\/",Fields,">).*",sep="")

	FieldIndex <- lapply(Patterns, function(pattern) grep(pattern, StripLines))
	FieldIndex[[1]] <- FieldIndex[[1]][1] # TAKE COUNT OF COMBINED QUERY

	Values <- lapply(1:length(Fields),
					function(i){
						result <- sapply(StripLines[FieldIndex[[i]]],
										function(x) sub(Patterns[i],"\\2",x),USE.NAMES=FALSE)
as.vector(result)
})


names(Values) <- Fields
Values$Count <- as.numeric(Values$Count)
Values$RetMax <- as.numeric(Values$RetMax)
Values$RetStart <- as.numeric(Values$RetStart)

Values
}


SplitIDs <- function(ids){

	if(length(ids)>200){
		group <- rep(1:ceiling(length(ids)/200),each=200)
		group <- group[1:length(ids)]
		split(ids,group)
	}
	else{
		list(ids)
	}

}


EUtilsGet <- function(x, type="efetch", db="pubmed"){

	if(class(x)[1]=="EUtilsSummary"){
		query <- x@querytranslation
		x <- x@PMID
	}
	else{
		query <- ""
	}
	
	IDList <- SplitIDs(x)
	Result <- lapply(IDList, EUtilsSubGet, type = type, db = db)
	Result <- unlist(Result, recursive=FALSE)

	if(type=="efetch"&db=="pubmed"){
		
		records <- Medline(Result, query)
	}
			
	authors<-vector()
	for(n in 1:length(records@Author)){
	  try(tmp<-paste0(records@Author[[n]]$LastName,", ",records@Author[[n]]$Initials), silent=TRUE)
	  try(tmp2<-paste(tmp, collapse="; "), silent=TRUE)
	  
	  if(length(tmp2)>0){
	    authors<-append(authors,tmp2)
	  }
	  
	  else{
	    authors<-append(authors,NA)
	  }
	}
	
	country <- vector()
	for(n in 1:length(records@Country)){
	  try(tmp<-paste0(records@Country[[n]]), silent=TRUE)
	  try(tmp2<-paste(tmp, collapse="; "), silent=TRUE)
	  
	  if(length(tmp2)>0){
	    country<-append(country,tmp2)
	  }
	  
	  else{
	    country<-append(country,NA)
	  }
	}
	
	keywords<-vector()
	for(n in 1:length(records@Keywords)){
	  try(tmp<-paste0(records@Keywords[[n]]), silent=TRUE)
	  try(tmp2<-paste(tmp, collapse="; "), silent=TRUE)
	  
	  if(length(tmp2)>0){
	    keywords<-append(keywords,tmp2)
	  }
	  
	  else{
	    keywords<-append(keywords,NA)
	  }
	}
	
	
	pubmed_data <- data.frame("author" = authors,
	                          "keywords" = keywords,
	                          "abstract" = AbstractText(records),
	                          "author_country" = country,
	                          "title" = ArticleTitle(records),
	                          "pages" = MedlinePgn(records),
	                          "issue" = Issue(records), 
	                          "volume" = Volume(records),
	                          "year" = YearPubmed(records),
	                          "pmid" = PMID(records),
	                          "doi" = DOI(records),
	                          "issn" = ISSN(records),
	                          "journal" = MedlineTA(records))
	
	
	pubmed_data <- pubmed_data %>%
	  mutate(url= paste0("https://www.ncbi.nlm.nih.gov/pubmed/",
	                     PMID(records))) %>%
	  mutate_if(is.factor, as.character)
	
	pubmed_data[pubmed_data == "NA" ] <- NA
	
	pubmed_data$uid <- paste0("pubmed-", pubmed_data$pmid)
	pubmed_data$source <- "pubmed"
	
}


EUtilsSubGet <- function(ids, type="efetch", db="pubmed"){

	FetchURL <- EUtilsURL(type,db=db)
	IDStr <- collapse("id=",paste(ids,collapse=","))
	EUtilsFetch <- collapse(FetchURL,IDStr)	
	
	res <- readLines(collapse(EUtilsFetch,"&retmode=xml"), warn = FALSE, encoding = "UTF-8")	
	res <- readLines(collapse(EUtilsFetch,"&retmode=xml"), warn = FALSE, encoding = "UTF-8")
	res <- xml2::read_xml(collapse(EUtilsFetch,"&retmode=xml"), encoding = "UTF-8")	
	res <- xml2::as_list(res)[[1]]
	
	if(db == "pubmed"){
		res <- lapply(res, ParseMedline) # Return list for each artcile, missing elements are NULL
	}

res
}




