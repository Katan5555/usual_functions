# https://medium.com/@ketanrd.009/how-to-extract-pdf-tables-in-r-e994c0fe4e28

library(tabulizer)
library(dplyr)


# Location of the pdf file
location <- 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4583086/pdf/2029.pdf'

# Extract the table at page 4
out <- extract_tables(location,
                            output = "data.frame",
                            pages = c(4),
                            guess = FALSE
                            
)
 
final <- out[[1]]
final <- final[3:58, ]

write.csv(final, "/Users/miyang/Downloads/pdf_table_extracted.csv") 
