library(pdftools)

# Lire le PDF page par page
pages <- pdf_text("/home/jules/Documents/phd/redaction/PANASENKAVA_Veranika.pdf")
pages %>% length()

writeLines(pages[6:104], "/home/jules/Documents/phd/redaction/panasenkava.txt")

# Lire le PDF page par page
pages <- pdf_text("/home/jules/Documents/phd/redaction/KIM_Artem.pdf")
pages %>% length()

writeLines(pages[2:80], "/home/jules/Documents/phd/redaction/kim.txt")
