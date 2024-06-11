library(crossdome)

# Open my file
mouse_peptides_h2db <- read.csv("mouse_mhc_peptides_h2db.txt", header=FALSE)
mouse_peptides_h2kb <- read.csv("mouse_peptides_h2kb.txt", header=FALSE)

default_database <- crossdome::hla_database

# Add each mouse_peptide to the database.
h2db <- data.frame()
for (i in mouse_peptides_h2db){
  new_row <- data.frame(peptide_sequence = i, hla_allele = "H2DB", mhcflurry_affinity="", mhcflurry_presentation_score="", 
  mhcflurry_presentation_percentile="", immunogenicity="", peptide_length=9, resource="", n_resource=1)
  h2db <- rbind(h2db, new_row)
}

h2kb <- data.frame()
for (i in mouse_peptides_h2db){
  new_row <- data.frame(peptide_sequence = i, hla_allele = "H2KB", mhcflurry_affinity="", mhcflurry_presentation_score="", 
  mhcflurry_presentation_percentile="", immunogenicity="", peptide_length=9, resource="", n_resource=1)
  h2kb <- rbind(h2kb, new_row)
}

database_h2db <- rbind(default_database, h2db)
database_h2db_h2kb <- rbind(database_h2db, h2kb)
  


viral_peptides_h2db <- read.csv("viral_peptides_H2Db.txt", header=FALSE)
viral_peptides_h2kb <- read.csv("viral_peptides_H2KB.txt", header=FALSE)

for (i in unlist(viral_peptides_h2db)){
  custom_background <- database_h2db_h2kb[database_h2db_h2kb$hla_allele == "H2DB", 'peptide_sequence']
  custom_object <- new("xrBackground", allele = "H2DB", peptides = custom_background, stats = list("off-target"="", "database" = length(custom_background)))
  print(i)
  result <- cross_compose(query=i, background=custom_object)
  write.csv(result@result, file=paste("H2DB/", i, ".csv", sep=""))
}

for (i in unlist(viral_peptides_h2kb)){
  custom_background <- database_h2db_h2kb[database_h2db_h2kb$hla_allele == "H2KB", 'peptide_sequence']
  custom_object <- new("xrBackground", allele = "H2KB", peptides = custom_background, stats = list("off-target"="", "database" = length(custom_background)))
  print(i)
  result <- cross_compose(query=i, background=custom_object)
  write.csv(result@result, file=paste("H2KB/", i, ".csv", sep=""))
}