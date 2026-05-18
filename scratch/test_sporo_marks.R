library(spatstat.data)

data(sporophores)
cat("Dataset Summary:\n")
print(sporophores)

cat("\nSpecies (Marks) Summary:\n")
print(summary(marks(sporophores)))
