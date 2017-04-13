items <- read.csv("MCTs_20239_Items.csv", TRUE, ",")
test <- items[items$UserName == "MCTs05",]
Protein <- mean(test$PROT)
FAT <- mean(test$TFAT)
Cals <- mean(test$KCAL)
Carbs <- mean(test$CARB)

## @knitr part1
Daily_Prot <- by( test$PROT, test$RecordDayNo, mean)
Daily_Fat <- by( test$PROT, test$RecordDayNo, mean)
Daily_Cals <- by( test$PROT, test$RecordDayNo, mean)
Daily_Carbs <- by( test$PROT, test$RecordDayNo, mean)

## @knitr part2
plot(Daily_Prot,type = "b", col = "red", xlab = "Day", ylab = "Prot Amount",main = "Protein Intake")

