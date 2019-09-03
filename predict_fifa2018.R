stats2018 <- read.csv("C:\\Users\\doyle\\Documents\\QLS_Projects\\FIFA 2018 Statistics.csv")

odds <- seq(1,nrow(stats2018),2)
evens <- seq(2,nrow(stats2018),2)

#Make new data frame
matches2018 <- data.frame("TeamA" = stats2018$Team[odds],
                          "TeamB" = stats2018$Team[evens],
                          "GoalsA" = stats2018$GoalsScored[odds],
                          "GoalsB" = stats2018$GoalsScored[evens],
                          "PSOGoalsA" = stats2018$PSOGoals[odds],
                          "PSOGoalsB" = stats2018$PSOGoals[evens])

matches2018$Winner <- ifelse( #If team A scored more goals in regular time or PSO
  matches2018$GoalsA > matches2018$GoalsB | matches2018$PSOGoalsA > matches2018$PSOGoalsB,
  as.character(matches2018$TeamA), #Then assign TeamA to Winner; otherwise
  ifelse(matches2018$GoalsA < matches2018$GoalsB | matches2018$PSOGoalsA < matches2018$PSOGoalsB,
         as.character(matches2018$TeamB), #if teamB scored more, assign it to Winner; otherwise
         "Tie")) #Make it a tie

matches2018$Winner <- as.factor(matches2018$Winner)


