stats2018 <- read.csv("FIFA_2018_Statistics.csv")

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


## R data processing for sane people
library(tidyverse)
## load data
train <- read_csv("HistoricWorldCupMatches.csv") %>% select(Home_Team_Name, Home_Team_Goals, Away_Team_Name, Away_Team_Goals) %>% drop_na()
## goals scored at home
home.goals <- train %>%  select(Home_Team_Name,  Home_Team_Goals) %>% rename(Goals = Home_Team_Goals, Team = Home_Team_Name)
## goals scored away
away.goals <- train %>%  select(Away_Team_Name,  Away_Team_Goals) %>% rename(Goals = Away_Team_Goals, Team = Away_Team_Name)

## calculate goal and game stats for every team, both at home and away
goals <- rbind(home.goals, away.goals) %>% group_by(Team) %>% summarise(TotalGoals = sum(Goals), Games = n()) %>% drop_na()

## Calculate Win/Loss outcomes. Draws are randomized since we want to apply logistic regresson
data <- train %>% select(Home_Team_Name, Away_Team_Name, Home_Team_Goals, Away_Team_Goals) %>%
  mutate(Outcome = ifelse(Away_Team_Goals > Home_Team_Goals, 0, ifelse(Home_Team_Goals > Away_Team_Goals, 1, sample(c(0,1), 1))))

## calculate the total games and total goals of each team in each game
data1 <- merge(data, goals, by.x="Home_Team_Name", by.y="Team")
data2 <-merge(data1, goals, by.x="Away_Team_Name", by.y="Team") %>%
  rename(HomeTotalGoals = TotalGoals.x, HomeTotalGames=Games.x, AwayTotalGoals=TotalGoals.y, AwayTotalGames=Games.y)

## generalized linear models
library(glmm)

## model fit
model <- glm(Outcome ~ HomeTotalGames + HomeTotalGoals + AwayTotalGames + AwayTotalGoals, family=binomial(link='logit'),data=data2)
## barf... they don't predict sauat
summary(model)

data2$Outcome <- as.factor(data2$Outcome)
randomForest(Outcome ~ HomeTotalGames + HomeTotalGoals + AwayTotalGames + AwayTotalGoals, data=data2)

## Confusion matrix:
## 0   1 class.error
## 0 60 114  0.65517241
## 1 65 613  0.09587021
