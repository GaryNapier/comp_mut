

n <- 100
maj <- 0.9
samp <- c(0, 1)
probs <- c(maj, 1-maj)
probs_rev <- rev(probs)
dst_zero <- sample(samp, n, replace = T, prob = probs)
dst_one <- sample(samp, n, replace = T, prob = probs_rev)
mut_zero <- sample(samp, n, replace = T, prob = probs)
mut_one <- sample(samp, n, replace = T, prob = probs_rev)

x <- data.frame(dst = c(dst_zero, dst_one), has_mut = c(mut_zero, mut_one)) 

cor(x)
# dst has_mut
# dst     1.000   0.722
# has_mut 0.722   1.000

model <- glm(formula = dst ~ has_mut, data = x, family = binomial(link = "logit"))

summary(model)

cols <- ifelse(x$dst == 1, "red", "blue")
jit <- 0.5
plot(jitter(x$dst, jit), jitter(x$has_mut, jit), col = cols)
x_dst <- seq(0, 1, 0.01)
y_has_mut <- predict(model, list(has_mut = x_dst),type="response")
lines(x_dst, y_has_mut)

# ----------------

lineages <- c("lin1", "lin2", "lin3", "lin4")

lin_zero <- sample(lineages, n, replace = T)
lin_one <- sample(lineages, n, replace = T, prob = c(maj, rep((1-maj)/3, 3)))

x$lin <- c(lin_zero, lin_one)

x

model_2 <- glm(formula = dst ~ lin, data = x, family = binomial(link = "logit"))

summary(model_2)

model_3 <- glm(formula = dst ~ lin + has_mut, data = x,  family = binomial(link = "logit"))

summary(model_3)

anova(model_3, model, test='Chisq')


cols_df <- data.frame(lin = sort(unique(x$lin)), col = c('red', 'green', 'blue', 'orange'))
cols_df <- cols_df[match(x[,'lin'], cols_df$lin), ]
plot(jitter(x$dst, jit), jitter(x$has_mut, jit), col = cols_df$col)

# y_has_mut <- predict(model_3, list(has_mut = x_dst), type="response")
# y_has_mut <- predict(model_3, type = "response")

# lines(x_dst, c(y_has_mut[1], y_has_mut))


# library(ggplot2)
# 
# qplot(x = jitter(has_mut), y = jitter(dst), color = lin, data = x) 
#   # stat_smooth(method = "lm", se = FALSE, fullrange = TRUE)



# ------

# https://www.r-bloggers.com/2021/05/how-to-generate-correlated-data-in-r/

library(MASS)

set.seed(5)
# create the variance covariance matrix
sigma <- rbind(c(1, 0.8, 0.7), c(0.8, 1, 0.9), c(0.7,0.9,1))
# create the mean vector
mu <- c(10, 5, 2) 
# generate the multivariate normal distribution
df <- as.data.frame(MASS::mvrnorm(n = 100, mu = mu, Sigma = sigma))

model <- lm(V1 ~ V2, data = df)
summary(model)

model_2 <- lm(V1 ~ V3, data = df) 
summary(model_2)

model_3 <- lm(V1 ~ V2 + V3, data = df)
summary(model_3)


# ------

library(tidyverse)
library(GGally)

df <- df %>% mutate(MyBinary = ifelse(V1 > median(V1), 1 ,0))

df <- df %>% mutate(MyNoisyBinary = ifelse(V1 > median(V1), 
                                           sample(c(0,1), 
                                                  n(), 
                                                  replace = TRUE, 
                                                  p = c(0.25, 0.75)),
                                           sample(c(0,1), n(),  replace = TRUE, p=c(0.75, 0.25))))

df <- df %>% mutate(AgeGroup= case_when(V1 < quantile(V1, 0.25) ~ "Group 1",
                                    V1 < quantile(V1, 0.5) ~ "Group 2",
                                    V1 < quantile(V1, 0.75) ~ "Group 3",
                                    TRUE ~ "Group 4"))

ggpairs(df)










