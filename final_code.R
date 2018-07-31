library(rstan); library(tidyverse)
setwd(dir = "Documents/STAT 615/Final Project/")
load(file = "commonGraphCounts.rda")
load(file = "d2l2_gc.rda")
load(file = "d1c3_gc.rda")

     
# ============================== PREP DATA ==============================
gc_common <- gc_common[,-c(2,3)]
# order of freq in data: X4_112 (.24), X6_112 (.11), X5_120 (.05), X3_224 (.05), X4_98 (.04), X8_112 (.04), X7_120 (.039), X2_192 (.029), X5_225 (.029), X4_120 (.029)

d <- data.frame(aggregate(.~IdentityName, data = gc_common, FUN = sum), 
                     N_w = rowSums(aggregate(.~IdentityName, data = gc_common, FUN = sum)[,-1]))

d_4 <- data.frame(aggregate(.~IdentityName, data = gc_common, FUN = sum)[,c("IdentityName", "X3_224", "X4_112", "X5_120", "X6_112")], 
                  N_w = rowSums(aggregate(.~IdentityName, data = gc_common, FUN = sum)[,c("X3_224", "X4_112", "X5_120", "X6_112")]))

d_10 <- data.frame(aggregate(.~IdentityName, data = gc_common, FUN = sum)[,c("IdentityName", "X2_192", "X3_224", "X4_112", "X4_120", "X4_98", "X5_120", "X5_225", "X6_112", "X7_120", "X8_112")], 
                  N_w = rowSums(aggregate(.~IdentityName, data = gc_common, FUN = sum)[,c("X2_192", "X3_224", "X4_112", "X4_120", "X4_98", "X5_120", "X5_225", "X6_112", "X7_120", "X8_112")]))

qd_4 <- d2l2_gc[,c("X3_224", "X4_112", "X5_120", "X6_112")]
qd_10 <- d1c3_gc[,c("X2_192", "X3_224", "X4_112", "X4_120", "X4_98", "X5_120", "X5_225", "X6_112", "X7_120", "X8_112")]


# ============================== RSTAN ====================================
##################################################################################
#   int N_w; //number of total graphs per writer

model = " 
data {
  int W; //number of writers (N)
  int G; //number of graph types considered (K)
  int Y[W,G]; //multinomial observations

  real<lower = 0> a;
  real<lower = 0> b; // gamma prior parameters
}

parameters { 
  vector<lower=0>[G] alpha;   
  simplex[G] theta[W]; 
} 

model { 
  for (g in 1:G) {
    alpha[g] ~ gamma(a, b);  // or some other prior on alpha 
}
  for (w in 1:W) { 
    theta[w] ~ dirichlet(alpha); 
    Y[w] ~ multinomial(theta[w]); 
  } 
} 
"

a = 2; b = 0.1 #set prior parameters
it = 10000 # number of iterations
dat = list(Y=d_10[,2:11], G = 10, W = nrow(d_10)) 
m = stan_model(model_code=model)
fit10 = sampling(m, data = dat, iter = it, chains = 1)

sim10 = rstan::extract(fit10, permuted=FALSE)
dimnames(sim10)

sim10 <- sim10[,1,-101]
sim10_d = data.frame(iters = 1:5000, sim10)
names(sim10_d)[2:101] = colnames(sim10)[1:100]

sim10_m = reshape2::melt(sim10_d,id=1, measure=2:101)
qplot(iters, value, data=sim10_m[(sim10_m$variable %in% c("alpha[1]","alpha[2]", "theta[1,1]", "theta[4,9]")),], alpha=I(.3), geom='line') + facet_grid(variable~.,scale='free_y') #chains

sim10_thetas <- droplevels(sim10_m[!(sim10_m$variable %in% c("alpha[1]","alpha[2]","alpha[3]","alpha[4]","alpha[5]","alpha[6]","alpha[7]","alpha[8]", "alpha[9]", "alpha[10]")),])
writer <- c(rep(rep(1:9, each = 5000), times = 10))
graph_type <- rep(c("X2_192", "X3_224", "X4_112", "X4_120", "X4_98", "X5_120", "X5_225", "X6_112", "X7_120", "X8_112"), each = 5000*9)
sim10_thetas <- cbind(writer, graph_type, sim10_thetas)

####### New writers ###########
qd_w1 <- data.frame(value = as.numeric(qd_10[1,]), graph_type = colnames(qd_10))
qd_w2 <- data.frame(value = as.numeric(qd_10[2,]), graph_type = colnames(qd_10))
qd_w3 <- data.frame(value = as.numeric(qd_10[3,]), graph_type = colnames(qd_10))
qd_w4 <- data.frame(value = as.numeric(qd_10[4,]), graph_type = colnames(qd_10))
qd_w5 <- data.frame(value = as.numeric(qd_10[5,]), graph_type = colnames(qd_10))
qd_w6 <- data.frame(value = as.numeric(qd_10[6,]), graph_type = colnames(qd_10))
qd_w7 <- data.frame(value = as.numeric(qd_10[7,]), graph_type = colnames(qd_10))
qd_w8 <- data.frame(value = as.numeric(qd_10[8,]), graph_type = colnames(qd_10))
qd_w9 <- data.frame(value = as.numeric(qd_10[9,]), graph_type = colnames(qd_10))


###### Posterior Probabilities ##########
posteriorProbs10 <- function(newCounts){
  thetas <- array(dim = c(5000, 10, 9))
  dmult <- data.frame()
  probs <- c()
  
  for(i in 1:9){
    for(j in 1:10){
      thetas[,j,i] <- sim10_d[,as.character(paste0("theta[", i, ",", j, "]"))]
    }
    
    for(k in 1:5000){
      dmult[k,i] <- dmultinom(x = newCounts, prob = thetas[k,,i])
    }
  }
  
  for(i in 1:9){
    probs[i] <- colSums(dmult)[i]/sum(colSums(dmult))
  }
  
  return(list(probs, dmult))
}

# 3, 5
temp <- posteriorProbs10(newCounts = qd_w5 $value)
temp2 <- temp[[2]]
round(temp1, 4)

colnames(temp2) <- 1:9
temp2 <- as.data.frame(t(apply(temp2, MARGIN = 1, function(x){x/sum(x)})))
temp2.m <- temp2 %>% gather(key = "writer", value = "value")


w5_var <- ggplot(data = temp2.m) + 
  geom_density(aes(x = value, colour = writer)) + 
  coord_cartesian(ylim = c(0, 5)) + 
  xlab("Probability") + 
  ggtitle(expression(paste("Density Curves for Estimates of P(", 'Y'^"*", "|", pi[w], ", Y", ")")))
w5_var

gridExtra::grid.arrange(w5_var)

######### PLOTS ###############
sim10_alphas <- droplevels(sim10_m[(sim10_m$variable %in% c("alpha[1]","alpha[2]","alpha[3]","alpha[4]","alpha[5]","alpha[6]","alpha[7]","alpha[8]", "alpha[9]", "alpha[10]")),])
#writer <- c(rep(rep(1:9, each = 5000), times = 10))
graph_type <- rep(c("X2_192", "X3_224", "X4_112", "X4_120", "X4_98", "X5_120", "X5_225", "X6_112", "X7_120", "X8_112"), each = 5000)
sim10_alphas <- cbind(graph_type, sim10_alphas)
as <- sim10_alphas  %>% group_by(graph_type) %>% summarise(quantile(value, probs = .025))

ggplot(sim10_alphas) + 
  geom_density(aes(x = value, color = as.factor(graph_type))) + 
  labs(color = expression(paste("Graph Type (", g, ")"))) + 
  xlab(expression(paste(alpha[g], " values"))) + 
  ggtitle(expression(paste("Density Plots of Posterior Samples of ", alpha[g], "'s")))


ggplot(sim10_thetas) + 
  geom_density(aes(x = value, color = as.factor(graph_type))) + 
  facet_wrap(~writer, nrow = 9, strip.position="right") +
  labs(color = expression(paste("Graph Type (", g, ")"))) + 
  xlab(expression(paste(pi[wg], " values"))) +
  ggtitle(expression(paste("Density Plots of Posterior Samples of ", pi[w], "vectors (Faceted by Writer and Graph Type)"))) +
  theme(strip.text = element_text(size=10))
  

ggplot(sim10_thetas) + 
  geom_density(aes(x = value, color = as.factor(graph_type))) + 
  facet_grid(writer ~ as.factor(graph_type), scale = 'free_x') + 
  geom_vline(data = qd_w3, aes(xintercept = value/sum(qd_w3$value), color = as.factor(graph_type)), linetype = 2, show.legend = T) + 
  ggtitle(expression(paste("Density Plots of Posterior Samples of ", pi[w,g], "'s (Faceted by Writer and Graph Type)")), subtitle = "Questioned Writing (Truth = Writer 3)") +
  labs(color = expression(paste("Graph Type (", g, ")"))) +
  xlab(expression(paste(pi[g], " values"))) + 
  scale_x_continuous(breaks = c(0.05, 0.15, 0.25, 0.35))

  
