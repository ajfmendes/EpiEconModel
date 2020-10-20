### Author: AJ Mendes.
### Date: 30th July 2020.
### Title: Socially versus privately optimal control of livestock diseases: a case for integration of epidemiology and economics.

## Load libraries:
packages <- c("deSolve", "ggplot2", "extrafont", "tidyr", "DataCombine")
installed <- packages %in% installed.packages()
if(length(packages[!installed]) > 0) install.packages(packages[!installed])
lapply(packages, require, character.only=TRUE)
loadfonts(device = "win")
rm(packages,installed)


#######################################################################################
########     Scenario 1: Demonstrating the impact of farmer heterogeneity      ########
#######################################################################################

# set up populations
nfarms <- 1 # number of farms
n.per.farm <- 250 # number of animals per farm
prop.infected <- 0.1 # proportion of infected animals
N <- nfarms*n.per.farm # total animal population in the model
I_0<-rep(prop.infected*n.per.farm,times = nfarms) # initial infecteds per farm
S_0<-N-I_0 # the remaining are susceptible
initial.pop <- c(S=S_0,I=I_0) # vectorize initial populations

# Set up parameters
pi <- 0.2 # birth rate
miu <- pi # death rate = birth rate

# set up transmission rates
R0 <- n.per.farm/(n.per.farm*(1-prop.infected)) # calculate R0 of the disease knowing that it is in equilibrium (Re=1.0)
bw <- 0 # ratio between to within farms transmission rate is 0; we want on-farm dynamics only for now
within.beta <- ((R0*miu)/n.per.farm)*(1-bw) # within farm transmission rate
between.beta <- ((R0*miu)/n.per.farm)*(bw) # between farm transmission rate
betas <- matrix(rep(between.beta), nrow=nfarms,ncol=nfarms) # create matrix of betas (we only need a matrix in the second scenario but we use them here for consistency)
diag(betas) <- within.beta # betas for intra-farm (diagonals)
prms <- c(pi=pi, miu=miu, betas=betas) # vectorize prms

# Simulation times
start.time <- 0
end.time <- 30 # 25 time units of simulation after burning five time units (before awareness campaign).
step <- 1
times <- seq(start.time,end.time,step)

# Run simulation
final.df <- data.frame(matrix(,0,5)) # create empty data frame that will contain simulation results
xs <- 1:3 # x for each type of farmer
for (x in xs) {
  ode_f <- function(t, pop, prms){  # create ode function
    ns <- length(initial.pop)/nfarms # number of statuses (S and I)
    nfarms <- length(pop)/ns
    S    <- as.matrix(pop[1:nfarms]) # S matrix
    I    <- as.matrix(pop[(nfarms+1):(2*nfarms)]) # I matrix
    
    with(as.list(prms),{
      prev <- I/(S+I) # prevalence at each time step
      fh <- function(i) { (1/(1+(exp(-100*(i-0.025)))))} # fh: Highly responsive farmer
      fm <- function(i) { (1/(1+(exp(-100*(i-0.1)))))}   # fm: Moderately responsive farmer
      fs <- function(i) { (1/(1+(exp(-100*(i-0.175)))))} # fs: Slightly responsive farmer
      f3 <- list(fh,fm,fs) # list farmers' functions
      
      if (t < 5) {  # burn the first 5 time steps
        betas <- matrix(rep(between.beta), nrow=nfarms,ncol=nfarms)
        diag(betas) <- within.beta
      } else {      # after the first 5 time steps farmers take action by reducing the transmission rate based on their level of responsiveness
        betas <- matrix(rep(between.beta), nrow=nfarms,ncol=nfarms)
        diag(betas) <- within.beta-(within.beta*f3[[x]](prev))
      }
      
      # Define our differential equations
      dS_dt <-  pi*(S+I) - S*betas%*%(I) - miu*S
      dI_dt <-  S*betas%*%(I)  - miu*I
      
      # solutions
      output<-c(dS_dt,dI_dt)
      list(output,prev=prev) # population numbers and prevalence
    })
  }
  df <- as.data.frame(ode(initial.pop, times, ode_f, prms)) # get solutions in one data frame
  df$x <- x # keep record of the farmer type (1,2,3)
  final.df <- rbind(final.df,df) # combine data frames with output for each farmer type
}

# Plot farmers' functions
# Fig. 1. Panel A: control effort applied (percentage reduction in transmission rate - beta - within the farm in relation to the initial beta in endemic equilibrium) as a function of disease prevalence for three different farmer types. The vertical dash-dotted line indicates initial prevalence. 
x<-c(0,1.999);ymin2<-c(0,0);ymax2<-c(100,100);rib2 <- data.frame(x,ymin2,ymax2) # df for shaded area
fh <- function(i) { (1/(1+(exp(-100*(i-0.025)))))} # fh: Highly responsive farmer (as in ode_f)
fm <- function(i) { (1/(1+(exp(-100*(i-0.1)))))}   # fm: Moderately responsive farmer (as in ode_f)
fs <- function(i) { (1/(1+(exp(-100*(i-0.175)))))} # fs: Slightly responsive farmer (as in ode_f)
x <- seq(0,1,0.001);a <- data.frame(x=x,y=fh(x),z="a");b <- data.frame(x=x,y=fm(x),z="b");d <- data.frame(x=x,y=fs(x),z="d");data <- rbind(a,b,d);data$y <- data$y*100 # create dummy prevalences and solve farmers' functions
fig1A <- ggplot(data = data, mapping = aes(x=x))+ 
  geom_ribbon(data=rib2, aes(ymin=as.numeric(ymin2),ymax=as.numeric(ymax2)),alpha = 0.15) + 
  geom_vline(xintercept = 10, linetype=4, color = "grey25", size=0.5) +
  geom_line(data = data,mapping = aes(x=x*100, y=y, color=z, lty=factor(z)), size=0.75) +
  scale_colour_manual(name='Farmer type   ',values=c("#F8766D", "#7CAE00", "#00BFC4"),labels=c('Highly responsive','Moderately responsive','Slightly responsive'))+ scale_linetype_manual(name='Farmer type   ',values=c("solid", "dotted", "dashed"),labels=c('Highly responsive','Moderately responsive','Slightly responsive'))  + 
  ggtitle("(A)") + 
  xlab("Prevalence (%)") + 
  ylab(bquote(paste("Reduction in within-farm    ",beta, " (%)"))) + 
  xlim(0,25) + 
  theme_bw() + 
  theme(legend.position="bottom", legend.box = "horizontal", text=element_text(size=11,  family="Times New Roman"))
fig1A

# Plot simulation results
# Fig. 1. Panel B: disease prevalence over time within the three farms.
time<-c(0,30);ymin1<-c(0,0);ymax1<-c(1.999,1.999);rib1 <- data.frame(time,ymin1,ymax1) # df for shaded area
fig1B <- ggplot(data = final.df, mapping = aes(x = time-5)) +
  geom_ribbon(data=rib1, aes(ymin=as.numeric(ymin1),ymax=as.numeric(ymax1)),alpha = 0.15) +
  geom_line(data = final.df, mapping = aes(x=time-5, y=prev*100, color=factor(x), lty=factor(x)), size=0.75) +
  scale_colour_manual(name='Farmer type   ',values=c("#F8766D", "#7CAE00", "#00BFC4"),labels=c('Highly responsive','Moderately responsive','Slightly responsive')) +
  scale_linetype_manual(name='Farmer type   ',values=c("solid", "dotted", "dashed"),labels=c('Highly responsive','Moderately responsive','Slightly responsive')) +
  scale_x_continuous(breaks=seq(-5,25,5)) +
  ggtitle("(B)") +
  xlab("Time") +
  ylab("Prevalence (%)") +
  ylim(0,NA) +
  theme_bw() +
  theme(legend.position="bottom", legend.box = "horizontal", text=element_text(size=11,  family="Times New Roman"))
fig1B
rm(list=ls(all=TRUE)) # clear workspace


#######################################################################################
############           Scenario 2: Demonstrating externalities          ###############
#######################################################################################

# set up populations
nfarms <- 2 # number of farms
n.per.farm <- 250 # number of animals per farm
prop.infected <- 0.1 # proportion of infected animals
total.pop <- nfarms*n.per.farm # total animal population in the model
p <- rep(1/nfarms,times = nfarms) # vector with the proportion of farms
N <- total.pop*p # vector with the number of animals per farm
I_0 <- rep(prop.infected*n.per.farm,times = nfarms) # initial infecteds per farm
S_0 <- N-I_0 # the remaining are susceptible
initial.pop <- c(S=S_0,I=I_0) # vectorize initial populations

# Set up parameters
pi <- 0.2 # birth rate
miu <- pi # death rate = birth rate

# Simulation times
start.time <- 0
end.time <- 30 # 25 time units of simulation after burning five time units (before awareness campaign).
step <- 1
times <- seq(start.time,end.time,step)  

# Run simulation
xs <- c(1,2) # 2 loops, one without transmission between farms and one with
ci <- 50 # annual cost on one infected animal in monetary units
data <- data.frame(matrix(,0,6)) # create empty data frame that will contain simulation results
for (x in xs) {
  # set up transmission rates
  R0 <- n.per.farm/(n.per.farm*(1-prop.infected)) # calculate R0 of the disease knowing that it is in equilibrium (Re=1.0)
  if (x==1){ # only in the second loop we allow disease transmission between farms
    bw <- 0 } else {
      bw <- 0.1 } # the transmission rate between farms is 0% and 10% of the initial within-farm transmission rate.
  within.beta <- ((R0*miu)/n.per.farm)*(1-bw) # within farm transmission rate
  between.beta <- ((R0*miu)/n.per.farm)*(bw) # between farm transmission rate
  betas <- matrix(rep(between.beta), nrow=nfarms,ncol=nfarms) # create matrix of betas
  diag(betas) <- within.beta # betas for intra-farm (diagonals)
  prms <- c(pi=pi, miu=miu, betas=betas) # vectorize parameters
  
  
  ode_f <- function(t, pop, prms){  # create ode function for each type of farmer
    ns <- length(initial.pop)/nfarms # number of statuses (S and I)
    nfarms = length(pop)/ns
    S    <- as.matrix(pop[1:nfarms]) # S matrix
    I    <- as.matrix(pop[(nfarms+1):(2*nfarms)]) # I matrix
    
    with(as.list(prms),{
      prev <- I/(S+I) # prevalence at each time step
      fh <- function(i) { (1/(1+(exp(-100*(i-0.025)))))} # fh: Highly responsive farmer
      fb <- function(i) { (1/(1+(exp(-100*(i-0.175)))))} # fs: Slightly responsive farmer
      
      if (t < 5) { # burn the first 5 time steps
        betas <- matrix(rep(between.beta), nrow=nfarms,ncol=nfarms)
        diag(betas) <- within.beta
      } else {  # after the first 5 time steps farmers take action by reducing the transmission rate in each farm
        betas[1,1] <- within.beta-(within.beta*fh(prev[1]))
        betas[2,2] <- within.beta-(within.beta*fb(prev[2]))
        betas[1,2] <- between.beta # Actions taken by each farmer do not influence the between.beta
        betas[2,1] <- between.beta # Actions taken by each farmer do not influence the between.beta
      }
      
      # Define our differential equations
      dS_dt <-  pi*(S+I) - S*betas%*%(I) - miu*S
      dI_dt <-  S*betas%*%(I)  - miu*I
      
      # solutions
      output<-c(dS_dt,dI_dt)
      list(output,prev1=prev[1],prev2=prev[2]) # population numbers and prevalence over time
    })
  }
  df <- as.data.frame(ode(initial.pop, times, ode_f, prms)) # get solutions in one data frame
  df$loop <- x # keep record of the loop
  data <- rbind(data,df) # combine data frames
}

# Plot simulation results
# Prepare ribbons
rib1 <- cbind(subset(data,  loop == 1,select = c(time,I1)),subset(data,  loop == 2,select = c(I1)))
names(rib1)[2] <- "loop1"
names(rib1)[3] <- "loop2"
rib1$farm <- "A" 
rib2 <- cbind(subset(data,  loop == 1,select = c(time,I2)),subset(data,  loop == 2,select = c(I2)))
names(rib2)[2] <- "loop1"
names(rib2)[3] <- "loop2"
rib2$farm <- "B"
rib <- rbind(rib1,rib2)

# Prepare dataset for lines
data <- data[-c(2:3)]
data <- data %>% gather(farm, I, I1:I2)
data$lty <- 0
data$lty[data$loop==1&data$farm=="I1"] <- 1
data$lty[data$loop==2&data$farm=="I1"] <- 2
data$lty[data$loop==1&data$farm=="I2"] <- 3
data$lty[data$loop==2&data$farm=="I2"] <- 4
rib$lty<- 1
rib2<-rib
rib2$lty <- 2
rib<-rbind(rib,rib2)
rib$lty2[rib$farm=="A"&rib$lty=="1"] <- 1
rib$lty2[rib$farm=="A"&rib$lty=="2"] <- 2
rib$lty2[rib$farm=="B"&rib$lty=="1"] <- 3
rib$lty2[rib$farm=="B"&rib$lty=="2"] <- 4
rib$lty2 <- factor(rib$lty2, levels=c("1", "2", "3","4"), labels=c("1","3", "2", "4"))
data$lty <- factor(data$lty, levels=c("1", "3", "2","4"), labels=c("1","3", "2", "4"))
rib.s <- data.frame(time=c(0,30),ymin1=c(0,0),ymax1=c(1.999,1.999),lty=c("I1","I1"))

# Plot ribbons and lines
fig2 <- ggplot(data = data, mapping = aes(x = time-5, group = factor(lty))) + 
  geom_ribbon(data=rib.s, aes(ymin=as.numeric(ymin1),ymax=as.numeric(ymax1)),alpha = 0.15) +
  geom_ribbon(data=rib, aes(ymin=as.numeric(loop1/250*100),ymax=as.numeric(loop2/250*100),group=factor(lty2),fill=factor(lty2)),alpha = 0.15) + 
  scale_fill_manual(name='Farmer type (scenario)',values=c("#F8766D", "#F8766D", "#00BFC4", "#00BFC4"),labels=c(bquote(paste("Highly responsive (b",beta==0,")")),bquote(paste("Slightly responsive (b",beta==0,")")),bquote(paste("Highly responsive (b",beta>0,")")),bquote(paste("Slightly responsive (b",beta>0,")")))) +
  guides(fill=FALSE) +
  geom_line(data = data, mapping = aes(x=time-5, y=I/250*100, color=factor(lty), lty=factor(lty)), size=0.75) +
  scale_colour_manual(name='Farmer type (scenario)',values=c("#F8766D", "#00BFC4", "#F8766D", "#00BFC4"),labels=c(bquote(paste("Highly responsive (b",beta==0,")")),bquote(paste("Slightly responsive (b",beta==0,")")),bquote(paste("Highly responsive (b",beta>0,")")),bquote(paste("Slightly responsive (b",beta>0,")")))) +
  scale_linetype_manual(name='Farmer type (scenario)',values=c(1,3,2,4),labels=c(bquote(paste("Highly responsive (b",beta==0,")")),bquote(paste("Slightly responsive (b",beta==0,")")),bquote(paste("Highly responsive (b",beta>0,")")),bquote(paste("Slightly responsive (b",beta>0,")")))) +
  scale_y_continuous(sec.axis = sec_axis(~.*125, name = bquote(paste(italic("Ex post  "), "cost of disease")),breaks=c(0,313,781,1250),labels=c("0", "Low","Medium","High")),limits =c(0,I_0[1]/250*100)) +
  scale_x_continuous(breaks=seq(-5,25,5))+
  xlab("Time") +
  theme_bw() +
  ylab("Prevalence (%)")

# Add arrows to indicate differences once externalities are taken into account
t <- rep(data[26,1],4)
loss <- ci*c(data[57,6],data[26,6],data[119,6],data[88,6])
arrow1 <- data.frame(t=t,loss=loss)
arrow1$lty <- c("I1","I1","I2","I2")
t <- rep(data[21,1],4)
loss <- ci*c(data[52,6],data[21,6],data[114,6],data[83,6])
arrow2 <- data.frame(t=t,loss=loss)
arrow2$lty <- c("I1","I1","I2","I2")
t <- rep(data[31,1],4)
loss <- ci*c(data[62,6],data[31,6],data[124,6],data[93,6])
arrow3 <- data.frame(t=t,loss=loss)
arrow3$lty <- c("I1","I1","I2","I2")
fig2 <- fig2 + geom_point(data=arrow1,aes(x=t-5,y=loss/50/250*100,group = factor(lty)), size=NA) + geom_line(data=arrow1,aes(x=t-5,y=loss/50/250*100), size=0.25, lty="solid", arrow = arrow(length=unit(0.1,"cm"), ends="first", type = "open")) +
  geom_line(data=arrow2,aes(x=t-5,y=loss/50/250*100), size=0.25, lty="solid", arrow = arrow(length=unit(0.1,"cm"), ends="first", type = "open")) +
  geom_line(data=arrow3,aes(x=t-5,y=loss/50/250*100), size=0.25, lty="solid", arrow = arrow(length=unit(0.1,"cm"), ends="first", type = "open")) +
  theme(text=element_text(size=11,  family="Times New Roman"))
fig2
rm(list=ls(all=TRUE)) # clear workspace


#######################################################################################
############          Scenario 3: Privately optimal behavior           ################
#######################################################################################

opt_finder <- function(t, v, params){  # create ode function
  nstatus <- length(initial.pop)/nfarms 
  nfarms = length(v)/nstatus
  
  # create a matrix for each epi status
  S    <- as.matrix(v[1:nfarms])
  I    <- as.matrix(v[(nfarms+1):(2*nfarms)])
  
  with(as.list(params),{
    # N is the total cattle population in the village - a vector of length nfarm
    prev <- I/(S+I)
    
    opt.red <- k
    
    exante <- (exa.unit*100*opt.red)/((1+disc.rate)^t) # calculate the discounted cost of control per time step
    expost <- (I*exp.unit)/((1+disc.rate)^t) # calculate the discounted ex post cost of disease
    
    betas <- matrix(rep(between.beta), nrow=nfarms,ncol=nfarms)
    diag(betas) <- within.beta-(within.beta*opt.red)
    
    # Define differential equations
    dS_dt <-  pi*(S+I) - S*betas%*%(I) - miu*S
    dI_dt <-  S*betas%*%(I)  - miu*I
    
    out<-c(dS_dt,dI_dt)
    list(out,prev=prev,opt.red=opt.red,exante=exante,expost=expost) # population numbers, prevalence over time and control effort over time
  })
}

### Run simulation looping over different values of unit cost of control
# set up populations
nfarms <- 1 # number of farms
n.per.farm <- 250 # number of animals per farm
prop.infected <- 0.1 # proportion of infected animals
total.pop <- nfarms*n.per.farm # total animal population in the model
p <- rep(1/nfarms,times = nfarms) # vector with the proportion of farms
N <- total.pop*p # vector with the number of animals per farm
I_0<-rep(prop.infected*n.per.farm,times = nfarms) # initial infecteds per farm
S_0<-N-I_0 # the remaining are susceptibles
initial.pop <- c(S=S_0,I=I_0) # vectorise initial populations

# Set up parameters
pi <- 0.2 # birth rate
miu <- pi # death rate = birth rate

# set up transmission rates
R0 <- n.per.farm/(n.per.farm*(1-prop.infected)) # calcaulate R0 of the disease knowing that it is in equilibrium (Re=1.0)
bw <- 0 # ratio between to within farms trasmission rate is 0; we want on farm dynamics only for now.
within.beta <- ((R0*miu)/n.per.farm)*(1-bw) # within farm transmission rate
between.beta <- ((R0*miu)/n.per.farm)*(bw) # between farm transmission rate
betas <- matrix(rep(between.beta), nrow=nfarms,ncol=nfarms) # create matrices of betas
diag(betas) <- within.beta # betas for intra-farm (diagonals)
params <- c(pi=pi, miu=miu, betas=betas) # vectorise params

## Time horizon for optimisation review
start.time <- 0
end.time <- 1 # the planning horizon is 1 timestep here
big.step <- 1
times <- seq(start.time,end.time,big.step)  
n.simul <- 30 # number of time steps/horizons to loop over

## Econ parameters
eus <- seq(0.0,1.25,0.05) # plausible unit cost of control (price of reducing beta by 1%)
disc.rate<- 0.0
exp.unit <- 50 # cost of one infected animal

## Set up data frames needed
k<-0 # to simulate populations before awareness campaign - no control action
exa.unit <- 0 # to simulate populations before awareness campaign - exa.unit could be any
before.camp <- as.data.frame(ode(initial.pop, seq(0,5,1), opt_finder, params)) 
data <- data.frame(matrix(,0,8))
for (exa.unit in eus) { # looping over different values of unit cost of control
  initial.pop <- c(S=S_0,I=I_0)  # update initial populations
  final.populations <- data.frame(matrix(,0,7))
  ops <- seq(0,1.0,0.0001)
  ### Run simulation in which farmers review their control policy every times step to achieve the economic optimum   
  for (x in seq(0,n.simul,end.time)) {
    list <- list()
    sums <- c()
    for (k in ops) {
      list[[which(ops==k)]] <- as.data.frame(ode(initial.pop, times, opt_finder, params)) # list of solutions under alternative control efforts ops
      
      # calculate the costs of control and ex post cost of disease in between every two time steps
      a<-list[[which(ops==k)]]$expost
      b<-list[[which(ops==k)]]$exante
      r<-c()
      s<-c()
      for (q in 1:length(a)){
        r[q]<-mean(c(a[q],a[q+1]))
        s[q]<-mean(c(b[q],b[q+1]))
      }
      r<-r[!is.na(r)] # delete last NA
      s<-s[!is.na(s)] # delete last NA
      sums[which(ops==k)] <- sum(sum(s),sum(r)) # which ops yields the lowest total costs of disease?
    }
    
    final.populations <- rbind(final.populations, as.data.frame(list[which(sums==min(sums))])[2:length(times),]) # output final populations under optimised control effort
    
    initial.pop <- c(S=final.populations[x+end.time,2],I=final.populations[x+end.time,3]) # update initial populations to retart where the time horizon finished
    
  }
  final.populations<- rbind(before.camp,final.populations[1:(n.simul-5),])
  final.populations$time <- seq(0,length(final.populations[,2])-1,1) # correct time steps 
  out_exa <- final.populations # to export with each exa unit
  out_exa$exa.unit <- exa.unit
  data <- rbind(data,out_exa)
}

# Plot simulation results 
time<-c(0,n.simul); ymin1<-c(0,0); ymax1<-c(1.999,1.999); rib1 <- data.frame(time,ymin1,ymax1) # shaded area
# Plot prevalence over time for each unit cost of control
cols <-c("#004FE2","#19DB00","#FFD000","#FC4700")
leg.values <- c(0.0,0.40,0.85,1.25) 
colfunc<-colorRampPalette(cols)
hl <- data.frame(x=c(-5,25),y=c(0.4,0.4))
fig3 <- ggplot(data = data, mapping = aes(x = time-5)) + 
  geom_ribbon(data=rib1, aes(ymin=as.numeric(ymin1),ymax=as.numeric(ymax1)),alpha = 0.15) +
  geom_path(data = hl, mapping = aes(x=x, y=y, fill=NULL), colour="Black",size=0.5, linetype = "dotted") + 
  geom_line(data = data, mapping = aes(x=time-5, y=prev*100, color=factor(exa.unit)), size=0.5) + 
  scale_color_manual(name="Price per control unit (range): ",values=colfunc(length(eus)),breaks=leg.values,labels=c("Free","Low","Medium","High")) + 
  scale_x_continuous(breaks = seq(-5,25,5)) +
  xlab("Time") + 
  theme_bw() + 
  ylab("Prevalence (%)") +
  theme(text=element_text(size=11,  family="Times New Roman"))
fig3

# now carry out a sensitivity analysis with different costs per infected animals as well: this may take several minutes (Fig. 4)
final.tiles <-  data.frame(matrix(,0,3)) # final.tiles for geom_tile
inc.exp <- 3 # increments in cost per infected animal
inc.eus <- 0.05 # increments in price per control unit
exps <- seq(10,85,inc.exp)
for (exp.unit in exps){
  ### Run simulation looping over different values of unit cost of control
  ## Econ parameters
  eus <- seq(0.0,1.25,inc.eus) # plausible unit cost of control (price of reducing beta by 1%)
  disc.rate<- 0.0 # farmer reviews optimal effort at each time step, so discount rate is 0%.
  ## Set up data frames needed
  k<-0 # to simulate populations before awareness campaign - no control action
  exa.unit <- 0 # to simulate populations before awareness campaign - exa.unit could be any
  data <- data.frame(matrix(,0,8))
  for (exa.unit in eus) { # looping over different values of unit cost of control
    initial.pop <- c(S=S_0,I=I_0)  # update initial populations
    final.populations <- data.frame(matrix(,0,7))
    ops <- seq(0,1.0,0.01)
    ### Run simulation in which farmers review their control policy every times step to achieve the economic optimum   
    for (x in seq(0,n.simul,end.time)) {
      list <- list()
      sums <- c()
      for (k in ops) {
        list[[which(ops==k)]] <- as.data.frame(ode(initial.pop, times, opt_finder, params)) # list of solutions under alternative control efforts ops
        
        # calculate the costs of control and ex post cost of disease in between every two time steps
        a<-list[[which(ops==k)]]$expost
        b<-list[[which(ops==k)]]$exante
        r<-c()
        s<-c()
        for (q in 1:length(a)){
          r[q]<-mean(c(a[q],a[q+1]))
          s[q]<-mean(c(b[q],b[q+1]))
        }
        r<-r[!is.na(r)] # delete last NA
        s<-s[!is.na(s)] # delete last NA
        sums[which(ops==k)] <- sum(sum(s),sum(r)) # which ops yields the lowest total costs of disease?
      }
      
      final.populations <- rbind(final.populations, as.data.frame(list[which(sums==min(sums))])[2:length(times),]) # output final populations under optimised control effort
      
      initial.pop <- c(S=final.populations[x+end.time,2],I=final.populations[x+end.time,3]) # update initial populations to retart where the time horizon finished
      
    }
    final.populations<- rbind(before.camp,final.populations[1:(n.simul-5),])
    final.populations$time <- seq(0,length(final.populations[,2])-1,1) # correct time steps 
    out_exa <- final.populations # to export with each exa unit
    out_exa$exa.unit <- exa.unit
    data <- rbind(data,out_exa)
  }
  prev.tiles <- data[data$time==30,4]
  tiles <- as.data.frame(prev.tiles)
  tiles$eus <- eus
  tiles$exp <- exp.unit
  final.tiles <- rbind(final.tiles,tiles)
}

# prepare data for plotting, including grey frame
for (x in 1:length(final.tiles$prev.tiles)){
  if (final.tiles$prev.tiles[x] < 0.02) {
    final.tiles$shp[x] <- 1 # achieved socially optimal level of control
  } else {
    final.tiles$shp[x] <- 2 # did not achieve socially optimal level of control
  }
}
frames <- final.tiles[final.tiles$shp=="1",]
eust <- c()
for (x in exps) {
  eust[which(exps==x)] <- max(frames[frames$exp==x,]$eus)+(inc.eus/2)
}
df.frame<- data.frame(eus=eust,exp=exps+(inc.exp/2))
for (x in 1:(length(eus)*length(exps))) {
  if (df.frame$eus[x+1]!=df.frame$eus[x] & df.frame$exp[x+1]!=df.frame$exp[x]){
    df.frame <- InsertRow(df.frame, NewRow = c(df.frame$eus[x+1],df.frame$exp[x]), RowNum = x+1)
  }
} # expect a missing value error - it corresponds to the last row
df.frame[length(df.frame$eus)+1,] <- c(min(frames$eus)-(inc.eus/2),max(frames$exp)+(inc.exp/2))
df.frame[length(df.frame$eus)+1,] <- c(min(frames$eus)-(inc.eus/2),min(frames$exp)-(inc.exp/2))
df.frame[length(df.frame$eus)+1,] <- c(df.frame$eus[1],df.frame$exp[1]-inc.exp)
df.frame[length(df.frame$eus)+1,] <- c(df.frame$eus[1],df.frame$exp[1])

# plot tiles with grey frame
fig4 <- ggplot(data = final.tiles, aes(x=eus, y=exp, fill=prev.tiles*100))+ theme_bw()  +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#19DB00", high = "#FC4700", mid = "#FFD000", 
                       midpoint = 5, limit = c(0,10), space = "Lab", 
                       name=expression('Prevalence at 25'^th*' time step (%)')) + 
  labs(y = bquote(paste(italic("Ex post  "), "cost of one infected animal")), x = bquote(paste("Price per control unit"))) + 
  scale_x_continuous(breaks=c(0,0.40,0.85,1.25),labels=c("Free","Low", "Medium", "High")) +
  scale_y_continuous(breaks=c(10,46,85),labels=c("Low", "Medium", "High")) + 
  geom_path(data = df.frame, mapping = aes(x=eus, y=exp, fill=NULL),colour="#807E80",alpha=0.75) +
  theme(text=element_text(size=11,  family="Times New Roman"))
fig4
