###### R Bootcamp Code
packages<-c('ncdf4','tidyr','RColorBrewer','ggplot2','class','fields','R.matlab',
            'readxl','gdata','propagate','nlstools','stringi','plyr','XML','biomod2',
            'grid', 'pwr','formattable','dplyr','zoo','lattice','gridExtra','gam','tidyr','cowplot',
            'nls2','reshape','stringr','tidyverse','fishualize')
funlist<-lapply(packages, function(x){
  if(!require(x, character.only = TRUE)){
  install.packages(x, dependencies = TRUE)
  library(x, character.only = TRUE)}else
    require(x,character.only=T)})


###### Basic code example

### Real Data
P_putida_time=c(0,2,3,4,5,6,7,8,9)
P_putida_abundance=c(0.2162, 1.0808, 1.1867, 2.1595, 6.3756, 19.4565, 40.3202, 81.5114, 101.2929)*1e6

pseudo_raw<-ggplot()+
  geom_line(aes(x=P_putida_time,y=P_putida_abundance))+
  geom_point(aes(x=P_putida_time,y=P_putida_abundance))+
  scale_x_continuous(name=bquote(Time~(hours)))+
  scale_y_continuous(name=bquote(italic(Psuedmonas~pudita)~Concentration))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = c(0.2,0.8),
        text = element_text(size=16),
        axis.text = element_text(color='black'))

#### Save the output
ggsave('~/Desktop/raw_psuedo.tiff',pseudo_raw,dpi=300,width=6,height=6)

## Parameters
r= 1.35 #1.3
dt=1 
iter=8
P=vector(length=iter)
P[1]=0.2162*1E6
for (i in 2:iter){
  P[i]=P[i-1]+(dt*P[i-1]*r) 
}
time_vec <- 1:iter
#### Visualize the output
plot(time_vec,P, ty="b", ylab = "P. putida Abundance (cells/mL)", xlab = "Time (hours)")
lines(P_putida_time,P_putida_abundance,ty='b',col='red')
lines(time_vec,P,ty='b',col='blue')
lines(time_vec,P,ty='b',col='green')
lines(time_vec,P,ty='b',col='purple')
lines(time_vec,P,ty='b',col='pink')


### Base-R plotting
tiff(filename = '~/Desktop/combo_psuedo.tiff',dpi=300,width=4,height=4)
plot(c(1:iter),P,ty='b',main=paste('r=',r,sep=''))
lines(P_putida_time,P_putida_abundance,ty='b',col='red')
dev.off()

### Plotting using ggplot2
plot.pseudo<-ggplot()+
  geom_line(aes(x=time_vec,y=P,col='Model'))+
  geom_point(aes(x=time_vec,y=P,col='Model'))+
  geom_line(aes(x=P_putida_time,y=P_putida_abundance,col='Real'))+
  geom_point(aes(x=P_putida_time,y=P_putida_abundance,col='Real'),shape='square')+
  scale_x_continuous(name=bquote(Time~(hours)))+
  scale_y_continuous(name=bquote(italic(Psuedmonas~pudita)~Abundance~(cells/mL)))+
  scale_color_manual(name='',values=c('Real'='black','Model'='red'))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = c(0.2,0.8),
        text = element_text(size=16),
        axis.text = element_text(color='black'))

#### Save the output
ggsave('~/Desktop/combo_psuedo.tiff',plot.pseudo,dpi=300,width=6,height=6)


###### More complex example of community dynamics
n <- 1
dt <- 1/24 # time step in hours converted to days
years <- 5 # how many years to run the model
days <- 365 * years # how many days in total
iter <- days/dt # convert that to time steps
u <-  1.4 #c(1.4,0.6) # doublings per day
k <-  1 #c(1,0.3) # half-saturation constant
m <- 0.3 # fraction of death per day
ga <- 0.6 # predation fraction
ro <- 0.3 # integration of P into Z
yt <- 0.1 # fraction of natural Z death
d <- 0.8 # supply of nutrients from the deep
N.0 <- 0.8 # deep concentration of nutrients
# Set up your initial conditions (i.e. t = 0)
P.init <- 1
Z.init<- 0.5
N.init <- 1
# P <- vector(length=iter)*NA
# P[1]<-P.init
P <- matrix(NA,nrow=iter,ncol=n)
P[1,]<-P.init
Z <- vector(length=iter)*NA
Z[1] <- Z.init
N <- vector(length=iter)*NA
N[1] <- N.init


for (i in 2:iter){
  P[i,] <- P[i-1,] + dt * ((u * P[i-1,] * (N[i-1] / (N[i-1] + k))) - m * P[i-1,] - ga * P[i-1,] * Z[i-1])
  Z[i] <- Z[i-1] + dt * (sum(ro * ga * Z[i-1] * P[i-1,]) - yt * Z[i-1])
  N[i] <- N[i-1] + dt * (d*(N.0-N[i-1])-sum((u * P[i-1,] * (N[i-1] / (N[i-1] + k))) + (1-ro) * ga * P[i-1,] * Z[i-1] + m * P[i-1,],na.rm=T) + yt * Z[i-1])
}


com.mat<-data.frame(time=1:iter,P=P,Z=Z,N=N) [1:(iter/15),]
com.mat.melt<-melt(com.mat,'time')
ggplot(data=com.mat.melt)+
  geom_line(aes(x=time/(24*365),y=value,col=variable))+
  scale_x_continuous(name=bquote(Time~(years)))+
  scale_y_continuous(name=bquote(Concentration~(mmol~P~m^-3)))+
  scale_color_brewer(name=bquote(Variable),palette="Dark2",
                     labels=c('Phyto-1', 'Zoo.', 'Nut.'))+
  theme_bw()+
  theme(panel.grid = element_blank())


### 1) What is a value of r in the 'Simple Model' example that produces the closest fit line to the real data?


### 2) Using the complex model above, test the following hypotheses:
        ## a) The species of phytoplankton with the HIGHEST half-saturation constant will out-compete 
        
### 3) Walk through the steps of creating your own model based of off whatever system you want. 
### Steps.
#### a. Pick the system you want to study.
#### b. Draw out the box model for that system and fill in the arrows and connections that are relevant (Hint: Defining a specific, testable hypothesis is very useful)
#### c. Take a stab at writing out the equations that would work for your example. 
#### d. Set-up the model in R like above using the outline below. 
#### e. Run the model and test your hypothesis


dt <- 1/24 # time step in hours converted to days
years <- 5 # how many years to run the model
days <- 365 * years # how many days in total
iter <- days/dt # convert that to time steps
u <- 1.4 # doublings per day
k <- 0.5 # half-saturation constant
m <- 0.3 # fraction of death per day
ga <- 0.6 # predation fraction
ro <- 0.3 # integration of P into Z
yt <- 0.1 # fraction of natural Z death
d <- 0.8 # supply of nutrients from the deep
N.0 <- 0.8 # deep concentration of nutrients
# Set up your initial conditions (i.e. t = 0)
P.init <- 1
Z.init<- 0.5
N.init <- 1
# Create vectors/matricies to capture the output
P <- vector(length = iter)*NA
P[1]<-P.init
Z <- vector(length = iter)*NA
Z[1] <- Z.init
N <- vector(length = iter)*NA
N[1] <- N.init
# The model
for (i in 2:iter){
  P[i] <- P[i-1] + dt * ((u * P[i-1] * (N[i-1] / (N[i-1] + k))) - m * P[i-1]- ga * P[i-1] * Z[i-1])
  Z[i] <- Z[i-1] + dt * (ro * ga * Z[i-1] * P[i-1] - yt * Z[i-1])
  N[i] <- N[i-1] + dt * (d*(N.0-N[i-1])-(u * P[i-1] * (N[i-1] / (N[i-1] + k))) + (1-ro) * ga * P[i-1] * Z[i-1] + m * P[i-1] + yt * Z[i-1])
}
com.mat<-data.frame(time=1:iter,P=P,Z=Z,N=N) [1:(iter/15),]
com.mat.melt<-melt(com.mat,'time')
ggplot(data=com.mat.melt)+
  geom_line(aes(x=time/(24*365),y=value,col=variable))+
  scale_x_continuous(name=bquote(Time~(years)))+
  scale_y_continuous(name=bquote(Concentration~(mmol~P~m^-3)))+
  scale_color_brewer(name=bquote(Variable),palette="Dark2",
                     labels=c('Phyto-1','Zoo.','Nut.'))+
  theme_bw()+
  theme(panel.grid = element_blank())

###### Model Outline -----
##### Parameters

##### Initial values

##### Equations

##### Visualize
