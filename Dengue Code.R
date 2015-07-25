library(deSolve) #call up deSolve package
par(mfrow = c(1,1)) #set graph window
set.seed(10) #sets a basis for random calculations
t=seq(1,3650,1) #time sequence
Temp=21.5074496176725+6.6022833995137*cos(3.46682877131032-0.0172139654563135*t)+rnorm(sd=0.4636,3650) #temperature model based on averages and standard deviations from Florida
plot(Temp,xlab="Days", ylab="Temperature in Celsius", main="Temperature Model") #plot of temperature variation
Perc= 4.95372370490162 + 2.80054808107155*cos(2.11003353042212 + 0.616284848835621*(t[1:12])) #temperature model based on averages from Florida
plot(Perc,pch=16, main="Precipitation Model", ylab="Inches", xlab="Month")
P=0
P[1:31]=Perc[1]/31
P[32:59]=Perc[2]/28
P[60:90]=Perc[3]/31
P[91:120]=Perc[4]/30
P[121:151]=Perc[5]/31
P[152:181]=Perc[6]/30
P[182:212]=Perc[7]/31
P[213:243]=Perc[8]/31
P[244:273]=Perc[9]/30
P[274:304]=Perc[10]/31
P[305:334]=Perc[11]/30
P[335:365]=Perc[12]/31
for (i in 1:9){
  P[((i-1)*365+366):(366+364+(i-1)*365)]=P[1:365]
}
Max=max(P)
Min=min(P)
Pnorm=(P-Min)/(Max-Min) 
CC = 0 #setting up CC for delayed control 
for (i in 1:3650){ #loop for 3650 days
  if (i < 11){ #first ten days
    CC[i] = 0 #no control
  }else{ #after that
    CC[i] = 0.3 #30% control
  }
}

pars = c( l =0.00007995, #Human Death rate (and Birth rate)       
          a =.005, #bites per day per host
          ti=0.126, #intrinsic incubation rate             
          y =.08, #Human Recovery rate
          d =11.2, #oviposition rate
          ym =.19, #aquatic transition rate
          mua =.02, #aquatic mortality rate
          bm =.5, #Effective Contact Rate
          mu =.055, #vector mortality rate
          te =0.11, #extrinsic incubation rate
          b1= 95, #nilparous oviposition rate
          b2= 75, #parous oviposition rate
          mue=.05, #egg mortality rate
          sigma=.5,#sex ratio
          mul=.08, #Minimum Larval mortality
          kl=250000, #Standard carrying capacity of larvae per hectare
          kp=250000, #Standard carrying capacity of larvae per hectare
          mup=.03, #minimum Pupal mortality
          muem=.1, # emerging adult mortality
          mur=.08, # mortality rate due to searching behavior
          Te=10.4, #minimum temperature for egg development
          yAem=.4,#emergance rate
          yAh=.2, #transition rate to host seeking
          yAo=.2, #transition rate to ovoposition
          Tag=10, #minimum temperature for egg maturation
          TDDe=110, #degree*days necessary for hatching
          TDDag=77, #degree*days necessary for egg maturation
          area = 16834.92272 #area of florida in hectares
)

init.values = c(t= 1, #time
                Sh= 19000000, #susceptible hosts
                Eh = 0, #emerging hosts
                Ih = 0, #infectious hosts
                Rh =0, #removed hosts
                E= 9000000, #eggs of vector
                L = 100000, #larvae of vector
                P= 100000,  #pupae of vector
                aem= 100000, #adult emerging vector
                ah1= 100000, #adult host seeking nulliparous vector
                ag1= 100000, #adult gravid nulliparous vector
                ao1= 100000, #adult ovipositing nulliparous vector
                ah2= 100000, #adult host seeking parous vector
                ag2= 100000, #adult gravid nulliparous vector
                ao2= 100000, #adult ovipositing nulliparous vector
                Iv = 100000  #infectious vector
)
times = seq(1,3650, by=1) #time sequence for SEIR model

SEIR = function(time, y.values, parameters,ca,cm) { #creates an function base on time, y values, parameters, and the two control mechanisms
  with(as.list(c(y.values, parameters)), {
    dt.dt= 1 # differential equation for time
    dSh.dt = l*(Sh+Eh+Ih+Rh)-a*bm*Sh*Iv/(Sh+Eh+Ih+Rh)-l*Sh # differential equation for susceptible host
    dEh.dt = a*bm*Sh*Iv/(Sh+Eh+Ih+Rh) -(ti+l)*Eh # differential equation for emergin host
    dIh.dt = ti*Eh-(y+l)*Ih # differential equation for infectious host
    dRh.dt = y*Ih-l*Rh # differential equation for recovered host
    dE.dt = yAo*(b1*ao1+b2*ao2)-(mue+(Temp[t]-Te)/(TDDe))*E # differential equation for vector eggs
    dL.dt = E*((Temp[t]-Te)/(TDDe))-((exp(-Temp[t]/2)+ca+mul)*(1+L/(area*kl*(1+Pnorm[t])))+(-.0007*Temp[t]^2+.0392*Temp[t]-.3901))*L # differential equation for vector larvae
    dP.dt = (-.0007*Temp[t]^2+.0392*Temp[t]-.3901)*L-((exp(-Temp[t]/2)+mup)+(0.0008*Temp[t]^2-0.0051*Temp[t]+0.0319))*P # differential equation for vector pupae
    daem.dt = (0.0008*Temp[t]^2-0.0051*Temp[t]+0.0319)*P*sigma*exp(-muem*(1+P/(area*kp*(1+Pnorm[t]))))-((0.04417 + 0.00217*Temp[t])+yAem)*aem # differential equation for adult emerging vector
    dah1.dt =yAem*aem-((0.04417 + 0.00217*Temp[t])+mur+yAh+(a*bm*Ih))*ah1 # differential equation for adult host seeking nulliparous vector
    dag1.dt =yAh*ah1-((0.04417 + 0.00217*Temp[t])+((Temp[t]-Tag)/TDDag))*ag1 # differential equation for adult gravid nulliparous vector
    dao1.dt =((Temp[t]-Tag)/TDDag)*ag1-((0.04417 + 0.00217*Temp[t])+mur+yAo)*ao1 # differential equation for adult ovipositing nulliparous vector
    dah2.dt =yAo*(ao1+ao2)-((0.04417 + 0.00217*Temp[t])+mur+yAh+(a*bm*Ih))*ah2 # differential equation for adult host seeking parous vector
    dag2.dt =yAh*ah2-((0.04417 + 0.00217*Temp[t])+((Temp[t]-Tag)/TDDag))*ag2 # differential equation for adult gravid parous vector
    dao2.dt =(((Temp[t]-Tag)/TDDag))*ag2-((0.04417 + 0.00217*Temp[t])+mur+yAo)*ao2 # differential equation for adult ovipositing parous vector
    dIv.dt = a*bm*Ih*(ah1+ah2)-(mu+cm)*Iv  # differential equation for adult infectious vector
    return(list(c(dt.dt, dSh.dt, dEh.dt, dIh.dt, dRh.dt, dE.dt, dL.dt, dP.dt, daem.dt, dah1.dt, dag1.dt, dao1.dt, dah2.dt, dag2.dt, dao2.dt , dIv.dt))) #returns a list of the state variables
  })
}

  
out1 <-as.data.frame(ode(func = SEIR, y = init.values, 
                        parms = pars, times = times,ca=0,cm=0)) #model run without mosquito control
out2<-as.data.frame(ode(func = SEIR, y = init.values, 
                        parms = pars, times = times,ca=.3,cm=.3)) #model run with 30% mosquito control for larvicide and adulticide
out3<-as.data.frame(ode(func = SEIR, y = init.values, 
                        parms = pars, times = times,ca=0,cm=.3)) #model run 30% mosquito for adulticide only
out4<-as.data.frame(ode(func = SEIR, y = init.values, 
                        parms = pars, times = times,ca=.3,cm=0)) #model run 30% mosquito for larvicide only
out5<-as.data.frame(ode(func = SEIR, y = init.values, 
                        parms = pars, times = times,ca=CC[i],cm=CC[i])) #model run with delayed control bassed on CC



matplot(out1[1], out1[,c(7:11)], type = "l", xlab = "Time (Days)", 
        ylab = "Aedes", main = "Aedes stages",lwd = 3, ylim=c(0,3000000),
        lty = 1, xlim=c(0,3650), col=c(2:12)) #plot of vector stages, no control
legend("topright", c("Egg", "Larva", "Pupa","Adult Emerging","Adult Host1"),
       col = 2:12, lty = 1, lwd = 3) #legend for vector plot
matplot(out1[1], out1[,c(11:14)], type = "l", xlab = "Time (Days)", 
        ylab = "Aedes", main = "Aedes stages", lwd = 3,
        lty = 1, xlim=c(0,3650), col=c(1:5,8), ylim=c(0,15000)) #second plot of vector stages, no control
legend("topright", c("Adult Gravid1","Adult Ovoposting1","Adult Host2",
        "Adult Gravid2","Adult Ovoposting2","Infectious"), col = c(1:5,8), lty = 1, lwd = 3) #legend for second vector plot

matplot(out1[1], out1[5], type = "l", xlab = "Time (Days)", 
        ylab = "People", main = "No Control", lwd = 3,ylim=c(0,15000),
        lty = 1, xlim=c(0,3650), col=2) # plot for infectious humans, no control

matplot(out2[1], out2[5], type = "l", xlab = "Time (Days)", 
        ylab = "Infectious Humans", main = "30% Control", lwd = 3,ylim=c(0,15000),
        lty = 1, xlim=c(0,3650), col=2) #plot for infectious humans, 30% of each control

matplot(out3[1], out3[5], type = "l", xlab = "Time (Days)", 
        ylab = "Infectious Humans", main = "No Larvicide, 30% Pesticide", lwd = 3,ylim=c(0,15000),
        lty = 1, xlim=c(0,3650), col=2) #plot for infectious humans, 30% of pesticide control only

matplot(out4[1], out4[5], type = "l", xlab = "Time (Days)", 
        ylab = "Infectious Humans", main = "No Pesticide, 30% Larvicide", lwd = 3,ylim=c(0,15000),
        lty = 1, xlim=c(0,3650), col=2) #plot for infectious humans, 30% of larvicide control only

matplot(out5[1], out5[5], type = "l", xlab = "Time (Days)", 
        ylab = "Infectious Humans", main = "10 Day Delayed 30% Control", ylim=c(0,15000), lwd = 3,
        lty = 1, xlim=c(0,3650), col=2) #plot for infectious humans, 30% of each control implemented after 10 days


#Ro
a =0.005 #rate of mosquito bites per day
bm =0.5 #proportion of disease transmitting bites
mu =0.055 #infectious mosquito death rate
yy =0.08 #human recovery rate
ah1=100000 #inital number of host seeking nulliparous adult mosquitos
ah2=100000 #inital number of host seeking parous adult mosquitos
R0 = (a^2*bm^2*(ah1+ah2))/(mu*yy) # derived formula for R0
R0


library(rootSolve) #calling up the rootSolve package

rhs <- function(p)  {

  l =0.00007995 #Human Death rate (and Birth rate)       
  a =.005 #bites per day per vector per host
  ti=0.126 #intrinsic incubation rate             
  yy =.08 #Human Recovery rate
  d =11.2 #oviposition rate
  ym =.19 #aquatic transition rate
  mua =.02 #aquatic mortality rate
  ca = 0 #control effort rates
  cm= 0 #control effort rates
  C =1 #Vector Carrying Capacity
  bm =.5 #Effective Contact Rate
  mu =.055 #vector mortality rate
  te =0.11 #extrinsic incubation rate
  b1= 95 #nilparous oviposition rate
  b2= 75 #parous oviposition rate
  mue=.05 #egg mortality rate
  sigma=.5#sex ratio
  mul=.08 #Minimum Larval mortality
  kl=250000 #Standard carrying capacity of larvae per hectare
  kp=250000 #Standard carrying capacity of larvae per hectare
  mup=.03 #minimum Pupal mortality
  muem=.1 # emerging adult mortality
  mur=.08 # mortality rate due to searching behavior
  Te=10.4 #minimum temperature for egg development
  yAem=.4#emergance rate
  yAh=.2 #transition rate to host seeking
  yAo=.2 #transition rate to ovoposition
  Tag=10 #minimum temperature for egg maturation
  TDDe=110 #degree*days necessary for hatching
  TDDag=77 #degree*days necessary for egg maturation
  area = 16834.92272 #area of Florida
  
  # We declare that the state variables are passed to the function rhs in the vector p
  Sh<- p[1]  
  Eh <- p[2]
  Ih <- p[3]
  Rh <- p[4]
  E <- p[5]
  L <- p[6]
  P <- p[7]
  aem <- p[8]
  ah1<- p[9]
  ag1 <- p[10]
  ao1 <- p[11]
  ah2 <- p[12]
  ag2 <- p[13]
  ao2 <- p[14]
  Iv <- p[15]
  
  r <- numeric()  # r will hold the values returned by the rhs on input p
  
# the same differential equations are used as in the SEIR model and passed to the r vectors
  r[1] = l*(Sh+Eh+Ih+Rh)-a*bm*Sh*Iv/(Sh+Eh+Ih+Rh)-l*Sh
  r[2] = a*bm*Sh*Iv/(Sh+Eh+Ih+Rh) -(ti+l)*Eh
  r[3] = ti*Eh-(yy+l)*Ih
  r[4] = yy*Ih-l*Rh
  r[5] = yAo*(b1*ao1+b2*ao2)-(ca+mue+(Temp[t]-Te)/(TDDe))*E
  r[6] = E*((Temp[t]-Te)/(TDDe))-((exp(-Temp[t]/2)+mul)*(1+L/(area*kl*(1+Pnorm[t])))+(-.0007*Temp[t]^2+.0392*Temp[t]-.3901))*L
  r[7] = (-.0007*Temp[t]^2+.0392*Temp[t]-.3901)*L-((exp(-Temp[t]/2)+mup)+(0.0008*Temp[t]^2-0.0051*Temp[t]+0.0319))*P
  r[8] = (0.0008*Temp[t]^2-0.0051*Temp[t]+0.0319)*P*sigma*exp(-muem*(1+P/(area*kp*(1+Pnorm[t]))))-((0.04417 + 0.00217*Temp[t])+yAem)*aem 
  r[9] =yAem*aem-((0.04417 + 0.00217*Temp[t])+mur+yAh+(a*bm*Ih))*ah1
  r[10] =yAh*ah1-((0.04417 + 0.00217*Temp[t])+((Temp[t]-Tag)/TDDag))*ag1
  r[11] =((Temp[t]-Tag)/TDDag)*ag1-((0.04417 + 0.00217*Temp[t])+mur+yAo)*ao1
  r[12] =yAo*(ao1+ao2)-((0.04417 + 0.00217*Temp[t])+mur+yAh+(a*bm*Ih))*ah2
  r[13] =yAh*ah2-((0.04417 + 0.00217*Temp[t])+((Temp[t]-Tag)/TDDag))*ag2
  r[14] =(((Temp[t]-Tag)/TDDag))*ag2-((0.04417 + 0.00217*Temp[t])+mur+yAo)*ao2
  r[15] = a*bm*Ih*(ah1+ah2)-(mu+cm)*Iv
  
  return(r)  # Return the value of r
}

p0 <- c( Sh = 1, Eh=0, Ih = 0, Rh = 0, E=0,L=0,P=0,aem=0,ah1=0,ag1=0,ao1=0,ah2=0,ag2=0,ao2=0,Iv=0)  # An initial guess as to the solution

ans <- multiroot(f = rhs, start = p0) # This calls the solver

ans #gives the values for the equilibrium for each state variable 

Shstar = ans$root[1] # equilibrium value for susceptible humans
Ehstar = ans$root[2] # equilibrium value for emerging humans
Ihstar = ans$root[3] # equilibrium value for infectious humans
Rhstar = ans$root[4] # equilibrium value for recovered humans
Estar = ans$root[5]  # equilibrium value for vector eggs
Lhstar = ans$root[6] # equilibrium value for vector larvae
Phstar = ans$root[7] # equilibrium value for vector pupae
aemstar = ans$root[8] # equilibrium value for vector emerging adults
ah1star = ans$root[9]  # equilibrium value for vector adult host seeking nulliparous
ag1star = ans$root[10]  # equilibrium value for vector adult gravid nulliparous
ao1star = ans$root[11]  # equilibrium value for vector adult ovipositing nulliparous
ah2star = ans$root[12]  # equilibrium value for vector adult host seeking parous
ag2star = ans$root[13]  # equilibrium value for vector adult gravid parous
ao2star = ans$root[14]  # equilibrium value for vector adult ovipositing parous
Ivstar = ans$root[15]  # equilibrium value for vector adult infectious

mod <- function (t = 0, y = c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0), parms = NULL) #setting a function to 
{
  #Define all of the parameters within the fucntion
  l =0.00007995 #Human Death rate (and Birth rate)       
  a =.005 #bites per day per vector per host
  ti=0.126 #intrinsic incubation rate             
  yy =.08 #Human Recovery rate
  d =11.2 #oviposition rate
  ym =.19 #aquatic transition rate
  mua =.02 #aquatic mortality rate
  ca = 0 #control effort rates
  cm= 0 #control effort rates
  C =1 #Vector Carrying Capacity
  bm =.5 #Effective Contact Rate
  mu =.055 #vector mortality rate
  te =0.11 #extrinsic incubation rate
  b1= 95 #nilparous oviposition rate
  b2= 75 #parous oviposition rate
  mue=.05 #egg mortality rate
  sigma=.5#sex ratio
  mul=.08 #Minimum Larval mortality
  kl=250000 #Standard carrying capacity of larvae per hectare
  kp=250000 #Standard carrying capacity of larvae per hectare
  mup=.03 #minimum Pupal mortality
  muem=.1 # emerging adult mortality
  mur=.08 # mortality rate due to searching behavior
  Te=10.4 #minimum temperature for egg development
  yAem=.4#emergance rate
  yAh=.2 #transition rate to host seeking
  yAo=.2 #transition rate to ovoposition
  Tag=10 #minimum temperature for egg maturation
  TDDe=110 #degree*days necessary for hatching
  TDDag=77 #degree*days necessary for egg maturation
  area = 16834.92272 #area of Florida
  TTT=Temp[1] #temperature set at a single value to evaluate the equlibrium
  PPP=Pnorm[1] #precipitation set at a single value to evaluate the equlibrium
  
  #setting all the state variables at y components
  Sh <- y[1]
  Eh <- y[2]
  Ih <- y[3]
  Rh <- y[4]
  E <- y[5]
  L <- y[6]
  P <- y[7]
  aem <- y[8]
  ah1 <- y[9]
  ag1 <- y[10]
  ao1 <- y[11]
  ah2 <- y[12]
  ag2 <- y[13]
  ao2 <- y[14]
  Iv <- y[15]
  
  r <- numeric()
  
  #setting all of the differential equations as r components
  r[1] = l*(Sh+Eh+Ih+Rh)-a*bm*Sh*Iv/(Sh+Eh+Ih+Rh)-l*Sh
  r[2] = a*bm*Sh*Iv/(Sh+Eh+Ih+Rh) -(ti+l)*Eh
  r[3] = ti*Eh-(yy+l)*Ih
  r[4] = yy*Ih-l*Rh
  r[5] = yAo*(b1*ao1+b2*ao2)-(mue+(TTT-Te)/(TDDe))*E
  r[6] = E*((TTT-Te)/(TDDe))-((exp(-TTT/2)+ca+mul)*(1+L/(area*kl*(1+PPP)))+(-.0007*TTT^2+.0392*TTT-.3901))*L
  r[7] = (-.0007*TTT^2+.0392*TTT-.3901)*L-((exp(-TTT/2)+mup)+(0.0008*TTT^2-0.0051*TTT+0.0319))*P
  r[8] = (0.0008*TTT^2-0.0051*TTT+0.0319)*P*sigma*exp(-muem*(1+P/(area*kp*(1+PPP))))-((0.04417 + 0.00217*TTT)+yAem)*aem 
  r[9] =yAem*aem-((0.04417 + 0.00217*TTT)+mur+yAh+(a*bm*Ih))*ah1
  r[10] =yAh*ah1-((0.04417 + 0.00217*TTT)+((TTT-Tag)/TDDag))*ag1
  r[11] =((TTT-Tag)/TDDag)*ag1-((0.04417 + 0.00217*TTT)+mur+yAo)*ao1
  r[12] =yAo*(ao1+ao2)-((0.04417 + 0.00217*TTT)+mur+yAh+(a*bm*Ih))*ah2
  r[13] =yAh*ah2-((0.04417 + 0.00217*TTT)+((TTT-Tag)/TDDag))*ag2
  r[14] =(((Te-Tag)/TDDag))*ag2-((0.04417 + 0.00217*TTT)+mur+yAo)*ao2
  r[15] = a*bm*Ih*(ah1+ah2)-(mu+cm)*Iv
  
  return(r)
}

JJ <- jacobian.full(y = c(Shstar,
                          Ehstar,
                          Ihstar,
                          Rhstar,
                          Estar,
                          Lhstar,
                          Phstar,
                          aemstar,
                          ah1star,
                          ag1star,
                          ao1star,
                          ah2star,
                          ag2star,
                          ao2star,
                          Ivstar
                          ), func = mod)# Jacobian matrix of the equilibrium point
eigen(JJ) #finds the eigen values and vector of the jacobian in order to determine the stability of the equilibrium point