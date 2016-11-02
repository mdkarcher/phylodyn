# firstloop <- function(coal_factor, s, event, lengthout)
# {
#     grid <- seq(min(s),max(s),length.out=lengthout+1)
#     u <- diff(grid)
#     field <- grid[-1]-u/2
#     sgrid <- grid
#     event_new <- 0
#     time <- 0
#     where <- 1
#     E.factor <- 0
#     for (j in 1:lengthout)
#     {
#         count <- sum(s>sgrid[j] & s<=sgrid[j+1])
#         if (count>1)
#         {
#             points <- s[s>sgrid[j] & s<=sgrid[j+1]]
#             u <- diff(c(sgrid[j],points))
#             event_new <- c(event_new,event[(where):(where+count-1)])
#             time <- c(time,rep(field[j],count))
#             E.factor <- c(E.factor,coal_factor[where:(where+count-1)]*u)
#             where <- where+count
#             if (max(points)<sgrid[j+1])
#             {
#                 event_new <- c(event_new,0)
#                 time <- c(time,field[j])
#                 E.factor <- c(E.factor,coal_factor[where]*(sgrid[j+1]-max(points)))
#             }
#         }
#         if (count==1)
#         {
#             event_new <- c(event_new,event[where])
#             points <- s[s>sgrid[j] & s<=sgrid[j+1]]
#             if (points==sgrid[j+1])
#             {
#                 E.factor <- c(E.factor,coal_factor[where]*(sgrid[j+1]-sgrid[j]))
#                 time <- c(time,field[j])
#                 where <- where+1
#             }
#             else
#             {
#                 event_new <- c(event_new,0)
#                 E.factor <- c(E.factor,coal_factor[where]*(points-sgrid[j]))
#                 E.factor <- c(E.factor,coal_factor[where+1]*(sgrid[j+1]-points))
#                 time <- c(time,rep(field[j],2))
#                 where <- where+1
#             }
#         }
#         if (count==0)
#         {
#             event_new <- c(event_new,0)
#             E.factor <- c(E.factor,coal_factor[where]*(sgrid[j+1]-sgrid[j]))
#             time <- c(time,field[j])
#         }
#     }
#     return(list(time=time, event_new=event_new, E_factor=E.factor, grid=grid, field=field))
# }
# 
# secondloop <- function(time2, event_new2, E_factor2, lengthout, field)
# {
#     for (j in 1:lengthout)
#     {
#         count <- sum(time2==field[j])
#         if (count>1)
#         {
#             indic <- seq(1:length(event_new2))[time2==field[j]]
#             if (sum(event_new2[indic])==0)
#             {
#                 event_new2 <- event_new2[-indic[-1]]
#                 time2 <- time2[-indic[-1]]
#                 temp <- sum(E_factor2[indic])
#                 E_factor2[indic[1]] <- temp
#                 E_factor2 <- E_factor2[-indic[-1]]
#             }
#             #else {}
#         }
#     }
#     
#     return(list(time2=time2, event_new2=event_new2, E_factor2=E_factor2))
# }
# 
# calculate_moller_hetero <- function(coal.factor, s, event, lengthout,
#                                     prec_alpha = 0.01, prec_beta = 0.01,
#                                     log_zero = -100, alpha = NULL, beta = NULL)
# {
#     if (prec_alpha == 0.01 & prec_beta==0.01 & !is.null(alpha) & !is.null(beta))
#     {
#         prec_alpha <- alpha
#         prec_beta  <- beta
#     }
#     
#     fl <- firstloop(coal_factor = coal.factor, s = s, event = event,
#     lengthout = lengthout)
#     
#     sl <- secondloop(time2 = fl$time, event_new2 = fl$event_new,
#     E_factor2 = fl$E_factor, lengthout = lengthout,
#     field = fl$field)
#     
#     E_log = log(sl$E_factor2)
#     E_log[sl$E_factor2 == 0] = log_zero
#     
#     data <- list(y = sl$event_new2[-1], time = sl$time2[-1], E = E_log[-1])
#     formula <- y ~ -1 + f(time, model="rw1", hyper = list(prec = list(param = c(prec_alpha, prec_beta))), constr = FALSE)
#     mod4 <- INLA::inla(formula, family = "poisson", data = data, offset = E_log, control.predictor = list(compute=TRUE))
#     
#     return(list(result=mod4,grid=fl$grid,data=data,E=E_log, x = fl$field))
# }
# 
# calculate_pref <- function(coal.factor,s,event,lengthout,prec_alpha=0.01,prec_beta=0.01)
# {
#   grid <- seq(0,max(s),length.out=lengthout+1)
#   u <- diff(grid)
#   field <- grid[-1]-u/2
#   #   sgrid <- grid
#   #   event_new <- 0
#   #   time <- 0
#   #   where <- 1
#   #   E.factor <- 0
#   #   for (j in 1:lengthout)
#   #   {
#   #     count <- sum(s>sgrid[j] & s<=sgrid[j+1])
#   #     if (count>1)
#   #     {
#   #       points <- s[s>sgrid[j] & s<=sgrid[j+1]]
#   #       u <- diff(c(sgrid[j],points))
#   #       event_new <- c(event_new,event[(where):(where+count-1)])
#   #       time <- c(time,rep(field[j],count))
#   #       E.factor <- c(E.factor,coal.factor[where:(where+count-1)]*u)
#   #       where <- where+count
#   #       if (max(points)<sgrid[j+1])
#   #       {
#   #         event_new <- c(event_new,0)
#   #         time <- c(time,field[j])
#   #         E.factor <- c(E.factor,coal.factor[where]*(sgrid[j+1]-max(points)))
#   #       }
#   #     }
#   #     if (count==1)
#   #     {
#   #       event_new <- c(event_new,event[where])
#   #       points <- s[s>sgrid[j] & s<=sgrid[j+1]]
#   #       if (points==sgrid[j+1])
#   #       {
#   #         E.factor <- c(E.factor,coal.factor[where]*(sgrid[j+1]-sgrid[j]))
#   #         time <- c(time,field[j])
#   #         where <- where+1
#   #       }
#   #       else
#   #       {
#   #         event_new <- c(event_new,0)
#   #         E.factor <- c(E.factor,coal.factor[where]*(points-sgrid[j]))
#   #         E.factor <- c(E.factor,coal.factor[where+1]*(sgrid[j+1]-points))
#   #         time <- c(time,rep(field[j],2))
#   #         where <- where+1
#   #       }
#   #     }
#   #     if (count==0)
#   #     {
#   #       event_new <- c(event_new,0)
#   #       E.factor <- c(E.factor,coal.factor[where]*(sgrid[j+1]-sgrid[j]))
#   #       time <- c(time,field[j])
#   #     }
#   #
#   #   }
#   #   time2 <- time
#   #   event_new2 <- event_new
#   #   E.factor2 <- E.factor
#   #
#   #   for (j in 1:lengthout)
#   #   {
#   #     count <- sum(time2==field[j])
#   #     if (count>1)
#   #     {
#   #       indic <- seq(1:length(event_new2))[time2==field[j]]
#   #       if (sum(event_new2[indic])==0)
#   #       {
#   #         event_new2 <- event_new2[-indic[-1]]
#   #         time2 <- time2[-indic[-1]]
#   #         temp <- sum(E.factor2[indic])
#   #         E.factor2[indic[1]] <- temp
#   #         E.factor2 <- E.factor2[-indic[-1]]
#   #       }
#   #       #else {}
#   #     }
#   #   }
#   
#   #   E.factor2[1:34] <- rep(1,34)
#   #data <- list(y=event_new2[-1],event=event_new2[-1],time=time2[-1],E=log(E.factor2[-1]))
#   #formula <- y~-1+f(time,model="rw1",hyper=list(prec = list(param = c(prec_alpha, prec_beta))),constr=FALSE)
#   #mod4 <- inla(formula,family="poisson",data=data,offset=E,control.predictor=list(compute=TRUE),...)
#   
#   #   n1 <- length(event_new2[-1])
#   #   n2 <- length(field)
#   #   Y  <- matrix(NA,n1+n2,2)
#   #
#   #dd <- heterochronous.gp.stat(influenza.tree)
#   #s.time <- dd$sample.times
#   #n.sample <- dd$sampled.lineages
#   newcount <- rep(0,length(grid)-1)
#   
#   # MK: added to replace reading from influenza.tree above
#   samps = s[event==0]
#   
#   for (j in 1:length(newcount))
#   {
#     # MK: altered to samps version
#     newcount[j] <- sum(samps > grid[j] & samps <= grid[j+1])
#   }
#   newcount[1] <- newcount[1]+1
#   #newcount[(sum(grid<max(samps))+1):lengthout] <- NA
#   newcount[head(grid, -1) >= max(samps)] <- NA
#   
#   data.sampling2<-data.frame(y=newcount,time=field,E=log(diff(grid)))
#   formula.sampling2=y~1+f(time,model="rw1",hyper = list(prec = list(param = c(prec_alpha, prec_beta))),constr=FALSE)
#   
#   # MK: Added during functionizing
#   E = diff(s) * coal.factor
#   E[1]=1
#   
#   mod.sampling2 <- INLA::inla(formula.sampling2,family="poisson",data=data.sampling2,offset=E,control.predictor=list(compute=TRUE))
#   
#   return(list(result=mod.sampling2,grid=grid, data = data.sampling2))
# }
# 
# calculate_moller_hetero_pref <- function(coal.factor,s,event,lengthout,prec_alpha=0.01,prec_beta=0.01,beta1_prec = 0.001,log_zero=-100,alpha=NULL,beta=NULL)
# {
#     if (prec_alpha == 0.01 & prec_beta==0.01 & !is.null(alpha) & !is.null(beta))
#     {
#         prec_alpha = alpha
#         prec_beta  = beta
#     }
#     
#     grid <- seq(0,max(s),length.out=lengthout+1)
#     u <- diff(grid)
#     field <- grid[-1]-u/2
#     sgrid <- grid
#     event_new <- 0
#     time <- 0
#     where <- 1
#     E.factor <- 0
#     for (j in 1:lengthout)
#     {
#         count <- sum(s>sgrid[j] & s<=sgrid[j+1])
#         if (count>1)
#         {
#             points <- s[s>sgrid[j] & s<=sgrid[j+1]]
#             u <- diff(c(sgrid[j],points))
#             event_new <- c(event_new,event[(where):(where+count-1)])
#             time <- c(time,rep(field[j],count))
#             E.factor <- c(E.factor,coal.factor[where:(where+count-1)]*u)
#             where <- where+count
#             if (max(points)<sgrid[j+1])
#             {
#                 event_new <- c(event_new,0)
#                 time <- c(time,field[j])
#                 E.factor <- c(E.factor,coal.factor[where]*(sgrid[j+1]-max(points)))
#             }
#         }
#         if (count==1)
#         {
#             event_new <- c(event_new,event[where])
#             points <- s[s>sgrid[j] & s<=sgrid[j+1]]
#             if (points==sgrid[j+1])
#             {
#                 E.factor <- c(E.factor,coal.factor[where]*(sgrid[j+1]-sgrid[j]))
#                 time <- c(time,field[j])
#                 where <- where+1
#             }
#             else
#             {
#                 event_new <- c(event_new,0)
#                 E.factor <- c(E.factor,coal.factor[where]*(points-sgrid[j]))
#                 E.factor <- c(E.factor,coal.factor[where+1]*(sgrid[j+1]-points))
#                 time <- c(time,rep(field[j],2))
#                 where <- where+1
#             }
#         }
#         if (count==0)
#         {
#             event_new <- c(event_new,0)
#             E.factor <- c(E.factor,coal.factor[where]*(sgrid[j+1]-sgrid[j]))
#             time <- c(time,field[j])
#         }
#         
#     }
#     time2 <- time
#     event_new2 <- event_new
#     E.factor2 <- E.factor
#     
#     for (j in 1:lengthout)
#     {
#         count <- sum(time2==field[j])
#         if (count>1)
#         {
#             indic <- seq(1:length(event_new2))[time2==field[j]]
#             if (sum(event_new2[indic])==0)
#             {
#                 event_new2 <- event_new2[-indic[-1]]
#                 time2 <- time2[-indic[-1]]
#                 temp <- sum(E.factor2[indic])
#                 E.factor2[indic[1]] <- temp
#                 E.factor2 <- E.factor2[-indic[-1]]
#             }
#             #else {}
#         }
#     }
#     
#     #   E.factor2[1:34] <- rep(1,34)
#     #data <- list(y=event_new2[-1],event=event_new2[-1],time=time2[-1],E=log(E.factor2[-1]))
#     #formula <- y~-1+f(time,model="rw1",hyper=list(prec = list(param = c(prec_alpha, prec_beta))),constr=FALSE)
#     #mod4 <- inla(formula,family="poisson",data=data,offset=E,control.predictor=list(compute=TRUE),...)
#     
#     n1 <- length(event_new2[-1])
#     n2 <- length(field)
#     Y  <- matrix(NA,n1+n2,2)
#     
#     #dd <- heterochronous.gp.stat(influenza.tree)
#     #s.time <- dd$sample.times
#     #n.sample <- dd$sampled.lineages
#     newcount <- rep(0,length(grid)-1)
#     
#     # MK: added to replace reading from influenza.tree above
#     samps = s[-1][event==0]
#     
#     for (j in 1:length(newcount))
#     {
#         # MK: altered to samps version
#         newcount[j] <- sum(samps > grid[j] & samps <= grid[j+1])
#     }
#     newcount[1] <- newcount[1]+1
#     newcount[(sum(grid<max(samps))+1):lengthout] <- NA
#     
#     #data.sampling2<-data.frame(y=newcount,time=field,E=log(diff(grid)))
#     #formula.sampling2=y~1+f(time,model="rw1",hyper = list(prec = list(param = c(prec_alpha, prec_beta))),constr=FALSE)
#     
#     # MK: Added during functionizing
#     E = diff(s) * coal.factor
#     E[1]=1
#     
#     #print("Got here 1")
#     
#     beta0<-c(rep(0,n1),rep(1,n2))
#     newE<-rep(NA,n1+n2)
#     
#     #MK: added to resolve INLA failure
#     E.factor2.log = log(E.factor2)
#     E.factor2.log[E.factor2 == 0] = log_zero
#     newE[1:n1]<-E.factor2.log[-1]
#     #newE[1:n1]<-log(E.factor2[-1])
#     
#     newE[(n1+1):(n1+n2)]<-log(diff(grid))
#     megafield<-c(time2[-1],field)
#     
#     #print("Got here 2")
#     
#     Y[1:n1,1]<-event_new2[-1]
#     #print("Got here 3")
#     Y[(n1+1):(n2+n1),2]<-newcount
#     #print("Got here 4")
#     r<-c(rep(1,n1),rep(2,n2))
#     w<-c(rep(1,n1),rep(-1,n2))
#     # data.pref<-data.frame(Y=Y,beta0=beta0,r=r,megafield=megafield,E=newE,w=w)
#     # formula.pref.rep<-Y~-1+beta0+f(megafield,w,model="rw1",replicate=r,hyper=list(prec = list(param = c(.001, .001))),constr=FALSE)
#     ii<-c(time2[-1],rep(NA,n2))
#     jj<-c(rep(NA,n1),field)
#     
#     #print("Got here 5")
#     
#     formula.pref<-Y~-1+beta0+f(ii,model="rw1",hyper=list(prec = list(param = c(prec_alpha, prec_beta))),constr=FALSE)+
#     f(jj,w,copy="ii",fixed=FALSE,param=c(0,beta1_prec))
#     #print("Got here 6")
#     
#     # MK: Changed data.frame to list
#     data.pref<-list(Y=Y,beta0=beta0,ii=ii,jj=jj,E=newE,w=w)
#     
#     #print("Got here 7")
#     #print(data.pref)
#     mod.pref<-INLA::inla(formula.pref,family=c("poisson","poisson"),offset=E,data=data.pref,control.predictor=list(compute=TRUE))
#     
#     #print("Got here 8")
#     
#     return(list(result=mod.pref,data = data.pref, grid=grid, x = field))
# }
# 
# calculate.moller.hetero.pref = function(...)
# {
#     return(calculate_moller_hetero_pref(...))
# }

# gen_INLA_args_old = function(coal_times, s_times, n_sampled)
# {
#     n         = length(coal_times) + 1
#     data      = matrix(0, nrow=n-1, ncol=2)
#     data[,1]  = coal_times
#     s_times   = c(s_times,max(data[,1])+1)
#     data[1,2] = sum(n_sampled[s_times <= data[1,1]])
#     tt = length(s_times[s_times <= data[1,1]]) + 1
#     
#     for (j in 2:nrow(data))
#     {
#         if (data[j,1] < s_times[tt])
#         {
#             data[j,2] = data[j-1,2]-1
#         }
#         else
#         {
#             data[j,2] = data[j-1,2] - 1 + sum(n_sampled[s_times > data[j-1,1] & s_times <= data[j,1]])
#             tt = length(s_times[s_times <= data[j,1]]) + 1
#         }
#     }
#     
#     s = unique(sort(c(data[,1], s_times[1:length(s_times)-1])))
#     event1 = sort(c(data[,1], s_times[1:length(s_times)-1]), index.return=TRUE)$ix
#     n = nrow(data)+1
#     l = length(s)
#     event = rep(0, l)
#     event[event1<n] = 1
#     
#     y = diff(s)
#     
#     coal.factor = rep(0,l-1)
#     #indicator = rep(0,l-1) # redundant
#     
#     t = rep(0,l-1)
#     indicator = cumsum(n_sampled[s_times<data[1,1]])
#     indicator = c(indicator, indicator[length(indicator)]-1)
#     ini = length(indicator) + 1
#     for (k in ini:(l-1))
#     {
#         j = data[data[,1]<s[k+1] & data[,1]>=s[k],2]
#         if (length(j) == 0)
#         {
#             indicator[k] = indicator[k-1] + sum(n_sampled[s_times < s[k+1] & s_times >= s[k]])
#         }
#         if (length(j) > 0)
#         {
#             indicator[k] = j-1+sum(n_sampled[s_times < s[k+1] & s_times >= s[k]])
#         }
#     }
#     coal_factor = indicator*(indicator-1)/2
#     
#     
#     return(list(coal_factor=coal_factor, s=s, event=event, lineages=indicator))
# }

# calculate_moller_hetero_old <- function(coal.factor,s,event,lengthout,prec_alpha=0.01,prec_beta=0.01,log_zero=-100,alpha=NULL,beta=NULL)
# {
#     if (prec_alpha == 0.01 & prec_beta==0.01 & !is.null(alpha) & !is.null(beta))
#     {
#         prec_alpha = alpha
#         prec_beta  = beta
#     }
#     grid <- seq(0,max(s),length.out=lengthout+1)
#     u <- diff(grid)
#     field <- grid[-1]-u/2
#     sgrid <- grid
#     event_new <- 0
#     time <- 0
#     where <- 1
#     E.factor <- 0
#     for (j in 1:lengthout)
#     {
#         count <- sum(s>sgrid[j] & s<=sgrid[j+1])
#         if (count>1)
#         {
#             points <- s[s>sgrid[j] & s<=sgrid[j+1]]
#             u <- diff(c(sgrid[j],points))
#             event_new <- c(event_new,event[(where):(where+count-1)])
#             time <- c(time,rep(field[j],count))
#             E.factor <- c(E.factor,coal.factor[where:(where+count-1)]*u)
#             where <- where+count
#             if (max(points)<sgrid[j+1])
#             {
#                 event_new <- c(event_new,0)
#                 time <- c(time,field[j])
#                 E.factor <- c(E.factor,coal.factor[where]*(sgrid[j+1]-max(points)))
#             }
#         }
#         if (count==1)
#         {
#             event_new <- c(event_new,event[where])
#             points <- s[s>sgrid[j] & s<=sgrid[j+1]]
#             if (points==sgrid[j+1])
#             {
#                 E.factor <- c(E.factor,coal.factor[where]*(sgrid[j+1]-sgrid[j]))
#                 time <- c(time,field[j])
#                 where <- where+1
#             }
#             else
#             {
#                 event_new <- c(event_new,0)
#                 E.factor <- c(E.factor,coal.factor[where]*(points-sgrid[j]))
#                 E.factor <- c(E.factor,coal.factor[where+1]*(sgrid[j+1]-points))
#                 time <- c(time,rep(field[j],2))
#                 where <- where+1
#             }
#         }
#         if (count==0)
#         {
#             event_new <- c(event_new,0)
#             E.factor <- c(E.factor,coal.factor[where]*(sgrid[j+1]-sgrid[j]))
#             time <- c(time,field[j])
#         }
#         
#     }
#     time2 <- time
#     event_new2 <- event_new
#     E.factor2 <- E.factor
#     
#     for (j in 1:lengthout)
#     {
#         count <- sum(time2==field[j])
#         if (count>1)
#         {
#             indic <- seq(1:length(event_new2))[time2==field[j]]
#             if (sum(event_new2[indic])==0)
#             {
#                 event_new2 <- event_new2[-indic[-1]]
#                 time2 <- time2[-indic[-1]]
#                 temp <- sum(E.factor2[indic])
#                 E.factor2[indic[1]] <- temp
#                 E.factor2 <- E.factor2[-indic[-1]]
#             }
#             #else {}
#         }
#     }
#     
#     #E.factor2[E.factor2 == 0] = exp(-1e6)
#     E.factor2.log = log(E.factor2)
#     E.factor2.log[E.factor2 == 0] = log_zero
#     #print(E.factor2.log)
#     
#     #   E.factor2[1:34] <- rep(1,34)
#     data <- list(y=event_new2[-1],event=event_new2[-1],time=time2[-1],E=E.factor2.log[-1])
#     formula <- y~-1+f(time,model="rw1",hyper=list(prec = list(param = c(prec_alpha, prec_beta))),constr=FALSE)
#     mod4 <- INLA::inla(formula,family="poisson",data=data,offset=E,control.predictor=list(compute=TRUE))
#     
#     return(list(result=mod4,grid=grid,data=data,E=E.factor2.log, x = field))
# }

# calculate_pref_old = function(coal.factor,s,event,lengthout,prec_alpha=0.01,prec_beta=0.01)
# {
#     grid <- seq(0,max(s),length.out=lengthout+1)
#     u <- diff(grid)
#     field <- grid[-1]-u/2
#     sgrid <- grid
#     event_new <- 0
#     time <- 0
#     where <- 1
#     E.factor <- 0
#     for (j in 1:lengthout)
#     {
#         count <- sum(s>sgrid[j] & s<=sgrid[j+1])
#         if (count>1)
#         {
#             points <- s[s>sgrid[j] & s<=sgrid[j+1]]
#             u <- diff(c(sgrid[j],points))
#             event_new <- c(event_new,event[(where):(where+count-1)])
#             time <- c(time,rep(field[j],count))
#             E.factor <- c(E.factor,coal.factor[where:(where+count-1)]*u)
#             where <- where+count
#             if (max(points)<sgrid[j+1])
#             {
#                 event_new <- c(event_new,0)
#                 time <- c(time,field[j])
#                 E.factor <- c(E.factor,coal.factor[where]*(sgrid[j+1]-max(points)))
#             }
#         }
#         if (count==1)
#         {
#             event_new <- c(event_new,event[where])
#             points <- s[s>sgrid[j] & s<=sgrid[j+1]]
#             if (points==sgrid[j+1])
#             {
#                 E.factor <- c(E.factor,coal.factor[where]*(sgrid[j+1]-sgrid[j]))
#                 time <- c(time,field[j])
#                 where <- where+1
#             }
#             else
#             {
#                 event_new <- c(event_new,0)
#                 E.factor <- c(E.factor,coal.factor[where]*(points-sgrid[j]))
#                 E.factor <- c(E.factor,coal.factor[where+1]*(sgrid[j+1]-points))
#                 time <- c(time,rep(field[j],2))
#                 where <- where+1
#             }
#         }
#         if (count==0)
#         {
#             event_new <- c(event_new,0)
#             E.factor <- c(E.factor,coal.factor[where]*(sgrid[j+1]-sgrid[j]))
#             time <- c(time,field[j])
#         }
#         
#     }
#     time2 <- time
#     event_new2 <- event_new
#     E.factor2 <- E.factor
#     
#     for (j in 1:lengthout)
#     {
#         count <- sum(time2==field[j])
#         if (count>1)
#         {
#             indic <- seq(1:length(event_new2))[time2==field[j]]
#             if (sum(event_new2[indic])==0)
#             {
#                 event_new2 <- event_new2[-indic[-1]]
#                 time2 <- time2[-indic[-1]]
#                 temp <- sum(E.factor2[indic])
#                 E.factor2[indic[1]] <- temp
#                 E.factor2 <- E.factor2[-indic[-1]]
#             }
#             #else {}
#         }
#     }
#     
#     #   E.factor2[1:34] <- rep(1,34)
#     #data <- list(y=event_new2[-1],event=event_new2[-1],time=time2[-1],E=log(E.factor2[-1]))
#     #formula <- y~-1+f(time,model="rw1",hyper=list(prec = list(param = c(prec_alpha, prec_beta))),constr=FALSE)
#     #mod4 <- inla(formula,family="poisson",data=data,offset=E,control.predictor=list(compute=TRUE),...)
#     
#     n1 <- length(event_new2[-1])
#     n2 <- length(field)
#     Y  <- matrix(NA,n1+n2,2)
#     
#     #dd <- heterochronous.gp.stat(influenza.tree)
#     #s.time <- dd$sample.times
#     #n.sample <- dd$sampled.lineages
#     newcount <- rep(0,length(grid)-1)
#     
#     # MK: added to replace reading from influenza.tree above
#     samps = s[-1][event==0]
#     
#     for (j in 1:length(newcount))
#     {
#         # MK: altered to samps version
#         newcount[j] <- sum(samps > grid[j] & samps <= grid[j+1])
#     }
#     newcount[1] <- newcount[1]+1
#     newcount[(sum(grid<max(samps))+1):lengthout] <- NA
#     
#     data.sampling2<-data.frame(y=newcount,time=field,E=log(diff(grid)))
#     formula.sampling2=y~1+f(time,model="rw1",hyper = list(prec = list(param = c(prec_alpha, prec_beta))),constr=FALSE)
#     
#     # MK: Added during functionizing
#     E = diff(s) * coal.factor
#     E[1]=1
#     
#     mod.sampling2<-INLA::inla(formula.sampling2,family="poisson",data=data.sampling2,offset=E,control.predictor=list(compute=TRUE))
#     
#     return(list(result=mod.sampling2,grid=grid))
# }

# calculate_moller_hetero_pref_old <- function(coal.factor,s,event,lengthout,prec_alpha=0.01,prec_beta=0.01,beta1_prec = 0.001,log_zero=-100,alpha=NULL,beta=NULL)
# {
#     if (prec_alpha == 0.01 & prec_beta==0.01 & !is.null(alpha) & !is.null(beta))
#     {
#         prec_alpha = alpha
#         prec_beta  = beta
#     }
#     
#     grid <- seq(0,max(s),length.out=lengthout+1)
#     u <- diff(grid)
#     field <- grid[-1]-u/2
#     sgrid <- grid
#     event_new <- 0
#     time <- 0
#     where <- 1
#     E.factor <- 0
#     for (j in 1:lengthout)
#     {
#         count <- sum(s>sgrid[j] & s<=sgrid[j+1])
#         if (count>1)
#         {
#             points <- s[s>sgrid[j] & s<=sgrid[j+1]]
#             u <- diff(c(sgrid[j],points))
#             event_new <- c(event_new,event[(where):(where+count-1)])
#             time <- c(time,rep(field[j],count))
#             E.factor <- c(E.factor,coal.factor[where:(where+count-1)]*u)
#             where <- where+count
#             if (max(points)<sgrid[j+1])
#             {
#                 event_new <- c(event_new,0)
#                 time <- c(time,field[j])
#                 E.factor <- c(E.factor,coal.factor[where]*(sgrid[j+1]-max(points)))
#             }
#         }
#         if (count==1)
#         {
#             event_new <- c(event_new,event[where])
#             points <- s[s>sgrid[j] & s<=sgrid[j+1]]
#             if (points==sgrid[j+1])
#             {
#                 E.factor <- c(E.factor,coal.factor[where]*(sgrid[j+1]-sgrid[j]))
#                 time <- c(time,field[j])
#                 where <- where+1
#             }
#             else
#             {
#                 event_new <- c(event_new,0)
#                 E.factor <- c(E.factor,coal.factor[where]*(points-sgrid[j]))
#                 E.factor <- c(E.factor,coal.factor[where+1]*(sgrid[j+1]-points))
#                 time <- c(time,rep(field[j],2))
#                 where <- where+1
#             }
#         }
#         if (count==0)
#         {
#             event_new <- c(event_new,0)
#             E.factor <- c(E.factor,coal.factor[where]*(sgrid[j+1]-sgrid[j]))
#             time <- c(time,field[j])
#         }
#         
#     }
#     time2 <- time
#     event_new2 <- event_new
#     E.factor2 <- E.factor
#     
#     for (j in 1:lengthout)
#     {
#         count <- sum(time2==field[j])
#         if (count>1)
#         {
#             indic <- seq(1:length(event_new2))[time2==field[j]]
#             if (sum(event_new2[indic])==0)
#             {
#                 event_new2 <- event_new2[-indic[-1]]
#                 time2 <- time2[-indic[-1]]
#                 temp <- sum(E.factor2[indic])
#                 E.factor2[indic[1]] <- temp
#                 E.factor2 <- E.factor2[-indic[-1]]
#             }
#             #else {}
#         }
#     }
#     
#     #   E.factor2[1:34] <- rep(1,34)
#     #data <- list(y=event_new2[-1],event=event_new2[-1],time=time2[-1],E=log(E.factor2[-1]))
#     #formula <- y~-1+f(time,model="rw1",hyper=list(prec = list(param = c(prec_alpha, prec_beta))),constr=FALSE)
#     #mod4 <- inla(formula,family="poisson",data=data,offset=E,control.predictor=list(compute=TRUE),...)
#     
#     n1 <- length(event_new2[-1])
#     n2 <- length(field)
#     Y  <- matrix(NA,n1+n2,2)
#     
#     #dd <- heterochronous.gp.stat(influenza.tree)
#     #s.time <- dd$sample.times
#     #n.sample <- dd$sampled.lineages
#     newcount <- rep(0,length(grid)-1)
#     
#     # MK: added to replace reading from influenza.tree above
#     samps = s[-1][event==0]
#     
#     for (j in 1:length(newcount))
#     {
#         # MK: altered to samps version
#         newcount[j] <- sum(samps > grid[j] & samps <= grid[j+1])
#     }
#     newcount[1] <- newcount[1]+1
#     newcount[(sum(grid<max(samps))+1):lengthout] <- NA
#     
#     #data.sampling2<-data.frame(y=newcount,time=field,E=log(diff(grid)))
#     #formula.sampling2=y~1+f(time,model="rw1",hyper = list(prec = list(param = c(prec_alpha, prec_beta))),constr=FALSE)
#     
#     # MK: Added during functionizing
#     E = diff(s) * coal.factor
#     E[1]=1
#     
#     #print("Got here 1")
#     
#     beta0<-c(rep(0,n1),rep(1,n2))
#     newE<-rep(NA,n1+n2)
#     
#     #MK: added to resolve INLA failure
#     E.factor2.log = log(E.factor2)
#     E.factor2.log[E.factor2 == 0] = log_zero
#     newE[1:n1]<-E.factor2.log[-1]
#     #newE[1:n1]<-log(E.factor2[-1])
#     
#     newE[(n1+1):(n1+n2)]<-log(diff(grid))
#     megafield<-c(time2[-1],field)
#     
#     #print("Got here 2")
#     
#     Y[1:n1,1]<-event_new2[-1]
#     #print("Got here 3")
#     Y[(n1+1):(n2+n1),2]<-newcount
#     #print("Got here 4")
#     r<-c(rep(1,n1),rep(2,n2))
#     w<-c(rep(1,n1),rep(-1,n2))
#     # data.pref<-data.frame(Y=Y,beta0=beta0,r=r,megafield=megafield,E=newE,w=w)
#     # formula.pref.rep<-Y~-1+beta0+f(megafield,w,model="rw1",replicate=r,hyper=list(prec = list(param = c(.001, .001))),constr=FALSE)
#     ii<-c(time2[-1],rep(NA,n2))
#     jj<-c(rep(NA,n1),field)
#     
#     #print("Got here 5")
#     
#     formula.pref<-Y~-1+beta0+f(ii,model="rw1",hyper=list(prec = list(param = c(prec_alpha, prec_beta))),constr=FALSE)+
#     f(jj,w,copy="ii",fixed=FALSE,param=c(0,beta1_prec))
#     #print("Got here 6")
#     
#     # MK: Changed data.frame to list
#     data.pref<-list(Y=Y,beta0=beta0,ii=ii,jj=jj,E=newE,w=w)
#     
#     #print("Got here 7")
#     #print(data.pref)
#     mod.pref<-INLA::inla(formula.pref,family=c("poisson","poisson"),offset=E,data=data.pref,control.predictor=list(compute=TRUE))
#     
#     #print("Got here 8")
#     
#     return(list(result=mod.pref,grid=grid, x = field))
# }

# sampling_old = function(data, para, alg, setting, init, verbose=TRUE)
# {
#   # pass the data and parameters
#   lik_init = data$lik_init # f_offset = data$f_offset
#   Ngrid = lik_init$ng+1
#   alpha = para$alpha
#   beta = para$beta
#   invC = para$invC
#   rtEV = para$rtEV
#   EVC = para$EVC
#   cholC = para$cholC
#   
#   # MCMC sampling setting
#   stepsz = setting$stepsz
#   Nleap  = setting$Nleap
#   if (alg=='aMALA')
#     szkappa = setting$szkappa
#   
#   if (alg=="HMC" | alg == "splitHMC")
#   {
#     rand_leap = setting$rand_leap
#   }
#   
#   
#   # storage of posterior samples
#   NSAMP = setting$NSAMP
#   NBURNIN = setting$NBURNIN
#   samp = matrix(NA,NSAMP-NBURNIN,Ngrid) # all parameters together
#   acpi = 0
#   acpt = 0
#   
#   # initialization
#   theta = init$theta
#   u = init$u
#   du = init$du
#   
#   # start MCMC run
#   start_time = Sys.time()
#   cat('Running ', alg ,' sampling...\n')
#   for(Iter in 1:NSAMP)
#   {
#     if(verbose&&Iter%%100==0)
#     {
#       cat(Iter, ' iterations have been finished!\n' )
#       cat('Online acceptance rate is ',acpi/100,'\n')
#       acpi=0
#     }
#     
#     # sample the whole parameter
#     #tryCatch({res=switch(alg,
#     #                     HMC=eval(parse(text='HMC'))(theta,u,du,function(theta,grad=FALSE)U(theta,lik_init,invC,alpha,beta,grad),stepsz,Nleap,rand_leap),
#     #                     splitHMC=eval(parse(text='splitHMC'))(theta,u,du,function(theta,grad=FALSE)U_split(theta,lik_init,invC,alpha,beta,grad),rtEV,EVC,stepsz,Nleap,rand_leap),
#     #                     MALA=eval(parse(text='MALA'))(theta,u,du,function(theta,grad=FALSE)U(theta,lik_init,invC,alpha,beta,grad),stepsz),
#     #                     aMALA=eval(parse(text='aMALA'))(theta,u,function(theta,grad=FALSE)U_kappa(theta,lik_init,invC,alpha,beta,grad),function(theta)Met(theta,lik_init,invC),szkappa,stepsz),
#     #                     ESS=eval(parse(text='ESS'))(theta[-Ngrid],u,function(f)coal_loglik(lik_init,f),cholC/sqrt(theta[Ngrid])),
#     #                     stop('The algorithm is not in the list!'));
#     #          theta[1:(Ngrid-(alg=='ESS'))]=res$q;u=res[[2]];if(any(grepl(alg,c('HMC','splitHMC','MALA'))))du=res$du;
#     #          acpi=acpi+res$Ind}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
#     res=switch(alg,
#                HMC=eval(parse(text='HMC'))(theta,u,du,function(theta,grad=FALSE)U(theta,lik_init,invC,alpha,beta,grad),stepsz,Nleap,rand_leap),
#                splitHMC=eval(parse(text='splitHMC'))(theta,u,du,function(theta,grad=FALSE)U_split(theta,lik_init,invC,alpha,beta,grad),rtEV,EVC,stepsz,Nleap,rand_leap),
#                MALA=eval(parse(text='MALA'))(theta,u,du,function(theta,grad=FALSE)U(theta,lik_init,invC,alpha,beta,grad),stepsz),
#                aMALA=eval(parse(text='aMALA'))(theta,u,function(theta,grad=FALSE)U_kappa(theta,lik_init,invC,alpha,beta,grad),function(theta)Met(theta,lik_init,invC),szkappa,stepsz),
#                ESS=eval(parse(text='ESS'))(theta[-Ngrid],u,function(f)coal_loglik(lik_init,f),cholC/sqrt(theta[Ngrid])),
#                stop('The algorithm is not in the list!'));
#     theta[1:(Ngrid-(alg=='ESS'))]=res$q
#     u=res[[2]]
#     if(any(grepl(alg,c('HMC','splitHMC','MALA'))))du=res$du
#     acpi=acpi+res$Ind
#     # Gibbs sample kappa for ESS
#     if(alg=='ESS')
#       theta[Ngrid]=rgamma(1,alpha+(Ngrid-1)/2,beta+t(theta[-Ngrid])%*%invC%*%theta[-Ngrid]/2)
#     
#     # save posterior samples after burnin
#     if(Iter>NBURNIN)
#     {
#       samp[Iter-NBURNIN,]<-theta
#       acpt<-acpt+res$Ind
#     }
#     
#   }
#   stop_time = Sys.time()
#   time = stop_time-start_time
#   cat('\nTime consumed : ',time)
#   acpt = acpt/(NSAMP-NBURNIN)
#   cat('\nFinal Acceptance Rate: ',acpt,'\n')
#   
#   return(list(samp=samp,time=time,acpt=acpt))
# }
# 
# coal_lik_init_old = function(samp_times, n_sampled, coal_times, grid)
# {
#   ns = length(samp_times)
#   nc = length(coal_times)
#   ng = length(grid)-1
#   
#   if (length(samp_times) != length(n_sampled))
#     stop("samp_times vector of differing length than n_sampled vector.")
#   
#   if (length(coal_times) != sum(n_sampled) - 1)
#     stop("Incorrect length of coal_times: should be sum(n_sampled) - 1.")
#   
#   if (max(samp_times, coal_times) > max(grid))
#     stop("Grid does not envelop all sampling and/or coalescent times.")
#   
#   t = sort(unique(c(samp_times, coal_times, grid)))
#   l = rep(0, length(t))
#   
#   for (i in 1:ns)
#     l[t >= samp_times[i]] = l[t >= samp_times[i]] + n_sampled[i]
#   
#   for (i in 1:nc)
#     l[t >= coal_times[i]] = l[t >= coal_times[i]] - 1
#   
#   #print(l)
#   
#   if (sum((l < 1) & (t >= min(samp_times))) > 0)
#     stop("Number of active lineages falls below 1 after the first sampling point.")
#   
#   mask = l > 0
#   t = t[mask]
#   l = head(l[mask], -1)
#   
#   gridrep = rep(0, ng)
#   for (i in 1:ng)
#     gridrep[i] = sum(t > grid[i] & t <= grid[i+1])
#   
#   C = 0.5 * l * (l-1)
#   D = diff(t)
#   
#   y = rep(0, length(D))
#   y[t[-1] %in% coal_times] = 1
#   
#   rep_idx = cumsum(gridrep)
#   rep_idx = cbind(rep_idx-gridrep+1,rep_idx)
#   
#   return(list(t=t, l=l, C=C, D=D, y=y, gridrep=gridrep, ns=sum(n_sampled), nc=nc, ng=ng, rep_idx=rep_idx, args=list(samp_times=samp_times, n_sampled=n_sampled, coal_times=coal_times, grid=grid)))
# }
#

# #' Generate inhomogeneous coalescent with heterochronous samples.
# #' 
# #' @param sample a matrix with 2 columns. The first column contains the number 
# #'   of samples collected at the time defined in the second column.
# #' @param traj_inv a function returning the one over the effective population 
# #'   size.
# #' @param upper numeric. The maximum of the traj_inv function on the relevant 
# #'   interval.
# #' @param ... additional parameters to pass to the traj_inv function.
# #'   
# #' @return A list containing coalescent times \code{coal_times}, intercoalescent
# #'   times \code{intercoal_times}, and number of active lineages between
# #'   coalescent times \code{lineages}.
# #' @export
# #' 
# #' @examples
# #' coalgen_thinning_hetero(cbind(c(10,3), c(0,5)), unif_traj_inv, 0.01)
# coalgen_thinning_hetero <- function(sample,traj_inv,upper,...)
# {
#   # sample = is a matrix with 2 columns. The first column contains the number of samples collected at the time defined in the second column
#   # traj_inv = the effective population size function
#   # this works for heterochronous sampling
#   # assumes sample[1,1]>1
#   samp_times = sample[,2]
#   n_sampled = sample[,1]
#   
#   s=sample[1,2]
#   b <- sample[1,1]
#   n <- sum(sample[,1])-1
#   m <- n
#   nsample <- nrow(sample)
#   sample <- rbind(sample,c(0,10*max(sample,2))) # Bug?
#   out <- rep(0,n)
#   branches <- rep(0,n)
#   i <- 1
#   while (i<(nsample))
#   {
#     #if (b==1)
#     #{
#     #  break
#     #}
#     if (b<2)
#     {
#       b <- b+sample[i+1,1]
#       s <- sample[i+1,2]
#       i <- i+1
#     }
#     E <- stats::rexp(1,upper*b*(b-1)*.5)
#     if (stats::runif(1) <= traj_inv(E+s,...)/upper)
#     {
#       if ( (s+E)>sample[i+1,2])
#       {
#         b <- b+sample[i+1,1]
#         s <- sample[i+1,2]
#         i <- i+1
#       }
#       else
#       {
#         s <- s+E
#         out[m-n+1] <- s
#         branches[m-n+1] <- b
#         n <- n-1
#         b <- b-1
#       }
#     }
#     else
#     {
#       s <- s+E
#     }    
#   }
#   
#   while (b>1)
#   { 
#     E <- stats::rexp(1,upper*b*(b-1)*.5)
#     if (stats::runif(1)<=traj_inv(E+s,...)/upper)
#     {
#       s <- s+E
#       out[m-n+1] <- s
#       branches[m-n+1] <- b
#       n <- n-1
#       b <- b-1
#     }
#     else
#     {
#       s <- s+E
#     }
#   }
#   
#   return(list(intercoal_times=c(out[1],diff(out)), lineages=branches, 
#               coal_times=out, samp_times = samp_times, n_sampled = n_sampled))   
# }
# 
# coalgen_thinning_iso <- function(sample,traj_inv,upper=25,...)
# {
#   ###Need to add correction to "systematic" definition of upper bound
#   s=sample[2]
#   n <- sample[1]
#   out <- rep(0,n-1)
#   time <- 0
#   j <- n
#   while (j>1)
#   {
#     time <- time+stats::rexp(1,upper*j*(j-1)*.5)
#     if (stats::runif(1)<=1/(traj_inv(time,...)*upper))
#     {
#       out[n-j+1] <- time
#       j <- j-1
#     }
#   }
#   return(list(intercoal_times=c(out[1],diff(out)),lineages=seq(n,2,-1), coal_times=out))
# }
# 
# coalgen_transformation_hetero <- function(sample, trajectory,val_upper=10)
# {
#   # sample = is a matrix with 2 columns. The first column contains the number of samples collected at the time defined in the second column
#   # trajectory = one over the effective population size function
#   # this works for heterochronous sampling
#   # assumes sample[1,1]>1
#   s=sample[1,2]
#   b <- sample[1,1]
#   n <- sum(sample[,1])-1
#   m <- n
#   nsample <- nrow(sample)
#   sample <- rbind(sample,c(0,10))
#   out <- rep(0,n)
#   branches <- rep(0,n)
#   i <- 1
#   while (i<(nsample+1))
#   {
#     #if (b==1)
#     #{
#     #  break
#     #}
#     if (b<2)
#     {
#       b <- b+sample[i+1,1]
#       s <- sample[i+1,2]
#       i <- i+1
#     }
#     x <- stats::rexp(1)
#     f <- function(bran,u,x,s) .5*bran*(bran-1)*stats::integrate(trajectory, s, s+u)$value - x    
#     y <- stats::uniroot(f,bran=b,x=x,s=s,lower=0,upper=val_upper)$root
#     while ( (s+y)>sample[i+1,2])
#     {
#       #     f <- function(bran,u,x,s) .5*bran*(bran-1)*stats::integrate(trajectory, s, s+u)$value - x    
#       #     y <- stats::uniroot(f,bran=b,x=x,s=s,lower=0,upper=val_upper)$root
#       x <- x-.5*b*(b-1)*stats::integrate(trajectory,s,sample[i+1,2])$value  
#       b <- b+sample[i+1,1]
#       s <- sample[i+1,2]
#       i <- i+1
#       f <- function(bran,u,x,s) .5*bran*(bran-1)*stats::integrate(trajectory, s, s+u)$value - x    
#       y <- stats::uniroot(f,bran=b,x=x,s=s,lower=0,upper=val_upper)$root
#       if (i==nsample)
#       {
#         sample[nsample+1,2] <- 10*(s+y)
#       }
#     } 
#     
#     s <- s+y
#     out[m-n+1] <- s
#     branches[m-n+1] <- b
#     n <- n-1
#     b <- b-1
#     if (i==nsample)
#     {
#       sample[nsample+1,2] <- 10*(s+y)
#     }
#   }
#   
#   return(list(branches=c(out[1],diff(out)),lineages=branches))   
# }
