#' Power function boosting and waning model
#' 
#' Solve the model described in Teunis et al. 2016 (Epidemics) and used in Aiemjoy et al. 2022 (Lancet Microbe). There is an addition of a titer-dependent boosting term as used in Ranjeva et al. 2019 (Nature Communications). 
#' @param times vector of times to solve model over
#' @param pars vector of named model parameters
#' @details Model parameters:
#' \itemize{
#'  \item{"t_peak"}{time of peak antibody level post exposure}
#'  \item{"peak"}{peak antibody level on natural scale}
#'  \item{"k"}{antibody level dependent boosting term}
#'  \item{"y0"}{starting antibody level on natural scale -- Note this needs to start at 1}
#'  \item{"v"}{antibody decay rate}
#'  \item{"r"}{shape of antibody decay function}
#'  }
#' @family kinetics_models
#' @return vector of solved antibody levels
#' @export
#' @examples
#' N <- 25
#' r <- rlnorm(N, log(1.3),0.01) ## Shape parameter
#' v <- rexp(N,500) ## Decay rate
#' y0 <- 1 ## Start titre
#' k <- 0 ## Titre dependency
#' peak <- 10^(rnorm(N,5, 1.5)) ## Peak titre
#' 
#' times <- seq(0,365,by=0.1)
#' ys <- matrix(nrow=N,ncol=length(times))
#' for(i in 1:N){
#'   ys[i,] <- kinetics_power_function(times, c(t_peak=21,peak=peak[i],r=r[i],v=v[i],y0=1,k=0))
#' }
#' ys <- log10(ys)
#' 
#' plot(times,(ys[1,]),type='l',ylim=c(0,10))
#' for(i in 2:N){
#'   lines(times,(ys[i,]))
#' }
kinetics_power_function <- function(times, pars){
  t1 <- times <= pars["t_peak"]
  t2 <- times > pars["t_peak"]
  y <- numeric(length(times))
  peak <- pars["peak"]
  k <- pars["k"]
  y0 <- pars["y0"]
  r <- pars["r"]
  v <- pars["v"]
  t_peak <- pars["t_peak"]
  ## Adjust to titre ceiling effect
  peak <- (peak - 1) * exp(-y0*k) + 1
  mu <- (1/t_peak) * log(peak)
  
  peak_reached <- peak*y0
  y[t1] <- y0*exp(mu*times[t1])
  y[t2] <- peak_reached * (1 + v*(r - 1)*(peak_reached^(r-1))*(times[t2]-t_peak))^(-1/(r-1))
  y
}



#' Internal boost function for Neil's simple boosting/waning model
boost_function <- function(times, pars){
  a <- pars["a"]
  b <- pars["b"]
  t_peak <- pars["t_peak"]
  y <- ((a+b)/b) * exp(a*(times-t_peak))
  y
}
#' Internal wane function for Neil's simple boosting/waning model, with added asymptote
#' @export
wane_function_asymptote <- function(times, pars){
  a <- pars["a"]
  b <- pars["b"]
  c <- pars["c"]
  t_peak <- pars["t_peak"]
  y <- ((a+b)/a - c) * exp(-b*(times-t_peak)) + c
  y
}
#' Internal wane function for Neil's simple boosting/waning model
#' @export
wane_function <- function(times, pars){
  a <- pars["a"]
  b <- pars["b"]
  t_peak <- pars["t_peak"]
  y <- ((a+b)/a) * exp(-b*(times-t_peak))
  y
}

#' Simple boosting and waning kinetics
#' 
#' Implements the biphasic boosting and waning model from Singanayagam et al. 2021 (Lancet ID). A phenomonelogical assuming exponential growth and decline, where the weighting of the two components shifts from growth to decline, swapping over at the peak time. Can have an asymptote or not in the waning rate. NOTE: when using the asymptote waning function, the function *does not* peak at t_peak, but slightly after.
#' @inheritParams kinetics_power_function
#' @param boost pointer to the function used to calculate boosting (defaults to boost_function)
#' @param wane pointer to the function used to calculate waning (defaults to wane_function_asymptote)
#' @details Model parameters:
#' \itemize{
#'  \item{"t_peak"}{time of peak antibody level post exposure}
#'  \item{"ymax"}{peak antibody level on natural scale}
#'  \item{"a"}{boosting rate}
#'  \item{"b"}{waning rate}
#'  \item{"c"}{asymptote of waning (can set to 0, or use wane_function instead of wane_function_asymptote to decline to 0)}
#'  }
#' @family kinetics_models
#' @return vector of antibody levels
#' @export
#' @examples
#' pars<- c(a=0.3,b=0.05,c=0.75,ymax=100,t_peak=18)
#' times <- seq(0,100,by=1)
#' y <- kinetics_simple(times,pars)
#' plot(times, y, type='l')
kinetics_simple <- function(times, pars, boost=boost_function, wane=wane_function_asymptote){
  y_boost <- boost(times, pars)
  y_wane <- wane(times, pars)
  pars["ymax"]*y_boost*y_wane/(y_boost + y_wane)
}


#' Internal boosting function for Ranjeva et al.
#' @export
f_rise <- function(times, pars){
  y0 <- pars["y0"]
  peak <- pars["peak"]
  k <- pars["k"]
  peak <- exp(-y0*pars["k"])*peak
  a <- (1/pars["t_peak"]) * log(peak + 1)
  exp(a*(times)) + y0 - 1
}
#' Internal waning function for Ranjeva et al.
#' @export
f_wane <- function(times, pars){
  y_peak <- f_rise(pars["t_peak"], pars)
  y_baseline <- pars["y0"] + pars["frac"]*pars["peak"]*exp(-pars["y0"]*pars["k"])
  w <- pars["w"]
  (y_peak - y_baseline)*exp(-w*(times-pars["t_peak"]))
}

#' Boosting and waning model with ceiling effect and asymptoting boost
#' 
#' Implements the antibody kinetics model from Ranjeva et al. 2019 (Nature Communications). Key features of this model are a short-term, titer-dependent boost that asymptotes towards a peak value rather than continuing to grow exponentially. The waning phase is biphasic -- an initial proportion of the boost is lost and the asymptotes to a steady state.
#' @inheritParams kinetics_simple
#' @param boost pointer to the antibody boosting function (defaults to f_rise)
#' @param wane pointer to the antibody waning function (defaults to f_wane)
#' #' @details Model parameters:
#' \itemize{
#'  \item{"t_peak"}{time of peak antibody level post exposure}
#'  \item{"peak"}{peak antibody level on natural scale}
#'  \item{"k"}{antibody level dependent boosting term}
#'  \item{"y0"}{starting antibody level on natural scale}
#'  \item{"w"}{antibody decay rate}
#'  \item{"frac"}{proportion of the intial boost lost during short-term waning}
#'  }
#' @family kinetics_models
#' @export
#' @return vector of antibody levels
#' @examples
#' pars <- c(y0=1,peak=80,k=0,t_peak=28,w=0.008,frac=0.5)
#' times <- seq(0,1000,by=0.1)
#' y <- kinetics_ranjeva(times,pars)
#' plot(times, y, type='l')
kinetics_ceiling <- function(times, pars, boost=f_rise, wane=f_wane){
  y0 <- pars["y0"]
  y_baseline <- y0 + pars["frac"]*pars["peak"]*exp(-y0*pars["k"])
  y <- numeric(length(times))
  t1 <- times < pars["t_peak"]
  t2 <- times >= pars["t_peak"]
  y[t1] <- boost(times[t1], pars)
  y[t2] <- y_baseline + wane(times[t2], pars)
  y
}


#' Boosting and waning model based on within-host ODEs
#' 
#' Implements the analytical solution to the ODE model described in Pelleau et al. 2021 (JID).
#' @inheritParams kinetics_power_function
#' @details Model parameters:
#' \itemize{
#'  \item{"A"}{initial antibody level}
#'  \item{"beta"}{growth rate of antibody secreting cells}
#'  \item{"rho"}{proportion of B cells which are short-lived}
#'  \item{"delta"}{time at which antibody levels start increasing}
#'  \item{"c_l"}{decay rate of long-lived B cells}
#'  \item{"c_s"}{decay rate of short-lived B cells}
#'  \item{"r"}{antibody decay rate}
#'  }
#' @family kinetics_models
#' @export
#' @return vector of antibody levels
#' @examples
#' pars <- c("A"=0,"beta"=100,"rho"=0.25,"c_s"=0.15,"c_l"=0.0075,"delta"=0,"r"=0.001)
#' times <- seq(0,1000,by=0.1)
#' y <- kinetics_ode(times,pars)
#' plot(times, y, type='l')
kinetics_simple_ode <- function(times, pars){
  A <- pars["A"]
  beta <- pars["beta"]
  rho <- pars["rho"]
  c_s <- pars["c_s"]
  c_l <- pars["c_l"]
  delta <- pars["delta"]
  r <- pars["r"]
  
  short_term <- exp(-c_s*(times-delta))- exp(-r*(times-delta))
  long_term <- exp(-c_l*(times-delta))- exp(-r*(times-delta))
  
  y <- A + beta * ( rho * short_term / (r-c_s) + (1-rho)*long_term/(r-c_l))
  y
}


#' Gamma function kinetics
#' 
#' Very simple scaled gamma function to represent kinetics. Used in Zhao et al. 2018, but not particular interpretable otherwise.
#' @inheritParams kinetics_power_function
#' @details Model parameters:
#' \itemize{
#'  \item{"y0"}{initial antibody level}
#'  \item{"peak"}{parameter to scale function in y-axis}
#'  \item{"shape"}{gamma distribution shape}
#'  \item{"scale"}{gamma distribution scale}
#'  }
#' @family kinetics_models
#' @export
#' @return vector of antibody levels
#' @examples
#' pars <- c("y0"=0,"peak"=300,"shape"=1.2,"scale"=50)
#' times <- seq(0,1000,by=0.1)
#' y <- kinetics_gamma(times,pars)
#' plot(times, y, type='l')
kinetics_gamma <- function(times, pars){
  pars["y0"] + pars["peak"]*dgamma(times, shape=pars["shape"],scale=pars["scale"])
}


#' Piecewise linear model of boosting and biphasic waning
#' 
#' Simple piecewise linear model of linear antibody boosting on log scale followed by linear biphasic waning. Model as in Hay et al. 2019 (PLOS Comp Biol).
#' @inheritParams kinetics_power_function
#' @details Model parameters:
#' \itemize{
#'  \item{"t_i"}{Timing of infection}
#'  \item{"tp"}{Time of peak antibody level post infection}
#'  \item{"mu"}{Magnitude of boost}
#'  \item{"dp"}{Proportion of boost lost in short-term waning}
#'  \item{"ts"}{Duration of short term waning phase}
#'  \item{"m"}{Long-term antibody waning rate}
#'  \item{"lower_bound"}{Lower bound of allowable antibody levels}
#'  \item{"y0"}{True starting titer level}
#'  }
#' @family kinetics_models
#' @export
#' @return vector of antibody levels
#' @examples
#' pars <- c("tp"=21,"mu"=8,"dp"=0.2,"ts"=50,"m"=0.0075,"lower_bound"=0, "y0"=0)
#' times <- seq(0,100,by=0.1)
#' y <- kinetics_piecewise(times,pars)
#' plot(times, y, type='l')
kinetics_piecewise <- function(times, pars){
  t_i <- 0 # pars["t_i"]
  tp <- pars["tp"]
  mu <- pars["mu"]
  dp <- pars["dp"]
  ts <- pars["ts"]
  m <- pars["m"]
  lower_bound <- pars["lower_bound"]
  eff_y0 <- pars["y0"]
  
  y <- numeric(length(times))
  i <- 1
  ## Loops through all times and calculate titre based on time relative to time of infection
  for(t in times){
    tmp <- 0
    if(t <= t_i) tmp = 0
    else if(t > t_i & t <= (t_i + tp)) tmp = (mu/tp)*(t-t_i)
    else if(t > (tp+t_i) & t <=(ts + t_i+tp)) tmp = ((-(dp*mu)/ts)*(t) + ((mu*dp)/ts)*(t_i+tp) + mu)
    else tmp = (-m*(t)+m*(t_i+tp+ts)+(1-dp)*mu)
    y[i] <- tmp + eff_y0
    if(y[i] < lower_bound) y[i] = lower_bound
    i <- i + 1
  }
  y
}

