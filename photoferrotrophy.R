
# =================================== #
# anoxygenic photosynthesis model 1.0 #
# =================================== #

# Chris Reinhard
# Model for exploring competition between 
# oxygenic and anoxygenic photosynthesizers 
# in a 1-D water column
#
# first published in: Ozaki et al., 2019 [doi:10.1038/s41467-019-10872-z]
#

# ======== load required packages ========== #
require(ReacTran)
require(deSolve)
require(rootSolve)
require(marelac)
require(tictoc)

tic("runtime")
print("running...")

# ========================================== #
# === physical/biogeochemical parameters === #
# ========================================== #

# --- boundary conditions --- #
t           <- 25                  # seawater temp (ºC)
PO4_top     <- 0.0                 # [PO4], top boundary (µmol/kg)  
PO4_bottom  <- 1.0                 # [PO4], bottom boundary (µmol/kg)  
Fe_top      <- 0.0                 # [Fe2+], top boundary (µmol/kg)   
rFeP        <- 300                 # deepwater Fe/P ratio
Fe_bottom   <- PO4_bottom*rFeP     # [Fe2+], bottom boundary (µmol/kg)  
O2_top      <- 1.0                 # [O2], top boundary (µmol/kg)
O2_bottom   <- 0.0                 # [O2], bottom boundary (µmol/kg)

# --- transport parameters --- #
Kv          <- 1e-4                # vertical turbulent diffusivity (m2/s)
w           <- - ( 0.5/(86400) )   # upwelling velocity (m/s)

kO2         <- gas_transfer(t=t, u10=7, species=c("O2"), 
                            method=c("Nightingale"), 
                            Schmidt=gas_schmidt(t=t, species=c("O2")))   # piston velocity for O2 gas exchange (m/s)

# --- light penetration/limitation parameters --- #
I0          <- 1300                # incident radiation at the surface (µEinst/m2/s)
k_I         <- 15                  # attenuation length scale (m)
Km_I_ox     <- 98                  # light-limited Km, oxic photosynthesizers (µEinst/m2/s)
Km_I_anox   <- 1                   # light-limited Km, photoferrotrophs (µEinst/m2/s)

# --- nutrient growth parameters --- #
Vmax_ox         <- 0.004           # max growth rate, oxygenic photosynthesizers (µmolC/kg/s)
Vmax_anox       <- 0.002           # max growth rate, photoferrotrophs (µmolC/kg/s)
Km_P_ox         <- 0.015           # phosphorus Km, oxygenic photosynthesizers (µmol/kg)
Km_P_anox       <- 0.005           # phosphorus Km, photoferrotrophs (µmol/kg)
Km_Fe_anox      <- 11              # iron Km, photoferrotrophs (µmol/kg)
rCP_ox          <- 106             # C/P ratio, oxygenic photosynthesizers
rCP_anox        <- 106             # C/P ratio, photoferrotrophs
rFeC            <- 4               # Fe/C ratio, photoferrotrophs

# --- iron oxidation parameters --- #
temp         <- t + 273.15                          # seawater temp (K)
I            <- 0.7                                  # seawater ionic strength
k0           <- 21.56 - (1545/temp)                  # temp-corrected rate constant
kFe          <- ( (10^(k0-(3.29*(I^0.5))+(1.52*I)))
                 / 60) / (1e6)                       # combined rate constant (kg3/µmol3/s)
pH           <- 8                                    # seawater pH
OH           <- (10^(-(14-pH))) * (1e6)              # seawater [OH-] (µmol/kg)
nO2          <- 1                                    # reaction order w.r.t. O2
nOH          <- 2                                    # reaction order w.r.t. OH-

# --- P scavenging parameters --- #
K_d          <- 0.025             # distribution coefficient for P on Fe(III)

  
# =============================================== #
# === set up model domain and grid properties === #
# =============================================== #

L          <- 500                                           # size of water column (m)
N          <- 100000                                        # number of grid layers
grid       <- setup.grid.1D(x.up=0, x.down=N, L=L, N=N)     # sets up finite element grid
light.grid <- setup.prop.1D(func=p.exp, grid=grid, 
                            y.0=I0, y.inf=0, x.att=k_I)     # sets up light penetration grid


# ========================================================== #
# === function defining system of differential equations === #
# ========================================================== #

anox_photo    <- function(t, Conc, parms = NULL)
{
  # --- pass initial concentrations --- #
  PO4  <-  Conc[1:N]
  Fe   <-  Conc[(N+1):(2*N)]
  O2   <-  Conc[(2*N+1):(3*N)]
  
  # --- diffusion/advection functions --- #
  tranPO4   <- tran.1D(C=PO4, C.up=NULL, C.down=PO4_bottom,
                       flux.up=NULL, flux.down=NULL, D=Kv, v=w, dx=grid)$dC
  
  tranFe    <- tran.1D(C=Fe, C.up=NULL, C.down=Fe_bottom,
                      flux.up=NULL, flux.down=NULL, D=Kv, v=w, dx=grid)$dC
  
  tranO2    <- tran.1D(C=O2, C.up=O2_top, C.down=NULL,
                       flux.up=-(kO2*(O2-O2_top)), flux.down=NULL, 
                       a.bl.up=NULL, a.bl.down=NULL, 
                       D=Kv, v=w, dx=grid)$dC

  
  # --- reaction terms --- #
  ox_photo   <- ( Vmax_ox  
                 * (light.grid$mid/(light.grid$mid+Km_I_ox)) 
                 * (PO4/(PO4+Km_P_ox)) )                             # rate of oxygenic photosynthesis
  anox_photo <- ( Vmax_anox 
                 * (light.grid$mid/(light.grid$mid+Km_I_anox)) 
                 * (Fe/(Fe+Km_Fe_anox)) * (PO4/(PO4+Km_P_anox)) )    # rate of photoferrotrophy
  Fe_ox      <- kFe * Fe * (O2^nO2) * (OH^nOH)                       # Fe(II) oxidation rate
  P_scav     <- Fe_ox * (K_d*PO4)                                    # equilibrium P scavenging term
  
  # --- differential equations --- #
  dPO4  <- tranPO4 - (1/rCP_ox)*ox_photo - (1/rCP_anox)*anox_photo - P_scav
  dFe   <- tranFe  - rFeC*anox_photo     - Fe_ox
  dO2   <- tranO2  + ox_photo            - (0.25*Fe_ox)
  
  # --- other outputs --- #
  return(list(c(dPO4=dPO4,dFe=dFe,dO2=dO2),Fe_ox=Fe_ox,ox=ox_photo,anox=anox_photo))    
}


# ======================= #
# === solve equations === #
# ======================= #

# initial guess at state variables
PO4        <- rep(PO4_bottom,N)
Fe         <- rep(Fe_bottom,N)
O2         <- rep(O2_bottom,N)

Conc       <- c(PO4,Fe,O2)

# --- run model --- #
sol <- steady.1D(y=Conc, func=anox_photo, parms=NULL, nspec=3, positive=TRUE, atol=1e-10)

# ======================= #
# === post-processing === #
# ======================= #

# --- collect/calculate some results --- #
depth      <- grid$x.mid
light      <- light.grid$mid
PO4_t      <- sol$y[1:N]
Fe_t       <- sol$y[(N+1) : (2*N)]
O2_t       <- sol$y[(2*N+1) : (3*N)]
FeP        <- Fe_t/PO4_t
ox         <- ( Vmax_ox 
               * (light.grid$mid/(light.grid$mid+Km_I_ox))
               *((PO4_t/(PO4_t+Km_P_ox))) ) * (86400)         # rates of oxygenic photosynthesis (µmolC/kg/d)
anox       <- ( Vmax_anox 
               * (light.grid$mid/(light.grid$mid+Km_I_anox))
               *(Fe_t/(Fe_t+Km_Fe_anox))
               *(PO4_t/(PO4_t+Km_P_anox)) ) * (86400)            # rates of photoferrotrophy (µmolC/kg/d)
max_ox     <- max(ox)                                            # peak rate of oxygenic photosynthesis
max_anox   <- max(anox)                                          # peak rate of photoferrotrophy
Fe_ox      <- ( kFe * Fe_t * (O2_t^nO2) * (OH^nOH) ) * (86400)   # rate of Fe(II) oxidation (µmol/kg/d)

# --- function for integrating photosynthesis curves --- #
auc <-
  function(x, y, from = min(x), to = max(x), type=c("linear", "spline"), ...) 
  {
    type <- match.arg(type)
    
    if (length(x) != length(y)) 
      stop("x and y must have the same length")
    if (length(unique(x)) < 2) 
      return(NA)
    
    if (type=="linear") {    
      values <- approx(x, y, xout = sort(unique(c(from, to, x[x > from & x < to]))), ...)
      res <- 0.5 * sum(diff(values$x) * (values$y[-1] + values$y[-length(values$y)]))
    } else {
      res <- integrate(splinefun(x, y, method="natural"), lower=from, upper=to)$value
    }
    res
  }

# --- calculate column-integrated rates of photosynthesis
tot_anox  <- auc(grid$x.mid, anox)
tot_ox    <- auc(grid$x.mid, ox)
f_anox    <- tot_anox/(tot_anox+tot_ox)
f_ox      <- tot_ox/(tot_anox+tot_ox)

# --- calculate column-integrated rate of Fe oxidation
tot_Fe_ox <- auc(grid$x.mid, Fe_ox)


# --- compile and dump some output --- #
out_df<-as.data.frame(cbind(depth,light,PO4_t,Fe_t,O2_t,FeP,ox,anox,Fe_ox))
write.csv(out_df, file = "~/_output/photosynth.out.csv")

# --- plot results --- #
par(mfrow=c(1,4))

matplot(light, depth, xlab="PAR (µEinst/m2/s)", ylab="depth [m]", type="l", ylim=c(L,0), col="blue", lty=1, lwd=3)
matplot(PO4_t, depth, xlab="[PO4] (µmol/kg)", ylab="depth [m]", type="l", xlim=c(0,max(PO4_t)), ylim=c(L,0), col="black", lty=1, lwd=3)
matplot(Fe_t, depth, xlab="[Fe(II)] (µmol/kg)", ylab="depth [m]", type="l", xlim=c(0,max(Fe_t)), ylim=c(L,0), col="tan", lty=1, lwd=3)
matplot(ox, depth, xlab="photosynthesis", ylab="depth [m]", type="l", xlim=c(0,max(max(ox),max(anox))), ylim=c(L,0), col="darkgreen", lty=1, lwd=3)
lines(anox, depth, xlab="photosynthesis", ylab="depth [m]", type="l", xlim=c(0,max(max(ox),max(anox))), ylim=c(L,0), col="darkred", lty=1, lwd=3)


# --- some diagnostics --- #
attr(sol,"precis")    # quotes solution precision
attr(sol,"steady")    # tracks whether inversion achieves steady state (TRUE/FALSE)
tot_ox                # column-integrated oxygenic photosynthesis
tot_anox              # column-integrated photoferrotrophy
f_ox                  # oxygenic fraction of total photosynthesis 
f_anox                # photoferrotrophic fraction of total photosynthesis

print("...finishing")
toc()

