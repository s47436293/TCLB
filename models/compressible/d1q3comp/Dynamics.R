

AddDensity( name="g[0]", dx= 0, dy= 0, group="g")
AddDensity( name="g[1]", dx= 1, dy= 0, group="g")
AddDensity( name="g[2]", dx=-1, dy= 0, group="g")



#Fields

AddField(name="rho",   dx=c(-1,1), dy=c(0,0), group="macro",   comment="Density")
AddField(name="u",     dx=c(-1,1), dy=c(0,0), group="macro",   comment="Particle velocity")
AddField(name="theta", dx=c(-1,1), dy=c(0,0), group="macro",   comment="Lattice temperature")
AddField(name="s",     dx=c(-2,2), dy=c(0,0), group="macro",   comment="Entropy")
AddField(name="pixx",  dx=c( 0,0), dy=c(0,0), group="moments", comment="Second Hermite moment of g, cached from Macro")
AddField(name="a1xx",  dx=c( 0,0), dy=c(0,0), group="tensor",  comment="Second-order off-equilibrium moment")

#stages

AddStage("Macro", "CalcMacro",
    load = DensityAll$group %in% c("g", "macro", "tensor"),
    save = Fields$group %in% c("macro", "moments"))

AddStage("Collide", "CalcCollision",
    load = DensityAll$group %in% c("macro", "moments"),
    save = Fields$group %in% c("g", "tensor"))

AddStage("InitMacro", "SetMacro",
    save = Fields$group %in% c("macro", "moments"))

AddStage("InitCollide", "CalcInitCollision",
    load = DensityAll$group == "macro",
    save = Fields$group %in% c("g", "tensor"))

AddAction("Init",      c("InitMacro", "InitCollide"))
AddAction("Iteration", c("Macro", "Collide"))



#Physical properties
AddSetting(name="Viscosity", default=0.0, unit="kg/m/s", comment="dynamic viscosity mu")
AddSetting(name="gamma_g",   default=1.4,  comment="heat capacity ratio")
AddSetting(name="Pr",        default=0.71, comment="Prandtl number")
AddSetting(name="T_r", default = 1.0, unit = "K", comment ="Reference temperature, lattice theta = T/T_r" )

#Scheme
AddSetting(name="sigma",   default=0.4, comment="a1 blend: sigma*a1_pr + (1-sigma)*a1_fd")
AddSetting(name="D_dim",   default=1.0, comment="spatial dimension in kappa=(D+2)/D")
AddSetting(name="chi",     default=0.33333333333333,  comment="MUSCL parameter")

#Initial condition - zonal: set per named zone painted in the XML
AddSetting(name="InitRho",   default=1.0, zonal=T, unit="kg/m3", comment="initial density")
AddSetting(name="T_inf", default=1.0, zonal=T, unit="K",     comment="initial temperature, theta = T_inf/T_r")
AddSetting(name="InitU",     default=0.0, zonal=T, unit="m/s",   comment="initial velocity")



#Boundary Condition Nodes
AddNodeType(name="WOutflow",   group="BOUNDARY")
AddNodeType(name="EOutflow",   group="BOUNDARY")
AddNodeType(name="WBuffer",    group="BUFFER")
AddNodeType(name="EBuffer",    group="BUFFER")

#Outputs

AddQuantity(name="Rho")
AddQuantity(name="U", vector=TRUE)
AddQuantity(name="Theta")
AddQuantity(name="Entropy")
AddQuantity(name="Pressure")
AddQuantity(name="A1xx")