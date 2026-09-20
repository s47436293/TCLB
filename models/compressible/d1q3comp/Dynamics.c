

//Lattice Sound Speed Constants for D1Q3 Lattice
#define CS2 (1.0/3.0)
#define CS4 (CS2*CS2)

// Boundary safe Field readign functions

CudaDeviceFunction int ClampOffset(int d) {
    switch (NodeType & NODE_BOUNDARY) {
    case NODE_WOutflow:
        if (d < 0) d = 0;
        break;
    case NODE_EOutflow:
        if (d > 0) d = 0;
        break;
    }
    switch (NodeType & NODE_BUFFER) {
    case NODE_WBuffer:
        if (d < -1) d = -1;
        break;
    case NODE_EBuffer:
        if (d > 1) d = 1;
        break;
    }
    return d;
}

CudaDeviceFunction real_t sAt(int d) {
    switch (ClampOffset(d)) {
        case -2: return s(-2,0);
        case -1: return s(-1,0);
        case  1: return s( 1,0);
        case  2: return s( 2,0);
        default: return s( 0,0);
    }
}

CudaDeviceFunction real_t rhoAt(int d) {
    switch (ClampOffset(d)) {
        case -1: return rho(-1,0);
        case  1: return rho( 1,0);
        default: return rho( 0,0);
    }
}

CudaDeviceFunction real_t uAt(int d) {
    switch (ClampOffset(d)) {
        case -1: return u(-1,0);
        case  1: return u( 1,0);
        default: return u( 0,0);
    }
}

CudaDeviceFunction real_t thetaAt(int d) {
    switch (ClampOffset(d)) {
        case -1: return theta(-1,0);
        case  1: return theta( 1,0);
        default: return theta( 0,0);
    }
}

//Thermodynmic Equation Helper Functions

CudaDeviceFunction real_t CV() { return CS2/(gamma_g - 1.0); }
CudaDeviceFunction real_t CP() { return gamma_g * CV(); }

CudaDeviceFunction real_t SFromState(real_t rho_, real_t theta_) {
    return CV() * log(theta_ / pow(rho_, gamma_g - 1.0));
}
CudaDeviceFunction real_t ThetaFromS(real_t rho_, real_t s_) {
    return pow(rho_, gamma_g - 1.0) * exp(s_/CV());
}

CudaDeviceFunction real_t CalcPressure(real_t rho_, real_t theta_) {
    return rho_ * CS2 * theta_;
}

//Correction Term Helper Functions

CudaDeviceFunction real_t A1FD_xx(real_t tau_t, real_t p_loc, real_t dudx) {
    return -tau_t * p_loc * (2.0 - 2.0/D_dim) * dudx;
}

CudaDeviceFunction real_t CalcE2xx(real_t rho_, real_t theta_, real_t dudx_) {
    real_t p_loc = CalcPressure(rho_, theta_);
    return p_loc * ( (D_dim + 2.0)/D_dim - gamma_g ) * dudx_;
}

CudaDeviceFunction real_t CalcE1xx(real_t G0, real_t Gp, real_t Gm, real_t sgn) {
    real_t dGb = G0 - Gm;
    real_t dGf = Gp - G0;
    return 0.5*(1.0 + sgn)*dGb + 0.5*(1.0 - sgn)*dGf;
}

CudaDeviceFunction real_t Gamma(real_t rho_, real_t theta_, real_t u_) {
    return rho_ * u_ * (1.0 - theta_ - u_*u_);
}

CudaDeviceFunction void Equilibrium(real_t rho_, real_t theta_, real_t u_, real_t fe_[3]) {
    real_t u2 = u_*u_;
    real_t A  = u2/(2.0*CS4) - u2/(2.0*CS2) + (theta_ - 1.0)*(1.0 - CS2)/(2.0*CS2);
    fe_[0] = (2.0/3.0) * rho_ * (1.0 - u2*0.5/CS2 - 0.5*(theta_ - 1.0));
    fe_[1] = (1.0/6.0) * rho_ * (1.0 + u_/CS2 + A);
    fe_[2] = (1.0/6.0) * rho_ * (1.0 - u_/CS2 + A);
}

// Common velocity gradient, used by A1FD, E2 and Phi. 
CudaDeviceFunction real_t DUDX() {
    return 0.5 * ( uAt(1) - uAt(-1) );
}

//Entropy Equation Update Helper Functions

CudaDeviceFunction real_t VanAlbada(real_t num, real_t den){
    if (fabs(den) < 1e-30) return 0.0;
    real_t r = num/den;
    if (r <= 0.0) return 0.0;
    return 2.0*r / (1.0 + r*r);
}

/* One MUSCL face value, Eq. (4.3.35a/b).
   sL  = s at the cell the slope is built from
   dm  = delta s one half-index below that cell
   dp  = delta s one half-index above that cell
   dir = +1 for a left-biased (L) state, -1 for a right-biased (R) state */
CudaDeviceFunction real_t cellInterfaceFunc(real_t sL, real_t dm, real_t dp, real_t dir){
    real_t phi = VanAlbada(dm, dp);
    if (dir > 0.0)
        return sL + dir*0.25*phi*((1.0-chi)*dm + (1.0+chi)*dp);
    else
        return sL + dir*0.25*phi*((1.0-chi)*dp + (1.0+chi)*dm);
}

CudaDeviceFunction real_t AdvectEntropy(real_t u_o){
    real_t sm2 = sAt(-2);
    real_t sm1 = sAt(-1);
    real_t s0  = sAt( 0);
    real_t sp1 = sAt( 1);
    real_t sp2 = sAt( 2);

    real_t d_3mh = sm1 - sm2;
    real_t d_mh  = s0  - sm1;
    real_t d_ph  = sp1 - s0;
    real_t d_3ph = sp2 - sp1;

    real_t SL_ph = cellInterfaceFunc(s0,  d_mh,  d_ph,   1.0);
    real_t SR_ph = cellInterfaceFunc(sp1, d_ph,  d_3ph, -1.0);
    real_t SL_mh = cellInterfaceFunc(sm1, d_3mh, d_mh,   1.0);
    real_t SR_mh = cellInterfaceFunc(s0,  d_mh,  d_ph,  -1.0);

    real_t S_ph = (u_o >= 0.0) ? SL_ph : SR_ph;
    real_t S_mh = (u_o >= 0.0) ? SL_mh : SR_mh;

    return u_o * (S_ph - S_mh);
}

CudaDeviceFunction real_t Fourier_term(real_t lam){
    return lam * ( thetaAt(1) - 2.0*thetaAt(0) + thetaAt(-1) );
}

CudaDeviceFunction real_t PhiTerm(real_t tau_, real_t tau_t_){
    return -(tau_/tau_t_) * a1xx(0,0) * DUDX();
}

//Initilisation Step


//Initial State from the XML inputs (Algorithm Step 1)

CudaDeviceFunction void SetMacro() {
    real_t theta_n = T_inf / T_r;
    rho   = InitRho;
    theta = theta_n;
    u     = InitU;
    s     = SFromState(InitRho, theta_n);
}

CudaDeviceFunction void CalcInitCollision() {
    const real_t Hproj[3] = {-1.0, 0.5, 0.5};

    real_t rho_n   = rhoAt(0);
    real_t u_n     = uAt(0);
    real_t theta_n = thetaAt(0);

    //Step 2 Equillibrium, written to storage
    Equilibrium(rho_n, theta_n, u_n, feq);

    real_t dudx  = DUDX();
    real_t p_loc = CalcPressure(rho_n, theta_n);
    real_t tau   = Viscosity / p_loc;
    real_t tau_t = tau + 0.5;
    real_t omega = 1.0 / tau_t;

    // step 3 - A1 from finite Differences
    a1xx = A1FD_xx(tau_t, p_loc, dudx);

    //Psi
    real_t E2xx = CalcE2xx(rho_n, theta_n, dudx);

    real_t Gm  = Gamma(rhoAt(-1), thetaAt(-1), uAt(-1));
    real_t G0  = Gamma(rho_n,     theta_n,     u_n    );
    real_t Gp  = Gamma(rhoAt( 1), thetaAt( 1), uAt( 1));
    real_t sgn = (real_t)((u_n > 0.0) - (u_n < 0.0));

    real_t E = CalcE1xx(G0, Gp, Gm, sgn) + E2xx;

    Exx = E;

    //Step 5 just for first iteration so can transition into iteration
    for (int k = 0; k < 3; k++)
        g[k] = feq[k] + Hproj[k] * ( (1.0 - omega)*a1xx + 0.5*E );
}

//Boundary Conditions
/* Zero-gradient outflow.  gg holds the streamed populations of this node;
   the one that arrived from outside the domain is replaced with the
   node's own equilibrium, taken straight from the stored feq.
      g[1] moves in +x, pulled from x-1: unknown at the WEST edge.
      g[2] moves in -x, pulled from x+1: unknown at the EAST edge. */

CudaDeviceFunction void WOutflow(real_t gg[3], const real_t feq_[3]) {
    gg[1] = feq_[1];
}

CudaDeviceFunction void EOutflow(real_t gg[3], const real_t feq_[3]) {
    gg[2] = feq_[2];
}


//Iteration Stage 1 Steps 5 and 7 of algorithimn 

CudaDeviceFunction void CalcRelaxation() {
    const real_t Hproj[3] = {-1.0, 0.5, 0.5};

    //level n values
    real_t rho_o   = rhoAt(0);
    real_t u_o     = uAt(0);
    real_t theta_o = thetaAt(0);
    real_t s_o     = sAt(0);

    real_t mu    = Viscosity;
    real_t p_loc = CalcPressure(rho_o, theta_o);
    real_t tau   = mu / p_loc;
    real_t tau_t = tau + 0.5;
    real_t omega = 1.0 / tau_t;
    real_t lam   = mu * CP() / Pr;

    //Entropy, step 7 of algorithm
    real_t inv_rhotheta = 1.0 / (rho_o * theta_o);

    real_t s_n = s_o
               - AdvectEntropy(u_o)
               + inv_rhotheta * Fourier_term(lam)
               + inv_rhotheta * PhiTerm(tau, tau_t);

    //Populations, step 5 of algorithm, Eq. (4.3.23)
    real_t E = Exx(0,0);
    for (int k = 0; k < 3; k++)
        g[k] = feq[k] + Hproj[k] * ( (1.0 - omega)*a1xx(0,0) + 0.5*E );

    s = s_n;
}

//Stage 2 step 6 of algorithm

CudaDeviceFunction void CalcMoments() {
    //Boundary condition fixing population coming from outside domain.
    //feq is still at level n here - PsiEQ has not run yet this iteration.
    real_t gg[3];
    gg[0] = g[0];
    gg[1] = g[1];
    gg[2] = g[2];

    switch (NodeType & NODE_BOUNDARY) {
    case NODE_WOutflow:
        WOutflow(gg, feq);
        break;
    case NODE_EOutflow:
        EOutflow(gg, feq);
        break;
    }

    real_t rho_n = gg[0] + gg[1] + gg[2];

    rho = rho_n;
    u   = ( gg[1] - gg[2] ) / rho_n;
}

//Stage 3 Step 8 of algorithm

CudaDeviceFunction void CalcTheta() {
    theta = ThetaFromS(rhoAt(0), sAt(0));
}

//Stage 4 Step 9 and 10 of algorithm

CudaDeviceFunction void CalcEqnPsi() {
    real_t rho_n   = rhoAt(0);
    real_t u_n     = uAt(0);
    real_t theta_n = thetaAt(0);

    real_t E2xx = CalcE2xx(rho_n, theta_n, DUDX());

    real_t Gm  = Gamma(rhoAt(-1), thetaAt(-1), uAt(-1));
    real_t G0  = Gamma(rho_n,     theta_n,     u_n    );
    real_t Gp  = Gamma(rhoAt( 1), thetaAt( 1), uAt( 1));
    real_t sgn = (real_t)((u_n > 0.0) - (u_n < 0.0));

    Exx = CalcE1xx(G0, Gp, Gm, sgn) + E2xx;   //step 9

    Equilibrium(rho_n, theta_n, u_n, feq);    //step 10, straight into storage
}

//Stage 5 step 11 of algorithm

CudaDeviceFunction void CalcA1() {
    const real_t Hproj[3] = {-1.0, 0.5, 0.5};
    const real_t Hxx[3]   = {-CS2, 1.0 - CS2, 1.0 - CS2};

    real_t rho_n   = rhoAt(0);
    real_t u_n     = uAt(0);
    real_t theta_n = thetaAt(0);

    real_t dudx  = DUDX();
    real_t p_loc = CalcPressure(rho_n, theta_n);
    real_t tau_t = Viscosity / p_loc + 0.5;

    real_t gg[3];
    gg[0] = g[0];
    gg[1] = g[1];
    gg[2] = g[2];

    switch (NodeType & NODE_BOUNDARY) {
    case NODE_WOutflow:
        WOutflow(gg, feq);
        break;
    case NODE_EOutflow:
        EOutflow(gg, feq);
        break;
    }

    real_t E = Exx(0,0);

    real_t a1_pr = gg[1] + gg[2] - CS2*( gg[0] + gg[1] + gg[2] );
    for (int k = 0; k < 3; k++)
        a1_pr -= Hxx[k] * ( feq[k] - 0.5*Hproj[k]*E );

    real_t a1_fd = A1FD_xx(tau_t, p_loc, dudx);

    a1xx = sigma*a1_pr + (1.0 - sigma)*a1_fd;
}

//Outputs

CudaDeviceFunction real_t getRho()      { return rho(0,0); }
CudaDeviceFunction real_t getTheta()    { return theta(0,0); }
CudaDeviceFunction real_t getEntropy()  { return s(0,0); }
CudaDeviceFunction real_t getPressure() { return CalcPressure(rho(0,0), theta(0,0)); }
CudaDeviceFunction real_t getA1xx()     { return a1xx(0,0); }

CudaDeviceFunction vector_t getU() {
    vector_t ret;
    ret.x = u(0,0);
    ret.y = 0.0;
    ret.z = 0.0;
    return ret;
}

CudaDeviceFunction float2 Color() {
    float2 ret;
    ret.x = 0;
    ret.y = 1;
    return ret;
}