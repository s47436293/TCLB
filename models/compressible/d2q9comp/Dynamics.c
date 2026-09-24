

//Lattice Sound Speed Constants for D2Q9 Lattice
#define CS2 (1.0/3.0)
#define CS4 (CS2*CS2)

// Boundary safe Field reading functions

CudaDeviceFunction int ClampY(int d) {
    switch (NodeType & NODE_YBOUNDARY) {
    case NODE_SOutflow: if (d <  0) d =  0; break;
    case NODE_NOutflow: if (d >  0) d =  0; break;
    case NODE_SBuffer:  if (d < -1) d = -1; break;
    case NODE_NBuffer:  if (d >  1) d =  1; break;
    case NODE_Wall: d = 0; break;
    }
    return d;
}
CudaDeviceFunction int ClampX(int d) {
    switch (NodeType & NODE_XBOUNDARY) {
    case NODE_WOutflow: if (d <  0) d =  0; break;
    case NODE_EOutflow: if (d >  0) d =  0; break;
    case NODE_WBuffer:  if (d < -1) d = -1; break;
    case NODE_EBuffer:  if (d >  1) d =  1; break;
    }
    return d;
}
CudaDeviceFunction real_t uxAtY(int d) {   // ux at y-offsets, ClampY
    switch (ClampY(d)) {
        case -1: return ux(0,-1);
        case  1: return ux(0, 1);
        default: return ux(0, 0);
    }
}
CudaDeviceFunction real_t uyAtX(int d) {   // uy at x-offsets, ClampX
    switch (ClampX(d)) {
        case -1: return uy(-1,0);
        case  1: return uy( 1,0);
        default: return uy( 0,0);
    }
}

/* Anything that streamed into a Wall node is sent straight back
   the way it came.   1<->3  2<->4  5<->7  6<->8                    */
CudaDeviceFunction void BounceBack() {
    real_t t;
    t = g[1]; g[1] = g[3]; g[3] = t;
    t = g[2]; g[2] = g[4]; g[4] = t;
    t = g[5]; g[5] = g[7]; g[7] = t;
    t = g[6]; g[6] = g[8]; g[8] = t;
}
CudaDeviceFunction real_t sX(int d) {
    switch (ClampX(d)) {
        case -2: return s(-2,0);
        case -1: return s(-1,0);
        case  1: return s( 1,0);
        case  2: return s( 2,0);
        default: return s( 0,0);
    }
}

CudaDeviceFunction real_t sY(int d) {
    switch (ClampY(d)) {
        case -2: return s(0,-2);
        case -1: return s(0,-1);
        case  1: return s(0, 1);
        case  2: return s(0, 2);
        default: return s(0, 0);
    }
}

CudaDeviceFunction real_t rhoAtX(int d) {
    switch (ClampX(d)) {
        case -1: return rho(-1,0);
        case  1: return rho( 1,0);
        default: return rho( 0,0);
    }
}


CudaDeviceFunction real_t rhoAtY(int d) {
    switch (ClampY(d)) {
        case -1: return rho(0,-1);
        case  1: return rho(0,1);
        default: return rho( 0,0);
    }
}

CudaDeviceFunction real_t uxAt(int d) {
    switch (ClampX(d)) {
        case -1: return ux(-1,0);
        case  1: return ux( 1,0);
        default: return ux( 0,0);
    }
}
CudaDeviceFunction real_t uyAt(int d) {
    switch (ClampY(d)) {
        case -1: return uy(0,-1);
        case  1: return uy( 0,1);
        default: return uy( 0,0);
    }
}

CudaDeviceFunction real_t thetaAtX(int d) {
    switch (ClampX(d)) {
        case -1: return theta(-1,0);
        case  1: return theta( 1,0);
        default: return theta( 0,0);
    }
}

CudaDeviceFunction real_t thetaAtY(int d) {
    switch (ClampY(d)) {
        case -1: return theta(0,-1);
        case  1: return theta(0,1);
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


//Relaxation Helping Functions
CudaDeviceFunction void Equilibrium(real_t rho_, real_t theta_, real_t ux_, real_t uy_, real_t feq_[9]) {
    const real_t w[9]  = {4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.};
    const real_t cx[9] = {0, 1, 0, -1,  0, 1, -1, -1,  1};
    const real_t cy[9] = {0, 0, 1,  0, -1, 1,  1, -1, -1};
    real_t tc  = (theta_ - 1.0) * CS2;
    real_t axx = ux_*ux_ + tc;
    real_t ayy = uy_*uy_ + tc;
    real_t axy = ux_*uy_;
    real_t axxy  = (ux_*ux_ + tc) * uy_;                              /* NEW */
    real_t ayyx  = (uy_*uy_ + tc) * ux_;                              /* NEW */
    real_t axxyy = ux_*ux_*uy_*uy_ + tc*(ux_*ux_ + uy_*uy_) + tc*tc;  /* NEW */

    for (int k = 0; k < 9; k++) {
        real_t Hxx = cx[k]*cx[k] - CS2;
        real_t Hyy = cy[k]*cy[k] - CS2;
        real_t Hxy = cx[k]*cy[k];
        feq_[k] = w[k] * rho_ * ( 1.0
               + ( cx[k]*ux_ + cy[k]*uy_ ) / CS2
               + ( Hxx*axx + 2.0*Hxy*axy + Hyy*ayy ) / (2.0*CS4)
               + ( Hxx*cy[k]*axxy + Hyy*cx[k]*ayyx ) / (2.0*CS2*CS4)   /* NEW */
               + ( Hxx*Hyy*axxyy ) / (4.0*CS4*CS4) );                  /* NEW */
    }
}

CudaDeviceFunction void OffEquilibrium(real_t ux_, real_t uy_, real_t theta_, real_t axx, real_t axy, real_t ayy, real_t g1_[9]) {
    const real_t w[9]  = {4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.};
    const real_t cx[9] = {0, 1, 0, -1,  0, 1, -1, -1,  1};
    const real_t cy[9] = {0, 0, 1,  0, -1, 1,  1, -1, -1};
    real_t tc = CS2 * (theta_ - 1.0);

    real_t axxy  = uy_*axx + 2.0*ux_*axy;
    real_t ayyx  = ux_*ayy + 2.0*uy_*axy;
    real_t axxyy = 2.0*( ux_*ayyx + uy_*axxy )
                 + ( tc - ux_*ux_ )*ayy
                 + ( tc - uy_*uy_ )*axx
                 - 4.0*ux_*uy_*axy;

    for (int k = 0; k < 9; k++){
        real_t Hxx = cx[k]*cx[k] - CS2;
        real_t Hyy = cy[k]*cy[k] - CS2;
        real_t Hxy = cx[k]*cy[k];
        real_t Hxxy = cx[k]*cx[k] * cy[k] - CS2 * cy[k];
        real_t Hyyx = cy[k]*cy[k] * cx[k] - CS2 * cx[k];
        real_t Hxxyy = pow(cx[k], 2) * pow(cy[k], 2) - CS2 * (pow(cx[k], 2)  +  pow(cy[k],2)) + CS4;
        g1_[k] = w[k] * ((axx * Hxx + 2* axy * Hxy + ayy * Hyy)/(2 * CS4) + (axxy * Hxxy + ayyx * Hyyx) / ((2*CS4*CS2)) + axxyy*Hxxyy/(4*CS4 * CS4)) ;
    }
}


//Correction Term Helper Functions

CudaDeviceFunction real_t A1PR_xx(real_t gi_[9],real_t feq_[9], real_t psi_[9]){
    const real_t cx[9] = {0, 1, 0, -1,  0, 1, -1, -1,  1};
    real_t A1PRxx = 0.0;
    for (int k = 0; k < 9; k++) {
        real_t Hxx = cx[k]*cx[k] - CS2;
        A1PRxx += Hxx * (gi_[k]- feq_[k] + 0.5 * psi_[k]);
    }
    return A1PRxx;
}

CudaDeviceFunction real_t A1PR_xy(real_t gi_[9], real_t feq_[9], real_t psi_[9]){
    const real_t cx[9] = {0, 1, 0, -1,  0, 1, -1, -1,  1};
    const real_t cy[9] = {0, 0, 1,  0, -1, 1,  1, -1, -1};
    real_t A1PRxy = 0.0;
    for (int k = 0; k < 9; k++) {
        real_t Hxy = cx[k]*cy[k];
        A1PRxy += Hxy * (gi_[k]- feq_[k] + 0.5 * psi_[k]);
    }
    return A1PRxy;
}

CudaDeviceFunction real_t A1PR_yy(real_t gi_[9], real_t feq_[9], real_t psi_[9]){
    const real_t cy[9] = {0, 0, 1,  0, -1, 1,  1, -1, -1};
    real_t A1PRyy = 0.0;
    for (int k = 0; k < 9; k++) {
        real_t Hyy = cy[k]*cy[k] - CS2;
        A1PRyy += Hyy * (gi_[k]- feq_[k] + 0.5 * psi_[k]);
    }
    return A1PRyy;
}

CudaDeviceFunction real_t A1FD_aa(real_t tau_t, real_t p_loc, real_t FD1_, real_t FD2_) {
    return -tau_t * p_loc * (FD1_ - FD2_);
}

CudaDeviceFunction real_t A1FD_xy(real_t tau_t, real_t p_loc, real_t FD1_, real_t FD2_) {
    return -tau_t * p_loc * (FD1_ + FD2_);
}




//Calculating Psi Helper Functions

CudaDeviceFunction real_t Gamma(real_t rho_, real_t theta_, real_t u_) {
    return rho_ * u_ * (1.0 - theta_ - u_*u_);
}

CudaDeviceFunction real_t CalcE1ab(real_t G0, real_t Gp, real_t Gm, real_t sgn) {
    real_t dGb = G0 - Gm;
    real_t dGf = Gp - G0;
    return 0.5*(1.0 + sgn)*dGb + 0.5*(1.0 - sgn)*dGf;
}

CudaDeviceFunction real_t CalcE2ab(real_t rho_, real_t theta_, real_t dudx_, real_t dudy_) {
    real_t p_loc = CalcPressure(rho_, theta_);
    return p_loc * ( (D_dim + 2.0)/D_dim - gamma_g ) * (dudx_ + dudy_);
}


CudaDeviceFunction void CalcPsi(real_t psi_[9], real_t Exx, real_t Eyy){
    const real_t w[9]  = {4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.};
    const real_t cx[9] = {0, 1, 0, -1,  0, 1, -1, -1,  1};
    const real_t cy[9] = {0, 0, 1,  0, -1, 1,  1, -1, -1};
    for (int k = 0; k < 9; k++){
        real_t Hxx = cx[k]*cx[k] - CS2;
        real_t Hyy = cy[k]*cy[k] - CS2;
        psi_[k] = w[k]*(Hxx * Exx + Hyy * Eyy)/ (2 * CS4);
    }
}




// Common velocity gradient, used by A1FD, E2 and Phi. 
CudaDeviceFunction real_t DUXDX() { return 0.5 * ( uxAt(1)  - uxAt(-1)  ); }
CudaDeviceFunction real_t DUXDY() { return 0.5 * ( uxAtY(1) - uxAtY(-1) ); }
CudaDeviceFunction real_t DUYDX() { return 0.5 * ( uyAtX(1) - uyAtX(-1) ); }
CudaDeviceFunction real_t DUYDY() { return 0.5 * ( uyAt(1)  - uyAt(-1)  ); }


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

CudaDeviceFunction real_t AdvectLine(real_t sm2, real_t sm1, real_t s0, real_t sp1, real_t sp2, real_t un){
    real_t d_3mh = sm1 - sm2;
    real_t d_mh  = s0  - sm1;
    real_t d_ph  = sp1 - s0;
    real_t d_3ph = sp2 - sp1;

    real_t SL_ph = cellInterfaceFunc(s0,  d_mh,  d_ph,   1.0);
    real_t SR_ph = cellInterfaceFunc(sp1, d_ph,  d_3ph, -1.0);
    real_t SL_mh = cellInterfaceFunc(sm1, d_3mh, d_mh,   1.0);
    real_t SR_mh = cellInterfaceFunc(s0,  d_mh,  d_ph,  -1.0);

    real_t S_ph = (un >= 0.0) ? SL_ph : SR_ph;
    real_t S_mh = (un >= 0.0) ? SL_mh : SR_mh;

    return un * (S_ph - S_mh);
}

CudaDeviceFunction real_t AdvectEntropy(real_t ux_o, real_t uy_o){
    real_t s0 = s(0,0);
    return AdvectLine(sX(-2), sX(-1), s0, sX(1), sX(2), ux_o)
         + AdvectLine(sY(-2), sY(-1), s0, sY(1), sY(2), uy_o);
}

CudaDeviceFunction real_t Fourier_term(real_t lam){
    real_t t0 = theta(0,0);
    return lam * ( thetaAtX(1) - 2.0*t0 + thetaAtX(-1)
                 + thetaAtY(1) - 2.0*t0 + thetaAtY  (-1) );
}

CudaDeviceFunction real_t PhiTerm(real_t tau_, real_t tau_t_){
    return -(tau_/tau_t_) *   (a1xx(0,0) * DUXDX() + a1xy(0,0) * (DUXDY() + DUYDX()) + a1yy(0,0) * DUYDY());
}





//Initilisation Step

//Initial State from the XML inputs (Algorithm Step 1)

CudaDeviceFunction void SetMacro() {
    real_t theta_n = T_inf / T_r;
    rho   = InitRho;
    theta = theta_n;
    ux    = InitUx;
    uy    = InitUy;
    s     = SFromState(InitRho, theta_n);
}

CudaDeviceFunction void CalcInitCollision() {

    real_t rho_n   = rho(0,0);
    real_t ux_n = ux(0,0);
    real_t uy_n = uy(0,0);
    real_t theta_n = theta(0,0);

    //Step 2 Calculate Equillibrium 
    Equilibrium(rho_n, theta_n, ux_n, uy_n, feq);


    // step 3 - A1 from finite Differences
    real_t duxdx  = DUXDX();
    real_t duydy  = DUYDY();
    real_t duxdy  = DUXDY();
    real_t duydx  = DUYDX();


    real_t p_loc = CalcPressure(rho_n, theta_n);
    real_t tau   = Viscosity / p_loc;
    real_t tau_t = tau + 0.5;
    real_t omega = 1.0 / tau_t;


    real_t a1xx_n = A1FD_aa(tau_t, p_loc, duxdx, duydy);
    real_t a1yy_n = A1FD_aa(tau_t, p_loc, duydy, duxdx);
    real_t a1xy_n = A1FD_xy(tau_t, p_loc, duxdy, duydx);
    a1xx = a1xx_n;  a1yy = a1yy_n;  a1xy = a1xy_n;


    //Psi
    real_t E2ab= CalcE2ab(rho_n, theta_n, DUXDX(), DUYDY());

    //calculate E1xx First
    real_t GmX = Gamma(rhoAtX(-1), thetaAtX(-1),uxAt(-1));
    real_t G0X  = Gamma(rho_n,     theta_n,  ux_n);
    real_t GpX  = Gamma(rhoAtX(1), thetaAtX(1),uxAt(1));
    real_t sgnX = (real_t)((ux_n > 0.0) - (ux_n < 0.0));

    real_t E1xx = CalcE1ab(G0X, GpX, GmX, sgnX);    //step 9

    //Repeat for E1yy
    real_t GmY = Gamma(rhoAtY(-1), thetaAtY(-1),uyAt(-1));
    real_t G0Y  = Gamma(rho_n,     theta_n,  uy_n);
    real_t GpY  = Gamma(rhoAtY(1), thetaAtY(1),uyAt(1));
    real_t sgnY = (real_t)((uy_n > 0.0) - (uy_n < 0.0));

    real_t E1yy= CalcE1ab(G0Y, GpY, GmY, sgnY);    //step 9

    Exx = E1xx + E2ab;
    Eyy = E1yy + E2ab;

    real_t g1[9], psi[9];
    OffEquilibrium(ux_n, uy_n, theta_n, a1xx_n, a1xy_n, a1yy_n, g1);
    CalcPsi(psi, E1xx + E2ab, E1yy + E2ab);
    for (int k = 0; k < 9; k++)
        g[k] = feq[k] + (1.0 - omega)*g1[k] + 0.5*psi[k];

}

//Boundary Conditions
/* Zero-gradient outflow.  gg holds the streamed populations of this node;
   the one that arrived from outside the domain is replaced with the
   node's own equilibrium, taken straight from the stored feq.
      g[1] moves in +x, pulled from x-1: unknown at the WEST edge.
      g[2] moves in -x, pulled from x+1: unknown at the EAST edge. */


CudaDeviceFunction void NOutflow(real_t gg[9], const real_t feq_[9]) {
 gg[4] = feq_[4];
 gg[7] = feq_[7];
 gg[8] = feq_[8];
}

CudaDeviceFunction void EOutflow(real_t gg[9], const real_t feq_[9]) {
    gg[3] = feq_[3];
    gg[6] = feq_[6];
    gg[7] = feq_[7];
}

CudaDeviceFunction void SOutflow(real_t gg[9], const real_t feq_[9]) {
 gg[2] = feq_[2];
 gg[5] = feq_[5];
 gg[6] = feq_[6];
}
CudaDeviceFunction void WOutflow(real_t gg[9], const real_t feq_[9]) {
    gg[1] = feq_[1];
    gg[5] = feq_[5];
    gg[8] = feq_[8];
}




//Iteration Stage 1 Steps 5 and 7 of algorithimn 

CudaDeviceFunction void CalcRelaxation() {
    if ((NodeType & NODE_YBOUNDARY) == NODE_Wall) {
        BounceBack();          // reverses the populations that streamed in
        rho = rho(0,0);        // keep the other fields unchanged
        ux  = ux(0,0);
        uy  = uy(0,0);
        theta = theta(0,0);
        s   = s(0,0);
        return;
    }
    
    real_t rho_o   = rho(0,0);
    real_t ux_o    = ux(0,0);
    real_t uy_o    = uy(0,0);
    real_t theta_o = theta(0,0);
    real_t s_o     = s(0,0);

    real_t mu    = Viscosity;
    real_t p_loc = CalcPressure(rho_o, theta_o);
    real_t tau   = mu / p_loc;
    real_t tau_t = tau + 0.5;
    real_t omega = 1.0 / tau_t;
    real_t lam   = mu * CP() / Pr;

    //Entropy, step 7 of algorithm
    real_t inv_rhotheta = 1.0 / (rho_o * theta_o);

    real_t s_n = s_o
               - AdvectEntropy(ux_o, uy_o)
               + inv_rhotheta * Fourier_term(lam)
               + inv_rhotheta * PhiTerm(tau, tau_t);
    //Calculate Equillibrium Locally
    Equilibrium(rho_o, theta_o, ux_o, uy_o, feq);

    //Populations, step 5 of algorithm, Eq. (4.3.23)
    //Calculate Gi
    real_t g1[9];

    OffEquilibrium(ux_o , uy_o, theta_o, a1xx(0,0), a1xy(0,0), a1yy(0,0), g1);

    //Calculate Psi
    real_t psi[9];
    CalcPsi(psi, Exx(0,0), Eyy(0,0));

    for (int k = 0; k < 9; k++)
        g[k] = feq[k] + (1.0 - omega)* g1[k] + 0.5 * psi[k];

    s = s_n;
}

//Stage 2 step 6 of algorithm

CudaDeviceFunction void CalcMoments() {
    //Boundary condition fixing population coming from outside domain.
    //feq is still at level n here - PsiEQ has not run yet this iteration.
    


    switch (NodeType & NODE_XBOUNDARY) {
    case NODE_WOutflow: WOutflow(g, feq); break;
    case NODE_EOutflow: EOutflow(g, feq); break;
    }

    switch (NodeType & NODE_YBOUNDARY) {
    case NODE_SOutflow: SOutflow(g, feq); break;
    case NODE_NOutflow: NOutflow(g, feq); break;
    }


    real_t rho_n = g[0] + g[1] + g[2] + g[3] + g[4] + g[5] + g[6] + g[7] + g[8];

    rho = rho_n;
    ux  = (g[8] - g[7] - g[6] + g[5] - g[3] + g[1])/rho_n;
    uy = (g[2] - g[4] + g[5] + g[6] - g[7] - g[8])/rho_n;
}

//Stage 3 Step 8 of algorithm

CudaDeviceFunction void CalcTheta() {
    theta = ThetaFromS(rho(0,0), s(0,0));
}

//Stage 4 Step 9 and 10 of algorithm

CudaDeviceFunction void CalcEqnPsi() {
    real_t rho_n   = rho(0,0);
    real_t ux_n = ux(0,0);
    real_t uy_n = uy(0,0);
    real_t theta_n = theta(0,0);

    real_t E2ab= CalcE2ab(rho_n, theta_n, DUXDX(), DUYDY());

    //calculate E1xx First
    real_t GmX = Gamma(rhoAtX(-1), thetaAtX(-1),uxAt(-1));
    real_t G0X  = Gamma(rho_n,     theta_n,  ux_n);
    real_t GpX  = Gamma(rhoAtX(1), thetaAtX(1),uxAt(1));
    real_t sgnX = (real_t)((ux_n > 0.0) - (ux_n < 0.0));

    real_t E1xx = CalcE1ab(G0X, GpX, GmX, sgnX);    //step 9

    //Repeat for E1yy
    real_t GmY = Gamma(rhoAtY(-1), thetaAtY(-1),uyAt(-1));
    real_t G0Y  = Gamma(rho_n,     theta_n,  uy_n);
    real_t GpY  = Gamma(rhoAtY(1), thetaAtY(1),uyAt(1));
    real_t sgnY = (real_t)((uy_n > 0.0) - (uy_n < 0.0));

    real_t E1yy= CalcE1ab(G0Y, GpY, GmY, sgnY);    //step 9

    Exx = E1xx + E2ab;
    Eyy = E1yy + E2ab;

    //Calculate Equillibrium
    Equilibrium(rho_n, theta_n, ux_n, uy_n, feq);
}

//Stage 5 step 11 of algorithm

CudaDeviceFunction void CalcA1() {


    real_t rho_n = rho(0,0);
    real_t ux_n = ux(0,0);
    real_t uy_n = uy(0,0);
    real_t theta_n = theta(0,0);

    //Finite Differences for a1FD Component 
    real_t duxdx  = DUXDX();
    real_t duxdy = DUXDY();
    real_t duydy = DUYDY();
    real_t duydx = DUYDX();


    real_t p_loc = CalcPressure(rho_n, theta_n);
    real_t tau_t = Viscosity / p_loc + 0.5;

    //Calculate Psi
    real_t psi[9];
    CalcPsi(psi, Exx(0,0),  Eyy(0,0));


    switch (NodeType & NODE_XBOUNDARY) {
    case NODE_WOutflow: WOutflow(g, feq); break;
    case NODE_EOutflow: EOutflow(g, feq); break;
    }

    switch (NodeType & NODE_YBOUNDARY) {
    case NODE_SOutflow: SOutflow(g, feq); break;
    case NODE_NOutflow: NOutflow(g, feq); break;
    }




    /* step 11, Eq. (4.3.26): project the streamed populations */
    real_t a1pr_xx = A1PR_xx(g, feq, psi);
    real_t a1pr_xy = A1PR_xy(g, feq, psi);
    real_t a1pr_yy = A1PR_yy(g, feq, psi);

    real_t a1fd_xx = A1FD_aa(tau_t, p_loc, duxdx, duydy);
    real_t a1fd_yy = A1FD_aa(tau_t, p_loc, duydy, duxdx);
    real_t a1fd_xy = A1FD_xy(tau_t, p_loc, duxdy, duydx);

    /* Eq. (4.3.25) */
    a1xx = sigma*a1pr_xx + (1.0 - sigma)*a1fd_xx;
    a1xy = sigma*a1pr_xy + (1.0 - sigma)*a1fd_xy;
    a1yy = sigma*a1pr_yy + (1.0 - sigma)*a1fd_yy;

}

//Outputs

CudaDeviceFunction real_t getRho()      { return rho(0,0); }
CudaDeviceFunction real_t getTheta()    { return theta(0,0); }
CudaDeviceFunction real_t getEntropy()  { return s(0,0); }
CudaDeviceFunction real_t getPressure() { return CalcPressure(rho(0,0), theta(0,0)); }

CudaDeviceFunction vector_t getU() {
    vector_t ret;
    ret.x = ux(0,0);
    ret.y = uy(0,0);
    ret.z = 0.0;
    return ret;
}

CudaDeviceFunction float2 Color() {
    float2 ret;
    ret.x = 0;
    ret.y = 1;
    return ret;
}