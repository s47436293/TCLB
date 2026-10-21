CudaDeviceFunction void Init(){
    real_t u[2] = {Velocity, 0.0};
    real_t d = Density;
    SetEquilibrium(d,u);
}

CudaDeviceFunction void SetEquilibrium(real_t d, real_t u[2])
{
f[0] = ( 2. + ( -u[1]*u[1] - u[0]*u[0] )*3. )*d*2./9.;
f[1] = ( 2. + ( -u[1]*u[1] + ( 1 + u[0] )*u[0]*2. )*3. )*d/18.;
f[2] = ( 2. + ( -u[0]*u[0] + ( 1 + u[1] )*u[1]*2. )*3. )*d/18.;
f[3] = ( 2. + ( -u[1]*u[1] + ( -1 + u[0] )*u[0]*2. )*3. )*d/18.;
f[4] = ( 2. + ( -u[0]*u[0] + ( -1 + u[1] )*u[1]*2. )*3. )*d/18.;
f[5] = ( 1. + ( ( 1 + u[1] )*u[1] + ( 1 + u[0] + u[1]*3. )*u[0] )*3. )*d/36.;
f[6] = ( 1. + ( ( 1 + u[1] )*u[1] + ( -1 + u[0] - u[1]*3. )*u[0] )*3. )*d/36.;
f[7] = ( 1. + ( ( -1 + u[1] )*u[1] + ( -1 + u[0] + u[1]*3. )*u[0] )*3. )*d/36.;
f[8] = ( 1. + ( ( -1 + u[1] )*u[1] + ( 1 + u[0] - u[1]*3. )*u[0] )*3. )*d/36.;
}

CudaDeviceFunction void Run() {
// This defines the dynamics that we run at each node in the domain.
    switch (NodeType & NODE_BOUNDARY) {
    case NODE_Wall:
        BounceBack();
        break;
    }
    if ((NodeType & NODE_COLLISION) == NODE_BGK) 
    {
        CollisionBGK();
    }
}

CudaDeviceFunction float2 Color() {
  float2 ret;
  ret.x = 0;
  ret.y = 1;
  return ret;
}

CudaDeviceFunction void BounceBack() {
// Method to reverse distribution functions along the bounding nodes.
        real_t uf;
    uf = f[3];
    f[3] = f[1];
    f[1] = uf;
    uf = f[4];
    f[4] = f[2];
    f[2] = uf;
    uf = f[7];
    f[7] = f[5];
    f[5] = uf;
    uf = f[8];
    f[8] = f[6];
    f[6] = uf;
}

CudaDeviceFunction void CollisionBGK() {
// Here we perform a single relaxation time collision operation.
// We save memory here by using a single dummy variable
    real_t u[2], d, f_temp[9];
    d = getRho();
    // pu* = pu + rG
    u[0] = (( f[8]-f[7]-f[6]+f[5]-f[3]+f[1] )/d + GravitationX/omega );
    u[1] = ((-f[8]-f[7]+f[6]+f[5]-f[4]+f[2] )/d + GravitationY/omega );
    f_temp[0] = f[0];
    f_temp[1] = f[1];
    f_temp[2] = f[2];
    f_temp[3] = f[3];
    f_temp[4] = f[4];
    f_temp[5] = f[5];
    f_temp[6] = f[6];
    f_temp[7] = f[7];
    f_temp[8] = f[8];
    SetEquilibrium(d, u); //stores equilibrium distribution in f[0]-f[8]
    f[0] = f_temp[0] - omega*(f_temp[0]-f[0]);  
    f[1] = f_temp[1] - omega*(f_temp[1]-f[1]);
    f[2] = f_temp[2] - omega*(f_temp[2]-f[2]);
    f[3] = f_temp[3] - omega*(f_temp[3]-f[3]);  
    f[4] = f_temp[4] - omega*(f_temp[4]-f[4]);
    f[5] = f_temp[5] - omega*(f_temp[5]-f[5]);
    f[6] = f_temp[6] - omega*(f_temp[6]-f[6]);  
    f[7] = f_temp[7] - omega*(f_temp[7]-f[7]);
    f[8] = f_temp[8] - omega*(f_temp[8]-f[8]);
}

CudaDeviceFunction real_t getRho() {
// This function defines the macroscopic density at the current node.
    return f[8]+f[7]+f[6]+f[5]+f[4]+f[3]+f[2]+f[1]+f[0];
}

CudaDeviceFunction vector_t getU() {
// This function defines the macroscopic velocity at the current node.
    real_t d = f[8]+f[7]+f[6]+f[5]+f[4]+f[3]+f[2]+f[1]+f[0];
    vector_t u;
    // pv = pu + G/2
    u.x = (( f[8]-f[7]-f[6]+f[5]-f[3]+f[1] )/d + GravitationX*0.5 );
    u.y = ((-f[8]-f[7]+f[6]+f[5]-f[4]+f[2] )/d + GravitationY*0.5 );
    u.z = 0;
    return u;
}