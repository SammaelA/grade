#include "functions.h"
#include "LiteMath.h"
namespace diff_render
{
/*
static inline float VS_X(float V[3], const CamInfo& data) // same as VertexShader().x
{
  const float W    =  V[0] * data.mWVP[3] + V[1] * data.mWVP[7] + V[2] * data.mWVP[11] + data.mWVP[15]; 
  const float xNDC = (V[0] * data.mWVP[0] + V[1] * data.mWVP[4] + V[2] * data.mWVP[ 8] + data.mWVP[12])/W;
  return (xNDC*0.5f + 0.5f)*data.width;
}

static inline float VS_Y(float V[3], const CamInfo& data) // // same as VertexShader().y
{
  const float W    =   V[0] * data.mWVP[3] + V[1] * data.mWVP[7] + V[2] * data.mWVP[11] + data.mWVP[15]; 
  const float xNDC = -(V[0] * data.mWVP[1] + V[1] * data.mWVP[5] + V[2] * data.mWVP[ 9] + data.mWVP[13])/W;
  return (xNDC*0.5f + 0.5f)*data.height;
}
*/

void VS_X_grad(float V[3], const CamInfo &data, float _d_V[3]) 
{
    float _t0;
    float _t1;
    float _t2;
    float _t3;
    float _t4;
    float _t5;
    float _d_W = 0;
    float _t6;
    float _t7;
    float _t8;
    float _t9;
    float _t10;
    float _t11;
    float _t12;
    float _t13;
    float _d_xNDC = 0;
    float _t14;
    float _t15;
    _t1 = V[0];
    _t0 = data.mWVP[3];
    _t3 = V[1];
    _t2 = data.mWVP[7];
    _t5 = V[2];
    _t4 = data.mWVP[11];
    const float W = _t1 * _t0 + _t3 * _t2 + _t5 * _t4 + data.mWVP[15];
    _t8 = V[0];
    _t7 = data.mWVP[0];
    _t10 = V[1];
    _t9 = data.mWVP[4];
    _t12 = V[2];
    _t11 = data.mWVP[8];
    _t13 = (_t8 * _t7 + _t10 * _t9 + _t12 * _t11 + data.mWVP[12]);
    _t6 = W;
    const float xNDC = _t13 / _t6;
    _t15 = (xNDC * 0.5F + 0.5F);
    _t14 = data.width;
    float VS_X_return = _t15 * _t14;
    {
        float _r14 = 1 * _t14;
        float _r15 = _r14 * 0.5F;
        _d_xNDC += _r15;
        float _r16 = _t15 * 1;
    }
    {
        float _r6 = _d_xNDC / _t6;
        float _r7 = _r6 * _t7;
        _d_V[0] += _r7;
        float _r8 = _t8 * _r6;
        float _r9 = _r6 * _t9;
        _d_V[1] += _r9;
        float _r10 = _t10 * _r6;
        float _r11 = _r6 * _t11;
        _d_V[2] += _r11;
        float _r12 = _t12 * _r6;
        float _r13 = _d_xNDC * -_t13 / (_t6 * _t6);
        _d_W += _r13;
    }
    {
        float _r0 = _d_W * _t0;
        _d_V[0] += _r0;
        float _r1 = _t1 * _d_W;
        float _r2 = _d_W * _t2;
        _d_V[1] += _r2;
        float _r3 = _t3 * _d_W;
        float _r4 = _d_W * _t4;
        _d_V[2] += _r4;
        float _r5 = _t5 * _d_W;
    }
}

void VS_Y_grad(float V[3], const CamInfo &data, float _d_V[3]) 
{
    float _t0;
    float _t1;
    float _t2;
    float _t3;
    float _t4;
    float _t5;
    float _d_W = 0;
    float _t6;
    float _t7;
    float _t8;
    float _t9;
    float _t10;
    float _t11;
    float _t12;
    float _t13;
    float _d_xNDC = 0;
    float _t14;
    float _t15;
    _t1 = V[0];
    _t0 = data.mWVP[3];
    _t3 = V[1];
    _t2 = data.mWVP[7];
    _t5 = V[2];
    _t4 = data.mWVP[11];
    const float W = _t1 * _t0 + _t3 * _t2 + _t5 * _t4 + data.mWVP[15];
    _t8 = V[0];
    _t7 = data.mWVP[1];
    _t10 = V[1];
    _t9 = data.mWVP[5];
    _t12 = V[2];
    _t11 = data.mWVP[9];
    _t13 = -(_t8 * _t7 + _t10 * _t9 + _t12 * _t11 + data.mWVP[13]);
    _t6 = W;
    const float xNDC = _t13 / _t6;
    _t15 = (xNDC * 0.5F + 0.5F);
    _t14 = data.height;
    float VS_Y_return = _t15 * _t14;
    {
        float _r14 = 1 * _t14;
        float _r15 = _r14 * 0.5F;
        _d_xNDC += _r15;
        float _r16 = _t15 * 1;
    }
    {
        float _r6 = _d_xNDC / _t6;
        float _r7 = -_r6 * _t7;
        _d_V[0] += _r7;
        float _r8 = _t8 * -_r6;
        float _r9 = -_r6 * _t9;
        _d_V[1] += _r9;
        float _r10 = _t10 * -_r6;
        float _r11 = -_r6 * _t11;
        _d_V[2] += _r11;
        float _r12 = _t12 * -_r6;
        float _r13 = _d_xNDC * -_t13 / (_t6 * _t6);
        _d_W += _r13;
    }
    {
        float _r0 = _d_W * _t0;
        _d_V[0] += _r0;
        float _r1 = _t1 * _d_W;
        float _r2 = _d_W * _t2;
        _d_V[1] += _r2;
        float _r3 = _t3 * _d_W;
        float _r4 = _d_W * _t4;
        _d_V[2] += _r4;
        float _r5 = _t5 * _d_W;
    }
}

/*
inline float sign(float x) { return (x >= 0) ? 1 : -1; }

static inline float BarU( const float ray_pos[3], const float ray_dir[3], const float A[3], const float B[3], const float C[3])
{
  const float edge1X = B[0] - A[0];
  const float edge1Y = B[1] - A[1];
  const float edge1Z = B[2] - A[2];

  const float edge2X = C[0] - A[0];
  const float edge2Y = C[1] - A[1];
  const float edge2Z = C[2] - A[2];
  
  const float pvecZ = ray_dir[0]*edge2Y - ray_dir[1]*edge2X;
  const float pvecX = ray_dir[1]*edge2Z - ray_dir[2]*edge2Y;
  const float pvecY = ray_dir[2]*edge2X - ray_dir[0]*edge2Z;

  const float tvecX  = ray_pos[0] - A[0];
  const float tvecY  = ray_pos[1] - A[1];
  const float tvecZ  = ray_pos[2] - A[2];

  const float qvecZ  = tvecX*edge1Y - tvecY*edge1X;
  const float qvecX  = tvecY*edge1Z - tvecZ*edge1Y;
  const float qvecY  = tvecZ*edge1X - tvecX*edge1Z;

  const float e1dp   = edge1X*pvecX + edge1Y*pvecY + edge1Z*pvecZ;
  const float signv  = sign(e1dp); // put 1.0 to enable triangle clippin
  
  return signv*(qvecX*ray_dir[0] + qvecY*ray_dir[1] + qvecZ*ray_dir[2])/(e1dp*signv);
}

static inline float BarV(const float ray_pos[3], const float ray_dir[3], const float A[3], const float B[3], const float C[3])
{
  const float edge1X = B[0] - A[0];
  const float edge1Y = B[1] - A[1];
  const float edge1Z = B[2] - A[2];

  const float edge2X = C[0] - A[0];
  const float edge2Y = C[1] - A[1];
  const float edge2Z = C[2] - A[2];

  const float pvecZ = ray_dir[0]*edge2Y - ray_dir[1]*edge2X;
  const float pvecX = ray_dir[1]*edge2Z - ray_dir[2]*edge2Y;
  const float pvecY = ray_dir[2]*edge2X - ray_dir[0]*edge2Z;

  const float tvecX  = ray_pos[0] - A[0];
  const float tvecY  = ray_pos[1] - A[1];
  const float tvecZ  = ray_pos[2] - A[2];

  const float qvecZ  = tvecX*edge1Y - tvecY*edge1X;
  const float qvecX  = tvecY*edge1Z - tvecZ*edge1Y;
  const float qvecY  = tvecZ*edge1X - tvecX*edge1Z;

  const float e1dp   = edge1X*pvecX + edge1Y*pvecY + edge1Z*pvecZ;
  const float signv  = sign(e1dp); // put 1.0 to enable triangle clippin

  return signv*(tvecX*pvecX + tvecY*pvecY + tvecZ*pvecZ)/(signv*e1dp); 
}
*/


void BarV_grad(const float ray_pos[3], const float ray_dir[3], const float A[3], const float B[3], const float C[3], 
               float* _d_A, float* _d_B, float* _d_C) 
{
    float _d_edge1X = 0;
    float _d_edge1Y = 0;
    float _d_edge1Z = 0;
    float _d_edge2X = 0;
    float _d_edge2Y = 0;
    float _d_edge2Z = 0;
    float _t0;
    float _t1;
    float _t2;
    float _t3;
    float _d_pvecZ = 0;
    float _t4;
    float _t5;
    float _t6;
    float _t7;
    float _d_pvecX = 0;
    float _t8;
    float _t9;
    float _t10;
    float _t11;
    float _d_pvecY = 0;
    float _d_tvecX = 0;
    float _d_tvecY = 0;
    float _d_tvecZ = 0;
    float _t12;
    float _t13;
    float _t14;
    float _t15;
    float _d_qvecZ = 0;
    float _t16;
    float _t17;
    float _t18;
    float _t19;
    float _d_qvecX = 0;
    float _t20;
    float _t21;
    float _t22;
    float _t23;
    float _d_qvecY = 0;
    float _t24;
    float _t25;
    float _t26;
    float _t27;
    float _t28;
    float _t29;
    float _d_e1dp = 0;
    float _t30;
    float _d_signv = 0;
    float _t31;
    float _t32;
    float _t33;
    float _t34;
    float _t35;
    float _t36;
    float _t37;
    float _t38;
    float _t39;
    float _t40;
    float _t41;
    float _t42;
    const float edge1X = B[0] - A[0];
    const float edge1Y = B[1] - A[1];
    const float edge1Z = B[2] - A[2];
    const float edge2X = C[0] - A[0];
    const float edge2Y = C[1] - A[1];
    const float edge2Z = C[2] - A[2];
    _t1 = ray_dir[0];
    _t0 = edge2Y;
    _t3 = ray_dir[1];
    _t2 = edge2X;
    const float pvecZ = _t1 * _t0 - _t3 * _t2;
    _t5 = ray_dir[1];
    _t4 = edge2Z;
    _t7 = ray_dir[2];
    _t6 = edge2Y;
    const float pvecX = _t5 * _t4 - _t7 * _t6;
    _t9 = ray_dir[2];
    _t8 = edge2X;
    _t11 = ray_dir[0];
    _t10 = edge2Z;
    const float pvecY = _t9 * _t8 - _t11 * _t10;
    const float tvecX = ray_pos[0] - A[0];
    const float tvecY = ray_pos[1] - A[1];
    const float tvecZ = ray_pos[2] - A[2];
    _t13 = tvecX;
    _t12 = edge1Y;
    _t15 = tvecY;
    _t14 = edge1X;
    const float qvecZ = _t13 * _t12 - _t15 * _t14;
    _t17 = tvecY;
    _t16 = edge1Z;
    _t19 = tvecZ;
    _t18 = edge1Y;
    const float qvecX = _t17 * _t16 - _t19 * _t18;
    _t21 = tvecZ;
    _t20 = edge1X;
    _t23 = tvecX;
    _t22 = edge1Z;
    const float qvecY = _t21 * _t20 - _t23 * _t22;
    _t25 = edge1X;
    _t24 = pvecX;
    _t27 = edge1Y;
    _t26 = pvecY;
    _t29 = edge1Z;
    _t28 = pvecZ;
    const float e1dp = _t25 * _t24 + _t27 * _t26 + _t29 * _t28;
    _t30 = e1dp;
    const float signv = LiteMath::sign(_t30);
    _t33 = signv;
    _t35 = tvecX;
    _t34 = pvecX;
    _t37 = tvecY;
    _t36 = pvecY;
    _t39 = tvecZ;
    _t38 = pvecZ;
    _t32 = (_t35 * _t34 + _t37 * _t36 + _t39 * _t38);
    _t40 = _t33 * _t32;
    _t42 = signv;
    _t41 = e1dp;
    _t31 = (_t42 * _t41);
    float BarV_return = _t40 / _t31;
    {
        float _r31 = 1 / _t31;
        float _r32 = _r31 * _t32;
        _d_signv += _r32;
        float _r33 = _t33 * _r31;
        float _r34 = _r33 * _t34;
        _d_tvecX += _r34;
        float _r35 = _t35 * _r33;
        _d_pvecX += _r35;
        float _r36 = _r33 * _t36;
        _d_tvecY += _r36;
        float _r37 = _t37 * _r33;
        _d_pvecY += _r37;
        float _r38 = _r33 * _t38;
        _d_tvecZ += _r38;
        float _r39 = _t39 * _r33;
        _d_pvecZ += _r39;
        float _r40 = 1 * -_t40 / (_t31 * _t31);
        float _r41 = _r40 * _t41;
        _d_signv += _r41;
        float _r42 = _t42 * _r40;
        _d_e1dp += _r42;
    }
    {
        float _r24 = _d_e1dp * _t24;
        _d_edge1X += _r24;
        float _r25 = _t25 * _d_e1dp;
        _d_pvecX += _r25;
        float _r26 = _d_e1dp * _t26;
        _d_edge1Y += _r26;
        float _r27 = _t27 * _d_e1dp;
        _d_pvecY += _r27;
        float _r28 = _d_e1dp * _t28;
        _d_edge1Z += _r28;
        float _r29 = _t29 * _d_e1dp;
        _d_pvecZ += _r29;
    }
    {
        float _r20 = _d_qvecY * _t20;
        _d_tvecZ += _r20;
        float _r21 = _t21 * _d_qvecY;
        _d_edge1X += _r21;
        float _r22 = -_d_qvecY * _t22;
        _d_tvecX += _r22;
        float _r23 = _t23 * -_d_qvecY;
        _d_edge1Z += _r23;
    }
    {
        float _r16 = _d_qvecX * _t16;
        _d_tvecY += _r16;
        float _r17 = _t17 * _d_qvecX;
        _d_edge1Z += _r17;
        float _r18 = -_d_qvecX * _t18;
        _d_tvecZ += _r18;
        float _r19 = _t19 * -_d_qvecX;
        _d_edge1Y += _r19;
    }
    {
        float _r12 = _d_qvecZ * _t12;
        _d_tvecX += _r12;
        float _r13 = _t13 * _d_qvecZ;
        _d_edge1Y += _r13;
        float _r14 = -_d_qvecZ * _t14;
        _d_tvecY += _r14;
        float _r15 = _t15 * -_d_qvecZ;
        _d_edge1X += _r15;
    }
    _d_A[2] += -_d_tvecZ;
    _d_A[1] += -_d_tvecY;
    _d_A[0] += -_d_tvecX;
    {
        float _r8 = _d_pvecY * _t8;
        float _r9 = _t9 * _d_pvecY;
        _d_edge2X += _r9;
        float _r10 = -_d_pvecY * _t10;
        float _r11 = _t11 * -_d_pvecY;
        _d_edge2Z += _r11;
    }
    {
        float _r4 = _d_pvecX * _t4;
        float _r5 = _t5 * _d_pvecX;
        _d_edge2Z += _r5;
        float _r6 = -_d_pvecX * _t6;
        float _r7 = _t7 * -_d_pvecX;
        _d_edge2Y += _r7;
    }
    {
        float _r0 = _d_pvecZ * _t0;
        float _r1 = _t1 * _d_pvecZ;
        _d_edge2Y += _r1;
        float _r2 = -_d_pvecZ * _t2;
        float _r3 = _t3 * -_d_pvecZ;
        _d_edge2X += _r3;
    }
    {
        _d_C[2] += _d_edge2Z;
        _d_A[2] += -_d_edge2Z;
    }
    {
        _d_C[1] += _d_edge2Y;
        _d_A[1] += -_d_edge2Y;
    }
    {
        _d_C[0] += _d_edge2X;
        _d_A[0] += -_d_edge2X;
    }
    {
        _d_B[2] += _d_edge1Z;
        _d_A[2] += -_d_edge1Z;
    }
    {
        _d_B[1] += _d_edge1Y;
        _d_A[1] += -_d_edge1Y;
    }
    {
        _d_B[0] += _d_edge1X;
        _d_A[0] += -_d_edge1X;
    }
}

void BarU_grad(const float ray_pos[3], const float ray_dir[3], const float A[3], const float B[3], const float C[3], 
               float* _d_A, float* _d_B, float* _d_C) 
{
    float _d_edge1X = 0;
    float _d_edge1Y = 0;
    float _d_edge1Z = 0;
    float _d_edge2X = 0;
    float _d_edge2Y = 0;
    float _d_edge2Z = 0;
    float _t0;
    float _t1;
    float _t2;
    float _t3;
    float _d_pvecZ = 0;
    float _t4;
    float _t5;
    float _t6;
    float _t7;
    float _d_pvecX = 0;
    float _t8;
    float _t9;
    float _t10;
    float _t11;
    float _d_pvecY = 0;
    float _d_tvecX = 0;
    float _d_tvecY = 0;
    float _d_tvecZ = 0;
    float _t12;
    float _t13;
    float _t14;
    float _t15;
    float _d_qvecZ = 0;
    float _t16;
    float _t17;
    float _t18;
    float _t19;
    float _d_qvecX = 0;
    float _t20;
    float _t21;
    float _t22;
    float _t23;
    float _d_qvecY = 0;
    float _t24;
    float _t25;
    float _t26;
    float _t27;
    float _t28;
    float _t29;
    float _d_e1dp = 0;
    float _t30;
    float _d_signv = 0;
    float _t31;
    float _t32;
    float _t33;
    float _t34;
    float _t35;
    float _t36;
    float _t37;
    float _t38;
    float _t39;
    float _t40;
    float _t41;
    float _t42;
    const float edge1X = B[0] - A[0];
    const float edge1Y = B[1] - A[1];
    const float edge1Z = B[2] - A[2];
    const float edge2X = C[0] - A[0];
    const float edge2Y = C[1] - A[1];
    const float edge2Z = C[2] - A[2];
    _t1 = ray_dir[0];
    _t0 = edge2Y;
    _t3 = ray_dir[1];
    _t2 = edge2X;
    const float pvecZ = _t1 * _t0 - _t3 * _t2;
    _t5 = ray_dir[1];
    _t4 = edge2Z;
    _t7 = ray_dir[2];
    _t6 = edge2Y;
    const float pvecX = _t5 * _t4 - _t7 * _t6;
    _t9 = ray_dir[2];
    _t8 = edge2X;
    _t11 = ray_dir[0];
    _t10 = edge2Z;
    const float pvecY = _t9 * _t8 - _t11 * _t10;
    const float tvecX = ray_pos[0] - A[0];
    const float tvecY = ray_pos[1] - A[1];
    const float tvecZ = ray_pos[2] - A[2];
    _t13 = tvecX;
    _t12 = edge1Y;
    _t15 = tvecY;
    _t14 = edge1X;
    const float qvecZ = _t13 * _t12 - _t15 * _t14;
    _t17 = tvecY;
    _t16 = edge1Z;
    _t19 = tvecZ;
    _t18 = edge1Y;
    const float qvecX = _t17 * _t16 - _t19 * _t18;
    _t21 = tvecZ;
    _t20 = edge1X;
    _t23 = tvecX;
    _t22 = edge1Z;
    const float qvecY = _t21 * _t20 - _t23 * _t22;
    _t25 = edge1X;
    _t24 = pvecX;
    _t27 = edge1Y;
    _t26 = pvecY;
    _t29 = edge1Z;
    _t28 = pvecZ;
    const float e1dp = _t25 * _t24 + _t27 * _t26 + _t29 * _t28;
    _t30 = e1dp;
    const float signv = LiteMath::sign(_t30);
    _t33 = signv;
    _t35 = qvecX;
    _t34 = ray_dir[0];
    _t37 = qvecY;
    _t36 = ray_dir[1];
    _t39 = qvecZ;
    _t38 = ray_dir[2];
    _t32 = (_t35 * _t34 + _t37 * _t36 + _t39 * _t38);
    _t40 = _t33 * _t32;
    _t42 = e1dp;
    _t41 = signv;
    _t31 = (_t42 * _t41);
    float BarW_return = _t40 / _t31;
    {
        float _r31 = 1 / _t31;
        float _r32 = _r31 * _t32;
        _d_signv += _r32;
        float _r33 = _t33 * _r31;
        float _r34 = _r33 * _t34;
        _d_qvecX += _r34;
        float _r35 = _t35 * _r33;
        float _r36 = _r33 * _t36;
        _d_qvecY += _r36;
        float _r37 = _t37 * _r33;
        float _r38 = _r33 * _t38;
        _d_qvecZ += _r38;
        float _r39 = _t39 * _r33;
        float _r40 = 1 * -_t40 / (_t31 * _t31);
        float _r41 = _r40 * _t41;
        _d_e1dp += _r41;
        float _r42 = _t42 * _r40;
        _d_signv += _r42;
    }
    {
        float _r24 = _d_e1dp * _t24;
        _d_edge1X += _r24;
        float _r25 = _t25 * _d_e1dp;
        _d_pvecX += _r25;
        float _r26 = _d_e1dp * _t26;
        _d_edge1Y += _r26;
        float _r27 = _t27 * _d_e1dp;
        _d_pvecY += _r27;
        float _r28 = _d_e1dp * _t28;
        _d_edge1Z += _r28;
        float _r29 = _t29 * _d_e1dp;
        _d_pvecZ += _r29;
    }
    {
        float _r20 = _d_qvecY * _t20;
        _d_tvecZ += _r20;
        float _r21 = _t21 * _d_qvecY;
        _d_edge1X += _r21;
        float _r22 = -_d_qvecY * _t22;
        _d_tvecX += _r22;
        float _r23 = _t23 * -_d_qvecY;
        _d_edge1Z += _r23;
    }
    {
        float _r16 = _d_qvecX * _t16;
        _d_tvecY += _r16;
        float _r17 = _t17 * _d_qvecX;
        _d_edge1Z += _r17;
        float _r18 = -_d_qvecX * _t18;
        _d_tvecZ += _r18;
        float _r19 = _t19 * -_d_qvecX;
        _d_edge1Y += _r19;
    }
    {
        float _r12 = _d_qvecZ * _t12;
        _d_tvecX += _r12;
        float _r13 = _t13 * _d_qvecZ;
        _d_edge1Y += _r13;
        float _r14 = -_d_qvecZ * _t14;
        _d_tvecY += _r14;
        float _r15 = _t15 * -_d_qvecZ;
        _d_edge1X += _r15;
    }
    _d_A[2] += -_d_tvecZ;
    _d_A[1] += -_d_tvecY;
    _d_A[0] += -_d_tvecX;
    {
        float _r8 = _d_pvecY * _t8;
        float _r9 = _t9 * _d_pvecY;
        _d_edge2X += _r9;
        float _r10 = -_d_pvecY * _t10;
        float _r11 = _t11 * -_d_pvecY;
        _d_edge2Z += _r11;
    }
    {
        float _r4 = _d_pvecX * _t4;
        float _r5 = _t5 * _d_pvecX;
        _d_edge2Z += _r5;
        float _r6 = -_d_pvecX * _t6;
        float _r7 = _t7 * -_d_pvecX;
        _d_edge2Y += _r7;
    }
    {
        float _r0 = _d_pvecZ * _t0;
        float _r1 = _t1 * _d_pvecZ;
        _d_edge2Y += _r1;
        float _r2 = -_d_pvecZ * _t2;
        float _r3 = _t3 * -_d_pvecZ;
        _d_edge2X += _r3;
    }
    {
        _d_C[2] += _d_edge2Z;
        _d_A[2] += -_d_edge2Z;
    }
    {
        _d_C[1] += _d_edge2Y;
        _d_A[1] += -_d_edge2Y;
    }
    {
        _d_C[0] += _d_edge2X;
        _d_A[0] += -_d_edge2X;
    }
    {
        _d_B[2] += _d_edge1Z;
        _d_A[2] += -_d_edge1Z;
    }
    {
        _d_B[1] += _d_edge1Y;
        _d_A[1] += -_d_edge1Y;
    }
    {
        _d_B[0] += _d_edge1X;
        _d_A[0] += -_d_edge1X;
    }
}
}

