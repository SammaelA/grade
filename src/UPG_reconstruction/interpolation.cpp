#include "interpolation.h"

const std::vector<float> 
interpolation::create_A(const std::vector<LiteMath::float3>& X)
{
    std::vector<float>A;

    for (int point_ind = 0; point_ind < 64; point_ind++)
    {
        for (int i = 0; i < interpolation_power; i++)
        {
            for (int j = 0; j < interpolation_power; j++)
            {
                for (int k = 0; k < interpolation_power; k++)
                {
                    A.push_back(
                        std::pow(X[point_ind].x, i) * 
                        std::pow(X[point_ind].y, j) * 
                        std::pow(X[point_ind].z, k)
                    );
                }
            } 
        }
    }

    return A;
}

void 
interpolation::QR(const std::vector<float> &M, const size_t &size, std::vector<float> &Q, std::vector<float> &R)
{
    for (int ind = 0; ind < size; ind++)
    {
        std::vector<float> cur_vector(size, 0);

        for (int i = 0; i < size; i++)
        {
            cur_vector[i] = M[i * size + ind];
        }

        //  Substitute projects
        for (int j = 0; j < ind; j++)
        {
            //  calc proj
            float proj_coef = 0;

            for (int i = 0; i < size; i++)
            {
                proj_coef += cur_vector[i] * Q[i * size + j];
            }

            for (int i = 0; i < size; i++)
            {
                cur_vector[i] -= proj_coef * Q[i * size + j];
            }
        }

        //  norm new basis vector
        float size_vec = 0;

        for (int i = 0; i < size; i++)
        {
            size_vec += cur_vector[i] * cur_vector[i];
        }
        
        size_vec = std::sqrt(size_vec);

        for (int i = 0; i < size; i++)
        {
            cur_vector[i] /= size_vec;
        }

        for (int i = 0; i < size; i++)
        {
            Q[i * size + ind] = cur_vector[i];
        }
    }

    for (int i = 0; i < size; i++)
    {
        for (int j = i; j < size; j++)
        {
            float dot_product = 0;

            for (int k = 0; k < size; k++)
            {
                dot_product += M[k * size + j] * Q[k * size + i];
            }

            R[i * size + j] = dot_product;
        }
    }
}

float 
interpolation::matrix_norm(const std::vector<float> &A1, const std::vector<float> &A2)
{
    float eps = 0;
    int size = A1.size();

    for (int i = 0; i < size; i++)
    {
        eps += std::pow(A1[i] - A2[i], 2);
    }

    eps = std::sqrt(eps);

    return eps;
}

std::vector<float> 
interpolation::calc_qr_coefs(const std::vector<float> &Q, const std::vector<float> &R, const std::vector<float> &b)
{
    std::vector<float> x(interpolation_power * interpolation_power * interpolation_power, 0);
    std::vector<float> v(interpolation_power * interpolation_power * interpolation_power, 0);

    size_t size = b.size();

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            v[i] += Q[j * size + i] * b[j];
        }
    }

    //  back-substitution
    for (int i = size - 1; i >= 0; i--)
    {
        float s = 0;
        for (int j = i + 1; j < size; j++)
        {
            s += R[i * size + j] * x[j];
        }

        x[i] = (v[i] - s) / R[i * size + i];
    }

    return x;
}

std::vector<float> 
interpolation::mul_qr(const std::vector<float> &Q, const std::vector<float> &R, size_t size)
{
    std::vector<float> G(size * size, 0);

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            for (int k = 0; k < size; k++)
            {
                G[i * size + j] += Q[i * size + k] * R[k * size + j];
            }
        }
    }

    return G;
}

float 
interpolation::perform_interpolation(const std::vector<float> &coefs, const LiteMath::float3 &pos)
{
    float res = 0;
    int ind = 0;

    for (int i = 0; i < interpolation_power; i++)
    {
        for (int j = 0; j < interpolation_power; j++)
        {
            for (int k = 0; k < interpolation_power; k++)
            {
                res += coefs[ind] * std::pow(pos.x, i) * std::pow(pos.y, j) * std::pow(pos.z, k);
                ind++;
            }
        } 
    }

    return res;
}

void 
interpolation::householder_qr(const std::vector<float>& M, const size_t &size, std::vector<float> &Q, std::vector<float> &R)
{
    std::vector<float> I(size * size, 0);
    std::vector<float> P(size * size, 0);

    //  Init identity matrix
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            if (i == j)
            {
                I[i * size + j] = 1;
            }
        }
    }

    P = I;
    Q = I;
    R = M;

    float eps = 0.00000001;

    for (int i = 0; i < size; i++)
    {
        std::vector<float> v(size, 0), u(size, 0);

        for (int j = i; j < size; j++)
        {
            u[j] = R[j * size + i];
        }
        
        float u_size = 0;

        for (int j = 0; j < size; j++)
        {
            u_size += u[j] * u[j];
        }

        u_size = std::sqrt(u_size);

        float alpha = (u[i] < 0) ? u_size : -u_size;

        for (int j = 0; j < size; j++)
        {
            v[j] = (j == i) ? u[j] + alpha : u[j];
        }

        float v_size = 0;

        for (int j = 0; j < size; j++)
        {
            v_size += v[j] * v[j];
        }

        v_size = std::sqrt(v_size);

        if (v_size < eps) continue;

        for (int j = 0; j < size; j++)
        {
            v[j] /= v_size;
        }

        std::vector<float> v_mat(size * size, 0);

        for (int j = 0; j < size; j++)
        {
            for (int k = 0; k < size; k++)
            {
                v_mat[j * size + k] = v[j] * v[k];
            }
        }

        //  calc new P
        for (int j = 0; j < size; j++)
        {
            for (int k = 0; k < size; k++)
            {
                P[j * size + k] = I[j * size + k] - 2 * v_mat[j * size + k];
            }
        }

        std::vector<float> new_R(size * size, 0), new_Q(size * size, 0);

        for (int j = 0; j < size; j++)
        {
            for (int k = 0; k < size; k++)
            {
                for (int z = 0; z < size; z++)
                {
                    new_R[j * size + k] += P[j * size + z] * R[z * size + k];
                    new_Q[j * size + k] += Q[j * size + z] * P[z * size + k];
                }
            }
        }

        R = new_R;
        Q = new_Q;

        new_Q.clear();
        new_R.clear();
        v.clear();
        u.clear();
    }
}

void partialvec_copy (double * x, double * y, int length, int index) {
    int i, length5;

    length5 = length % 5;

    for(i = 0; i < length5; i++) {
        y[i] = x[i + index];
    }
    for(; i < length; i += 5) {
        y[i] = x[i + index];
        y[i + 1] = x[i + index + 1];
        y[i + 2] = x[i + index + 2];
        y[i + 3] = x[i + index + 3];
        y[i + 4] = x[i + index + 4];
    }
}

void scalar_div (double * x, double r, int length, double * y) {
    int i, length5;

    length5 = length % 5;

    for(i = 0; i < length5; i++) {
        y[i] = x[i]/r;
    }
    for(; i < length; i += 5) {
        y[i] = x[i]/r;
        y[i + 1] = x[i + 1]/r;
        y[i + 2] = x[i + 2]/r;
        y[i + 3] = x[i + 3]/r;
        y[i + 4] = x[i + 4]/r;
    }
}

double partialdot_product (double * x, double * y, int length, int index) {
    int i, length5;
    double sum = 0;

    length5 = length % 5;

    for(i = index; i < length5; i++) {
        sum += x[i] * y[i];
    }
    for(; i < length; i += 5) {
        sum += x[i] * y[i] + x[i + 1] * y[i + 1] + x[i + 2] * y[i + 2]
                           + x[i + 3] * y[i + 3] + x[i + 4] * y[i + 4];
    }

    return sum;
}

double subdot_product (double * x, double * y, int length, int index) {
    int i, length5;
    double sum = 0;

    length5 = length % 5;

    for(i = 0; i < length5; i++) {
        sum += x[i + index] * y[i];
    }
    for(; i < length; i += 5) {
        sum += x[i + index] * y[i] + x[i + index + 1] * y[i + 1] 
                                   + x[i + index + 2] * y[i + 2]
                                   + x[i + index + 3] * y[i + 3]
                                   + x[i + index + 4] * y[i + 4];
    }

    return sum;
}

void partialscalar_sub (double * x, double r, int length, 
                                              int index, double * y) 
{
    int i, length5;

    length5 = length % 5;

    for(i = 0; i < length5; i++) {
        y[i + index] -= r * x[i];
    }
    for(; i < length; i += 5) {
        y[i + index] -= r * x[i];
        y[i + index + 1] -= r * x[i + 1];
        y[i + index + 2] -= r * x[i + 2];
        y[i + index + 3] -= r * x[i + 3];
        y[i + index + 4] -= r * x[i + 4];
    }
}

void householder (double ** a, double ** v, int m, int n) {
    int i, j;
    double vnorm, vTa, vpartdot;

    for(i = 0; i < n; i++) {
        /* set v[i] equal to subvector a[i][i : m] */
        partialvec_copy(a[i], v[i], m - i, i);

        /* vpartdot = ||v[i]||^2 - v[i][0] * v[i][0]; since vpartdot 
           is unaffected by the change in v[i][0], storing this value 
           prevents the need to recalculate the entire norm of v[i] 
           after updating v[i][0] in the following step              */
        vpartdot = partialdot_product(v[i], v[i], m - i, 1);

        /* set v[i][0] = v[i][0] + sign(v[i][0]) * ||v[i]|| */
        if(v[i][0] < 0) {
            v[i][0] -= sqrt(v[i][0] * v[i][0] + vpartdot);
        }
        else {
            v[i][0] += sqrt(v[i][0] * v[i][0] + vpartdot);
        }

        /* normalize v[i] */
        vnorm = sqrt(v[i][0] * v[i][0] + vpartdot);
        scalar_div(v[i], vnorm, m - i, v[i]);
    
        for(j = i; j < n; j++) {
            /* set a[j][i:m] = a[j][i:m] - 2 * (v[i]^T a[j][i:m]) * v[i] */
            vTa = subdot_product(a[j], v[i], m - i, i);
            vTa *= 2;
            partialscalar_sub(v[i], vTa, m - i, i, a[j]);
        }
    }
}