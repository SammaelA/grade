#include "interpolation.h"

const std::vector<float> 
interpolation::create_A(const std::vector<LiteMath::float3>& X)
{
    std::vector<float>A;

    for (int x = 0; x < 8; x++)
    {
        for (int y = 0; y < 8; y++)
        {
            for (int i = 0; i < interpolation_power; i++)
            {
                for (int j = 0; j < interpolation_power; j++)
                {
                    for (int k = 0; k < interpolation_power; k++)
                    {
                        A.push_back(
                            std::pow(X[interpolation_power * y + x].x, i) * 
                            std::pow(X[interpolation_power * y + x].y, j) * 
                            std::pow(X[interpolation_power * y + x].z, k)
                        );
                    }
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