#pragma once

namespace qmc
{
  constexpr static int   QRNG_DIMENSIONS = 11;
  constexpr static int   QRNG_RESOLUTION = 31;
  constexpr static float INT_SCALE       = (1.0f / (float)0x80000001U);

  void init(unsigned int table[QRNG_DIMENSIONS][QRNG_RESOLUTION]);

  /**
  \brief get arbitrary 'Sobol/Niederreiter' quasi random float value in range [0,1]
  \param pos     - id of value in sequence (i.e. first, second, ... )
  \param dim     - dimention/coordinate id (i.e. x,y,z,w, ... )
  \param c_Table - some table previously initialised with hr_qmc::init function

  */
  float rndFloat(unsigned int pos, int dim, unsigned int *c_Table); ///< return

  /**
  \brief get arbitrary 'Sobol/Niederreiter' quasi random float value in range [s,e].
  \param gen - pointer to current generator state
  \param s   - low  boundary of generated random number
  \param e   - high boundary of generated random number
  */
  float rndFloatUniform(unsigned int pos, int dim, unsigned int *c_Table, float s, float e);

  /**
  \brief get arbitrary 'Sobol/Niederreiter' quasi random integer value in range [s,e].
  \param gen     - pointer to current generator state
  \param s   - low  boundary of generated random number
  \param e   - high boundary of generated random number
  */
  int rndIntUniform(unsigned int pos, int dim, unsigned int *c_Table, int a, int b);
  
  /**
  \brief generate 2D Hammersley samples
  \param result - 2n floats
  \param n      - number of samples
  */
  void planeHammersley(float *result, int n);
};

namespace prng
{
  struct uint2
  {
    unsigned int x;
    unsigned int y;
  };

  typedef struct RandomGenT
  {
    uint2 state;
  } RandomGen;

  RandomGen RandomGenInit(const int a_seed);

  /**
   \brief get next pseudo random float value in range [0,1].
   \param gen - pointer to current generator state
   */
  float rndFloat(RandomGen *gen);

  /**
   \brief get next pseudo random float value in range [s,e].
   \param gen - reference to current generator state
   \param s   - low  boundary of generated random number
   \param e   - high boundary of generated random number
   */
  float rndFloatUniform(RandomGen& gen, float s, float e);

  /**
   \brief get next pseudo random integer value in range [s,e].
   \param gen - reference to current generator state
   \param s   - low  boundary of generated random number
   \param e   - high boundary of generated random number
   */
  int rndIntUniform(RandomGen& gen, int a, int b);

  static inline unsigned int NextState(RandomGen *gen)
  {
    const unsigned int x = (gen->state).x * 17 + (gen->state).y * 13123;
    (gen->state).x       = (x << 13) ^ x;
    (gen->state).y      ^= (x << 7);
    return x;
  }

  static inline int mapRndFloatToInt(float a_val, int a, int b)
  {
    const float fa = (float) (a + 0);
    const float fb = (float) (b + 1);
    const float fR = fa + a_val * (fb - fa);

    const int res = (int) (fR);

    if (res > b)
      return b;
    else
      return res;
  }
};
