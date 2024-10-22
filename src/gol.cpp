#include "gol_seq.hpp"

#include <cstdint>
#include <cstring>
#include <utility>


/**
 * @brief Executes one step of the `gol` process
 *
 * @param g_m The input matrix.
 * @param s_m The output matrix.
 * @param nx The size of each row.
 * @param ny The size of each column.
 */
static void gol_kernel(bool *d_m, bool *s_m, std::size_t nx, std::size_t ny);


int gol_seq(bool *g_m, std::size_t nx, std::size_t ny, std::size_t ni)
{
  std::size_t sz = sizeof(*g_m) * nx * ny;

  /* Allocate temporary ping-pong buffer */
  bool * t_m = new bool[sz];
  if (t_m == nullptr) {
    return 1;
  }

  bool *dest_m = t_m;
  bool *src_m = g_m;

  /* Ping-pong results between buffers */
  for (std::size_t i=0; i<ni; ++i) {
    gol_kernel(dest_m, src_m, nx, ny);
    dest_m = std::exchange(src_m, dest_m);
  }

  /* Odd iterations need to move back to input buffer */
  if (ni%2 != 0) {
    std::memcpy(g_m, t_m, sizeof(sz));
  }

  delete[] t_m;

  return 0; //TODO: Add dedicated return values
}



static constexpr bool gol_rule(int s, bool p)
{
  if (s == 3) {
    return true;
  }

  return p && s == 4;
}


static constexpr std::size_t cycle(std::size_t i, int di, std::size_t ni)
{
  return ((di % ni) + ni) % ni;
}


static void gol_kernel(bool *d_m, bool *s_m, std::size_t nx, std::size_t ny)
{
  for (std::size_t j=0; j<ny; ++j) {
    for (std::size_t i=0; i<nx; ++i) {
      int s = 0;

      for (int dj=-1; dj<1; ++dj) {
        for (int di=-1; di<1; ++di) {
          std::size_t idx = cycle(j, dj, ny) * nx + cycle(i, di, nx);
          s += s_m[idx];
        }
      }

      std::size_t idx = j*nx+i;

      d_m[idx] = gol_rule(s, s_m[idx]);
    }
  }
}

