#include "defs.hpp"
#include "gol_seq.hpp"
#include "sfcurves.hpp"
#include "wrap_strategies.hpp"

#include <cstdio>
#include <cstring>

static void print2d(bool const *buf, int nx, int ny);

int main(void)
{
  GOL::Cell sbuff[] = {
    0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0
  };

  GOL::Cell cbuff[sizeof(sbuff)];

  std::memcpy(cbuff, sbuff, sizeof(sbuff));

  printf("Before:\n");
  print2d(cbuff, 7, 7);
  GOL::SEQ::gol_seq<GOL::SFC::Curve::ROW_MAJOR, GOL::WS::Strategy::GENERIC, GOL::WS::Strategy::GENERIC>(cbuff, 7, 7, 2);
  printf("After\n");
  print2d(cbuff, 7, 7);
  return std::memcmp(cbuff, sbuff, sizeof(sbuff));
}


static void print2d(bool const *buf, int nx, int ny)
{
  for (int j=0; j<ny; ++j) {
    for (int i=0; i<nx; ++i) {
      printf("%1d ", buf[i+nx*j]);
    }
    printf("\n");
  }
}

