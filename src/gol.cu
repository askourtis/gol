#include <utility>
#include <cstring>
#include <cstdlib>

#ifndef CUDA_NB
#define CUDA_NB 1
#endif

#ifndef CUDA_NT
#define CUDA_NT
#endif

#ifdef CUDA_RUN
#define CUDA_PARAMS <<<CUDA_NB, CUDA_NT>>>
#define HOST __host__
#define DEVICE __device__
#define THREAD_IDX() threadIdx
#define BLOCK_IDX()  blockIdx
#define BLOCK_DIM()  blockDim
#else
#define CUDA_PARAMS
#define HOST
#define DEVICE
#define THREAD_IDX() []() {struct {Size x; Size y; Size z;} oobj{1,1,1}; return oobj; }()
#define BLOCK_IDX()  []() {struct {Size x; Size y; Size z;} oobj{1,1,1}; return oobj; }()
#define BLOCK_DIM()  []() {struct {Size x; Size y; Size z;} oobj{1,1,1}; return oobj; }()
#endif




namespace GOL {
  using Index = int;
  using Size = int;

  enum class Cell : int {
    DEAD = 0,
    ALIVE = 1
  };

  enum class WrapBehaviour {
    NOWRAP,
    CUTOFF,
    GENERIC,
    POW2,
  };

  template<WrapBehaviour WB>
  Index wrapIndex(Index x, Size nx);

  template<>
  Index wrapIndex<WrapBehaviour::NOWRAP>(Index x, Size nx) {
    return x;
  }

  template<>
  Index wrapIndex<WrapBehaviour::CUTOFF>(Index x, Size nx) {
    if (x < 0 || x >= nx) {
      return -1;
    }
    return x;
  }

  template<>
  Index wrapIndex<WrapBehaviour::GENERIC>(Index x, Size nx) {
    if (x < 0 || x >= nx) {
      return ((x%nx)+nx)%nx;
    }
    return x;
  }

  template<>
  Index wrapIndex<WrapBehaviour::POW2>(Index x, Size nx) {
    return x & (nx-1);
  }


  enum class SpaceFillingCurve {
    ROW_MAJOR,
    COL_MAJOR,
    Z_ORDER,
  };

  template<SpaceFillingCurve SFC>
  Index reduceIndexDimensions(Index x, Index y, Size nx, Size ny);


  template<>
  Index reduceIndexDimensions<SpaceFillingCurve::ROW_MAJOR>(Index x, Index y, Size nx, Size ny) {
    return x + y * nx;
  }

  template<>
  Index reduceIndexDimensions<SpaceFillingCurve::COL_MAJOR>(Index x, Index y, Size nx, Size ny) {
    return reduceIndexDimensions<SpaceFillingCurve::ROW_MAJOR>(y, x, ny, nx);
  }


  template<SpaceFillingCurve SFC, WrapBehaviour WB_X, WrapBehaviour WB_Y>
  Index getCellIndex(Index x, Index y, Size nx, Size ny) {
    auto rx = wrapIndex<WB_X>(x, nx);
    if (rx < 0) {
      return -1;
    }

    auto ry = wrapIndex<WB_Y>(y, ny);
    if (ry < 0) {
      return -1;
    }
    return reduceIndexDimensions<SFC>(rx, ry, nx, ny);
  }

  template<SpaceFillingCurve SFC, WrapBehaviour WB_X, WrapBehaviour WB_Y>
  Cell getCell(Cell *board, Index x, Index y, Size nx, Size ny) {
    auto idx = getCellIndex<SFC, WB_X, WB_Y>(x, y, nx, ny);
    if (idx < 0) {
      return Cell::DEAD;
    }

    return board[idx];
  }


  template<SpaceFillingCurve SFC, WrapBehaviour WB_X, WrapBehaviour WB_Y>
  DEVICE HOST Size gol_stencil(Cell *board, Index x, Index y, Size nx, Size ny) {
    Size cnt = 0;
    for (Index yy=y-1; yy<=+1; ++yy) {
      for (Index xx=x-1; xx<=+1; ++xx) {
        if (xx == x && yy == y) {
          continue;
        }
        cnt += (Size)getCell<SFC, WB_X, WB_Y>(board, xx, yy, nx, ny);
      }
    }
    return cnt;
  }


  Cell gol_rules(Cell state, Size alive_count)
  {
    return Cell::ALIVE;
  }


  template<SpaceFillingCurve SFC, WrapBehaviour WB_X, WrapBehaviour WB_Y>
  void gol_kernel(Cell *dstBoard, Cell *srcBoard, Size nx, Size ny) {
    for (Index y=THREAD_IDX().y; y<ny; y+=BLOCK_DIM().y) {
      for (Index x=THREAD_IDX().x; x<nx; x+=BLOCK_DIM().x) {
        Size alive_count = gol_stencil<SFC, WB_X, WB_Y>(srcBoard, x, y, nx, ny);
        auto idx = getCellIndex<SFC, WB_X, WB_Y>(x, y, nx, ny);
        dstBoard[idx] = gol_rules(dstBoard[idx], alive_count);
      }
    }
  }

  template<SpaceFillingCurve SFC, WrapBehaviour WB_X, WrapBehaviour WB_Y>
  void gol_impl(Cell *dstBoard, Cell *srcBoard, Size nx, Size ny, Size k) {
    Size boardByteCount = nx*ny*sizeof(Cell);

    for (Size i=0; i<k; ++i) {
      gol_kernel<SFC, WB_X, WB_Y> CUDA_PARAMS (dstBoard, srcBoard, nx, ny);
      dstBoard = std::exchange(srcBoard, dstBoard);
    }

    if (k%2 == 1) {
      std::memcpy(dstBoard, srcBoard, boardByteCount);
    }
  }


  void *allocate(Size sz) {
    return malloc(sz);
  }

  void release(void *ptr) {
    free(ptr);
  }


  template<SpaceFillingCurve SFC, WrapBehaviour WB_X, WrapBehaviour WB_Y>
  void gol(Cell *board, Size nx, Size ny, Size k) {
    Size boardByteCount = nx*ny*sizeof(Cell);
    Cell *tboard = (Cell*)allocate(boardByteCount);

    gol_impl<SFC, WB_X, WB_Y>(tboard, board, nx, ny, k);

    release(tboard);
  }

}


int main() {
  GOL::Cell board[25] = {
    GOL::Cell::DEAD
  };
  GOL::gol<GOL::SpaceFillingCurve::COL_MAJOR, GOL::WrapBehaviour::GENERIC, GOL::WrapBehaviour::GENERIC> (board, 5, 5, 5);
}


