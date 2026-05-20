
#include <iomanip>
#include <iostream>
#include <vector>

template <typename T>
struct Array3D
{
    T* __restrict__ p = nullptr;
    std::size_t jstride = 0;
    std::size_t kstride = 0;
    int ibegin = 0;
    int jbegin = 0;
    int kbegin = 0;
    int iend = -1;
    int jend = -1;
    int kend = -1;

    Array3D (T* a_p, int a_ibegin, int a_jbegin, int a_kbegin,
             int a_iend, int a_jend, int a_kend)
        : p(a_p),
          jstride(a_iend-a_ibegin),
          kstride(jstride*(a_jend-a_jbegin)),
          ibegin(a_ibegin),
          jbegin(a_jbegin),
          kbegin(a_kbegin),
          iend(a_iend),
          jend(a_jend),
          kend(a_kend)
    {}

    inline __attribute__((always_inline))
    T& operator() (int i, int j, int k) const
    {
        return p[(i-ibegin)+(j-jbegin)*jstride+(k-kbegin)*kstride];
    }
};

void add (double* a, double b) { *a += b; }

void smoother (Array3D<double> const& phi, Array3D<double> const& rhs,
               Array3D<double> const& sig, double dx, double dy, double dz)
{
    double facx = (1.0/36.0) / (dx*dx);
    double facy = (1.0/36.0) / (dy*dy);
    double facz = (1.0/36.0) / (dz*dz);
    double fxyz = facx + facy + facz;
    double fmx2y2z   =     -facx + 2.0*facy + 2.0*facz;
    double f2xmy2z   =  2.0*facx -     facy + 2.0*facz;
    double f2x2ymz   =  2.0*facx + 2.0*facy -     facz;
    double f4xm2ym2z =  4.0*facx - 2.0*facy - 2.0*facz;
    double fm2x4ym2z = -2.0*facx + 4.0*facy - 2.0*facz;
    double fm2xm2y4z = -2.0*facx - 2.0*facy + 4.0*facz;

#pragma omp parallel for collapse(2)
    for (int k = 0; k < rhs.kend; ++k) {
    for (int j = 0; j < rhs.jend; ++j) {
    for (int i = 0; i < rhs.iend; ++i) {
        double s0 = -4.0*fxyz*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1)
                              +sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  ));
        double Ax = phi(i,j,k)*s0
            + fxyz*(phi(i-1,j-1,k-1)*sig(i-1,j-1,k-1)
                  + phi(i+1,j-1,k-1)*sig(i  ,j-1,k-1)
                  + phi(i-1,j+1,k-1)*sig(i-1,j  ,k-1)
                  + phi(i+1,j+1,k-1)*sig(i  ,j  ,k-1)
                  + phi(i-1,j-1,k+1)*sig(i-1,j-1,k  )
                  + phi(i+1,j-1,k+1)*sig(i  ,j-1,k  )
                  + phi(i-1,j+1,k+1)*sig(i-1,j  ,k  )
                  + phi(i+1,j+1,k+1)*sig(i  ,j  ,k  ))
            + fmx2y2z*(phi(i  ,j-1,k-1)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1))
                     + phi(i  ,j+1,k-1)*(sig(i-1,j  ,k-1)+sig(i,j  ,k-1))
                     + phi(i  ,j-1,k+1)*(sig(i-1,j-1,k  )+sig(i,j-1,k  ))
                     + phi(i  ,j+1,k+1)*(sig(i-1,j  ,k  )+sig(i,j  ,k  )))
            + f2xmy2z*(phi(i-1,j  ,k-1)*(sig(i-1,j-1,k-1)+sig(i-1,j,k-1))
                     + phi(i+1,j  ,k-1)*(sig(i  ,j-1,k-1)+sig(i  ,j,k-1))
                     + phi(i-1,j  ,k+1)*(sig(i-1,j-1,k  )+sig(i-1,j,k  ))
                     + phi(i+1,j  ,k+1)*(sig(i  ,j-1,k  )+sig(i  ,j,k  )))
            + f2x2ymz*(phi(i-1,j-1,k  )*(sig(i-1,j-1,k-1)+sig(i-1,j-1,k))
                     + phi(i+1,j-1,k  )*(sig(i  ,j-1,k-1)+sig(i  ,j-1,k))
                     + phi(i-1,j+1,k  )*(sig(i-1,j  ,k-1)+sig(i-1,j  ,k))
                     + phi(i+1,j+1,k  )*(sig(i  ,j  ,k-1)+sig(i  ,j  ,k)))
            + f4xm2ym2z*(phi(i-1,j,k)*(sig(i-1,j-1,k-1)+sig(i-1,j,k-1)+sig(i-1,j-1,k)+sig(i-1,j,k))
                       + phi(i+1,j,k)*(sig(i  ,j-1,k-1)+sig(i  ,j,k-1)+sig(i  ,j-1,k)+sig(i  ,j,k)))
            + fm2x4ym2z*(phi(i,j-1,k)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j-1,k)+sig(i,j-1,k))
                       + phi(i,j+1,k)*(sig(i-1,j  ,k-1)+sig(i,j  ,k-1)+sig(i-1,j  ,k)+sig(i,j  ,k)))
            + fm2xm2y4z*(phi(i,j,k-1)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1))
                       + phi(i,j,k+1)*(sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  )));

        phi(i,j,k) += (rhs(i,j,k) - Ax) / s0;
    }}}
}

template <typename fty> fty *__raptor_truncate_op_func(fty *, int, int, int);

int main (int argc, char* argv[])
{
    int nx = 128;
    int ny = 128;
    int nz = 128;

    double dx = 1.0/nx;
    double dy = 1.0/ny;
    double dz = 1.0/nz;

    // 3d array on nodal grid with one ghost cell
    std::vector<double> phi(std::size_t(nx+3)*(ny+3)*(nz+3), 0.0);
    Array3D<double> phi_a{phi.data(), -1, -1, -1, nx+1, ny+1, nz+1};

    // 3d array on nodal grid without ghost cells
    std::vector<double> rhs(std::size_t(nx+1)*(ny+1)*(nz+1), 0.5);
    Array3D<double> rhs_a{rhs.data(), 0, 0, 0, nx, ny, nz};

    // 3d array on cell grid with one ghost cell
    std::vector<double> sig(std::size_t(nx+2)*(ny+2)*(nz+2), 0.8);
    Array3D<double> sig_a{sig.data(), -1, -1, -1, nx, ny, nz};

    auto smoother_wrapper = __raptor_truncate_op_func(smoother, 64, 0, 32);
    smoother_wrapper(phi_a, rhs_a, sig_a, dx, dy, dz);

    std::cout << std::setprecision(17) << "dx = " << dx << '\n';

    auto add_wrapper = __raptor_truncate_op_func(add, 64, 0, 16);
    add_wrapper(&dx, 1.e-9);

    std::cout << "dx = " << dx << '\n';

    add(&dx, 1.e-9);

    std::cout << "dx = " << dx << '\n';

    return 0;
}
