import cffi

ffibuilder = cffi.FFI()
ffibuilder.cdef(
    # """int calculate_sf(double **, double **, double **,
    #                                 int *, int, int*, int,
    #                                 int**, double **, int,
    #                                 double**, double**, double**);"""
    """int calculate_atomistic_mcsh(double **, double **, double **,
                                    int *, int, int*, int,
                                    int**, double **, int, double **, int*, 
                                    double**, double**);
    """
)
ffibuilder.set_source(
    "simple_nn.features.MCSH._libmcsh",
    '#include "calculate_atomistic_mcsh.h"',
    sources=[
        "simple_nn/features/MCSH/calculate_atomistic_mcsh.cpp",
        "simple_nn/features/MCSH/atomistic_mcsh.cpp"
    ],
    source_extension=".cpp",
    include_dirs=["simple_nn/features/MCSH/"],
)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
