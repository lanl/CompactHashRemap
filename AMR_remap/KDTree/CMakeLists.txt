set(libkdtree_LIB_SRCS Bounds1d.c Bounds2d.c KDTree1d.c KDTree2d.c
    Globals2d.h  KDTree1d.h  KDTree2d.h Bounds1d.h  Bounds2d.h  Globals1d.h)


add_library(kdtree STATIC ${libkdtree_LIB_SRCS})

