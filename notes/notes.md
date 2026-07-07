Right now, we're sitting at a ~15 minute computation time for the Barnes-Hut
computation. My goal is to hopefully decrease this time down to ~10 seconds
(about an order of magnitude or two decrease).

**Big Goals**

1. Vectorize the Octree Implementation
2. Parallelize the Barnes-Hut computation
3. Generalize the computation for unique fish
4. Implement the fast multipole method
5. Implement in GPU computation with CUDA

From the "Many-Body Tree Methods in Physics" book, there are a TON of
optimizations I can make. My main goal is to try and implement in as many as
possible, such as vectorizing and parallelizing the Barnes-Hut implementation,
implementing in the Fast Multipole Method, and then also using CUDA for GPU
computing, as described by the two below papers:

### Reading List

1. Many-Body Tree Methods in Physics (Chapters 2, 4, and 7) by Gibbon and
   Pfalzner
1. An Efficient CUDA Implementation of the Tree-Based Barnes Hut n-Body
   Algorithm[^1] by Burtscher and Pingali
2. Introduction to Fast Multipole Methods[^2] by Chen
3. Real-time N-body Simulation with Barnes-Hut Algorithm and CUDA [^3] by Wu
4. A short course on fast multipole methods [^4] by Beatson and Greengard

[^1]: [ScienceDirect](https://www.sciencedirect.com/science/chapter/edited-volume/pii/B9780123849885000061)
[^2]: [UC Irvine](https://www.math.uci.edu/~chenlong/226/FMMsimple.pdf)
[^3]: [GitHub](https://hsin-hung.github.io/N-body-simulation/)
[^4]: [New York University](https://math.nyu.edu/~greengar/shortcourse_fmm.pdf)
