# attenuatingsphere
API for average weighted distance with applications such as updating an audio emitter along a path.

The code is organized so that it is self-contained. I did not spend much time trying to optimize. Parallelization or SIMD would be natural approaches to optimize such an algorithm. I used double precision as I had precision issues using float. The precision issues can be handled by clamping spread to between or equal to 0 and 1, checking for cases where the clipped line almost becomes a point, and checking for points very close to the sphere's center.

A full explanation of the math used to derive these functions can be found on my blog here: https://computingandrecording.wordpress.com/
