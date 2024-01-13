ShapedSlp:
===============
Author: Tomohiro I

### Download

To download all the necessary source codes:
```sh
git clone --recursive https://github.com/itomomoti/ShapedSlp.git
```

### Compile

```sh
mkdir build
cd build
cmake ..
make
```

### Usage

Executables for benchmarks. See help by running without options.

```sh
./SlpEncBuild
./SubstrBenchmark
```

### Example

```sh
# build encoding for method SelfShapedSlp_SdSd_Sd
./SlpEncBuild -i path_to_data/chr19x1000.fa -o ds_fname -e SelfShapedSlp_SdSd_Sd -f Bigrepair
# With "-e All", build encodings for all methods. The file name for each method is prefixed by the string given by -o option and suffixed by the method name
./SlpEncBuild -i path_to_data/chr19x1000.fa -o chr_ -e All -f Bigrepair
```

```sh
# substring
./SubstrBenchmark -i ds_fname -e SelfShapedSlp_SdSd_Sd -n 10000 -l 10 -j 1000000
# With "-e All", test all methods. The base name should be given by -i option.
./SubstrBenchmark -i chr_ -e All -n 10000 -l 10 --dummy_flag 1 -j 1000000
```

### This branch
Here we implemented the query LCP(T[i..], P[j..]) [1, 2, 3] via Karp-Rabin hashes.

TODO:
+ Polish tests (make them independent of the large files that won't fit to GitHub)
+ Get rid of the slow modular inverse

#### References
[1] Baláz, A., Gagie, T., Goga, A., Heumos, S., Navarro, G., Petescia, A. & Sirén, J. (2023). Wheeler maps. arXiv preprint arXiv:2308.09836.

[2] Ahmed, O., Baláž, A., Brown, N. K., Depuydt, L., Goga, A., Petescia, A., Zakeri, M., Fostier, J., Gagie, T., Langmead, B., Navarro, G. & Prezza, N. (2023). r-indexing without backward searching. arXiv preprint arXiv:2312.01359.

[3] Goga, A., Depuydt, L., Brown, N. K., Fostier, J., Gagie, T., & Navarro, G. (2023). Faster Maximal Exact Matches with Lazy LCP Evaluation. arXiv preprint arXiv:2311.04538.
