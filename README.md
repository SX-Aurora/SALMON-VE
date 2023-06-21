# SALMON-VE: Scalable Ab-initio Light-Matter simulator for Optics and Nanoscience optimized for SX-Aurora TSUBASA Vector Engine.

SALMON (http://salmon-tddft.jp/) is an open-source software based on first-principles time-dependent density functional theory
to describe optical responses and electron dynamics in matters induced by light electromagnetic fields.
This is the optimized version for SX-Aurora TSUBASA Vector Engine base on SALMON-v.2.1.0.

## How to make

Make the executable "salmon" and shared object "libvhcall.so" as follows. NOTE that put "libvhcall.so" into the run directory.

```
cd SALMON-VE/gnu_makefiles
make -f Makefile.sx # Executable "bin/salmon" is made.
cd ../libvhcall
make                # Shared object "libvhcall.so" is made.
```

## License

SALMON is available under Apache License version 2.0.

    Copyright 2017-2022 SALMON developers

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
