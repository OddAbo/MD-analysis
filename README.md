# READ ME

Codes for analyzing Molecular Dynamics program results

2020-03-21

Author: OddAbo

> See [**reference**](https://github.com/OddAbo/MD-analysis/tree/master/reference) for more computing details

---

## [**getarg.c**](https://github.com/OddAbo/MD-analysis/blob/master/src/getarg.c)

A simple example shows how you can allocate memory to array depending on input
arguments.
```shellscript
./GetArg  natom 1024  nstep 1  dt 0.001
```
---

## [**readxyz.c**](https://github.com/OddAbo/MD-analysis/blob/master/src/readxyz.c)

A simple example shows how you can read binary file formed by Fortran (access =
    'stream') in C language.


