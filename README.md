Code Information
================

This code contains twelve *functions*, two *main functions*, a *header file* and a *makefile* for the least-squares fitting of an ellipsoid to a large set of points. The names of the executable programs for the two main techniques are:

* **separation_in_groups**
* **sequential_adjustments**

The input files must include the Cartesian coordinates and the weight of every point, specifically (x, y, z, w). Next, we will provide an example of the text file format (Cartesian coordinates and weights).

```
7.0 22.0 31.0 1.0
7.0 19.0 28.0 1.0
9.0 23.0 31.0 1.0
```

The numbers must be separated by spaces and the text files are created, they must be converted into binary files so that they can be used by the program.

## Makefile instructions

To run the code, it is necessary to determine the path of the ellipsoid functions in the *makefile*. In the first line of the *makefile*, the path can be entered as a value in the **IDIR** parameter.

Example:

```makefile
IDIR = /home/myname/ellipsoid_functions
```

After determining the path, it is necessary to compile all programs by typing the following command in the Linux command line:

```bash
make all
```

By running this command, the two executable programs described above will be created.

Executing the code
------------------

To run the code **separation_in_groups**, type the following command in the Linux command line:

```bash
./separation_in_groups group1.bin group2.bin
```

If there are numerous files to be included and their names are numbered, the **xargs** command can be used as follows:

```bash
find /measurements_path -type f -name "filename*.bin" | xargs ./separation_in_groups
```

For example, include 40 files named group1.bin, group2.bin, ..., group40.bin, from the file "ellipsoid_points":

```bash
find /home/myname/ellipsoid_points -type f -name "group*.bin" | xargs ./separation_in_groups
```

---

To run the code **sequential_adjustments**, type the following command in the Linux command line:

```bash
./sequential_adjustments group1.bin group2.bin
```

*For this program to function properly, it is necessary to have at least two measurement files.

If there are numerous files to be included and their names are numbered, the **xargs** command can be used as follows:

```bash
find /measurements_path -type f -name "filename*.bin" | xargs ./sequential_adjustments
```

For example, include 40 files named group1.bin, group2.bin, ..., group40.bin, from the file "ellipsoid_points":

```bash
find /home/myname/ellipsoid_points -type f -name "group*.bin" | xargs ./sequential_adjustments
```

## Cleaning the code

To clean all the **.o** files (which are typically kept to avoid recompiling unchanged source files) and the executables, type the following command:

```bash
make clean
```

After cleaning the files, you can re-compile the entire program by typing the command again:

```bash
make all
```
