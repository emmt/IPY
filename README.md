IPY: Inverse Problems with Yorick
=================================

**IPY** is a package of tools to help solving inverse problems with Yorick.


Installation
------------

In short, building and installing the plug-in can be as quick as:
```
    cd BUILD_DIR
    SRC_DIR/configure
    make
    make install
```
where `BUILD_DIR` is the build directory (at your convenience) and `SRC_DIR` is
the source directory of the plug-in code.  The build and source directories
can be the same in which case, call `./configure` to configure for building.

Then, to use the plug-in, start the Yorick interpreter and type:
```
    #include "ipy.i"
```
More detailled explanations are given below.

0. You must have [Yorick](http://github.com/LLNL/yorick/) and
   [Yeti](https://github.com/emmt/Yeti) installed on your machine.

1. Unpack the plug-in code somewhere.

2. Configure for compilation.  The are two possibilities:

   * For an in-place build, go to the source directory of the plug-in code
     and run the configuration script:
     ```
         cd SRC_DIR
         ./configure
     ```
     To see the configuration options, call:
     ```
         ./configure --help
     ```

   * To compile in a different build directory, say `BUILD_DIR`, create the
     build directory, go to the build directory, and run the configuration
     script:
     ```
         mkdir -p BUILD_DIR
         cd BUILD_DIR
         SRC_DIR/configure
     ```
     where `SRC_DIR` is the path to the source directory of the plug-in code.
     To see the configuration options, call:
     ```
         SRC_DIR/configure --help
     ```

3. Compile the code:
   ```
       make
   ```

4. Install the plug-in in Yorick directories:
   ```
       make install
   ```

License
-------

IPY is open source sofware released under the GNU-GPL License
(see [./LICENSE.md](./LICENSE.md)).

