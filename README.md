##Compiling##
On Linux, compiling can be done by
```mkdir build
cd build
cmake ../
cmake --build .
```

In order to compile the test folders on a Linux system, the prerequisite `libgtest-dev` is needed. Move to the `src` folder and run the following:
```
g++ -c *.cpp
mv *.o ../tests
g++ -Wall -g -pthread *.o FILENAME -lgtest_main -lgtest -lpthread -o OUTPUT
```
where FILENAME, OUTPUT are replaced by the desired names.
