gcc -c -fPIC -shared vl/sift.c -o libsift.so
g++ -c example.cpp -o example.o
# add opencv flag
g++ example.o -o example -L. -lsift `pkg-config --cflags --libs opencv`
# remove internal files
rm -f example.o libsift.so
