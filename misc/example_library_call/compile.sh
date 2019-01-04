cd ../../
./compile.sh
cd misc/example_library_call
cp ../../deploy/*.h . 
scons program=interfacetest variant=optimized -j 8
cp ./optimized/interface_test . 
rm -rf ./optimized/
rm config.log
