mkdir run$1 
g++ -lm main.cpp -O3 -o solve
cp solve run$1/. 
cd run$1
./solve &
cd .. 
