clear
rm *.o
make
gdb --args ./stone_slide_main 10 
view ../output/console_output.txt
