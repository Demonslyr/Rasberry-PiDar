

make: 
	@echo ""
	@echo "----Makefile Description----"
	@echo "post: for the post processing version"
	@echo "no_post: for the non post processed version"
	@echo "info: for instructions"
	@echo "all: for all versions"
	@echo "clean: to erase executables"
	@echo ""

all: post no_post

no_post:
	g++ main_no_post.cpp -o Fourier_2D_no_post -lfftw3 -lrt -ffast-math -std=c++11 -O3 -march=native

post:
	g++ main_post.cpp -o Fourier_2D_post -lfftw3 -lrt -ffast-math -std=c++11 -O3 -march=native

info:
	@echo ""
	@echo "$$ make (no_)post {PC/Pi2/Cust} {filename} {header lines}"
	@echo ""
	@echo "{filename} is required only for {Cust}"
	@echo "only input the name of the file without the file extension"
	@echo "This must be a file in the Test_Data folder and must be a .csv file"
	@echo ""
	@echo "{header_lines} is the number of lines at the top of te .csv to ignore"
	@echo ""

clean:
	-rm Fourier_2D_no_post
	-rm Fourier_2D_post
