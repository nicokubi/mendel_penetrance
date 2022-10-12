# 2020.01.17.1107: Per Chris, installing mendel
# but have to build gcc with fortran etc. first

setup1st:
	@echo "Make sure you do the following first:"
	@echo "module load gcc/9.2.0"
	@echo "Next: make go"

go:clean
	gfortran -O3 mendela.f single.f -o single
	gfortran -O3 mendela.f poly.f -o poly

clean:
	\rm -rf single poly



