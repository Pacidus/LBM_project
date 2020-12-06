MAKEFLAGS += --no-print-directory

make: windtunnel.f90
	gfortran windtunnel.f90 -O3 -Wall -o prog
	make clean

Rmake: Rewindtunnel.f90
	gfortran Rewindtunnel.f90 -O3 -Wall -o prog
	make clean

lmake: windtunnel.f90
	gfortran windtunnel.f90 -O3 -mcmodel=large -o prog
	make clean

movie: ./anim/*.png
	@echo y | ffmpeg -r 25 -i "./anim/%05d.png" -loop 1 -i "./objet.pbm" \
	-filter_complex "color=white:s=10x10[th];\
	[th][0]scale2ref=w=iw:h=ih[th][v];[1][th][1][v]threshold" -qscale:v 0 \
	-c:v libvpx-vp9 -crf 0 rot.mp4

out: ./anim2/*.svg
	 @echo y | ffmpeg -r 30 -i "./anim2/%05d.svg" -qscale:v 0 -c:v libvpx-vp9 \
	 -vf scale=1080:1080 -crf 20 out.mp4

clean:
	rm *.mod

cleandta:
	rm ./P/*.csv
	rm ./Ux/*.csv
	rm ./Uy/*.csv
