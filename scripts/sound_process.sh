#!/bin/bash 

#-----------------------------------------------------------
# This script converts NIST SPHERE .sph files to .wav
# It normalizes the amplitude to -1dB, applies a lowpass
# at 3400Hz and resamples the output to 8000Hz
#
# This can is used to produce a training set
# for LSF vector quantization using TED-LIUM speech corpus.
#
# SP5WWP 07/04/2020
#-----------------------------------------------------------

#convert all .sph files to .wav
#do it if the .wav file doesn't exist already
#trim out the crap (20s on both ends)
#of course some windowing should be added
#to avoid sudden signal changes,
#but we can live with that
for f in *.sph
do
	if [ ! -f "${f%.*}.wav" ]
	then
		sox -t sph "$f" -b 16 -t wav "${f%.*}.wav" trim 20 -20
	fi
done

#normalize to -1dB
for f in *.wav
do
	sox --norm=-1 "$f" "${f%.*}_new.wav"
	mv "${f%.*}_new.wav" "$f"
done

#lowpass at 3.4kHz, 4 passes
#maybe sox would be better here
for f in *.wav
do
	for passes in {1..4}
	do
		ffmpeg -i "$f" -filter lowpass=3400 "${f%.*}_new.wav"
		mv "${f%.*}_new.wav" "$f"
	done
done

#resample to 8000Hz
for f in *.wav
do
	sox "$f" -r 8000 "${f%.*}_new.wav"
	mv "${f%.*}_new.wav" "$f"
done

#merge all pieces into one
sox $(ls *.wav) corpus.wav